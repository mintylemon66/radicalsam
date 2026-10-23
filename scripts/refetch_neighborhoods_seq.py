#!/usr/bin/env python3
"""
Refetch genomic neighborhoods WITH protein sequences and coordinates.

Why: the original pipeline stored only GenBank product-name TEXT, so every
DUF4229 call is a string match. This rebuilds neighborhoods carrying
(protein_id, product, coordinates, strand, translation) so that family
membership can be recalled from SEQUENCE via an HMM, and so that the
MqnE->DUF4229 genomic DISTANCE can be measured.

Window is deliberately WIDER (+/-25 kb, ~ +/-20 genes) than the original
+/-8 kb, to test whether DUF4229 sits just outside the original window in
genomes currently scored as negative.

Resumable: re-running skips accessions already present in the output.
Rate-limited to stay within NCBI's 10 req/s API-key allowance.
"""
import csv
import os
import re
import sys
import time
import threading
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "mqn_cooccur_FULL.tsv")
OUTDIR = os.path.join(ROOT, "audit")
OUT_TSV = os.path.join(OUTDIR, "neighbors_seq.tsv")
OUT_FAA = os.path.join(OUTDIR, "neighbors_seq.faa")
LOG = os.path.join(OUTDIR, "refetch.log")

WINDOW = 25000          # bp each side of the query CDS
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
WORKERS = 6
RATE = 9.0              # requests/sec ceiling (API key allows 10)

_keyfile = os.path.expanduser("~/.ncbi_api_key")
API_KEY = open(_keyfile).read().strip() if os.path.exists(_keyfile) else ""

_lock = threading.Lock()
_last = [0.0]


def throttle():
    with _lock:
        gap = 1.0 / RATE
        now = time.time()
        wait = _last[0] + gap - now
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()


def get(url, tries=4):
    for i in range(tries):
        throttle()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "radicalSAM-audit/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if i == tries - 1:
                return ""
            time.sleep(1.5 * (i + 1))
    return ""


def log(msg):
    with _lock:
        with open(LOG, "a") as fh:
            fh.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


# ---------- stage 1: UniProt accession -> (contig, embl_protein_id) ----------
def uniprot_batch(accs):
    """Batched UniProt lookup returning {acc: (contig, protein_id)}.

    The TSV projection of xref_embl omits the ProteinId, so we parse JSON
    cross-references, where the EMBL entry carries id=<contig> plus a
    ProteinId property = the EMBL/GenBank protein accession we need for
    the /coded_by lookup.
    """
    import json
    out = {}
    q = " OR ".join(f"accession:{a}" for a in accs)
    url = (UNIPROT + "?query=" + urllib.parse.quote(q)
           + "&fields=accession,xref_embl&format=json&size=500")
    txt = get(url)
    if not txt:
        return out
    try:
        data = json.loads(txt)
    except Exception:
        return out
    for entry in data.get("results", []):
        acc = entry.get("primaryAccession")
        if not acc:
            continue
        best = None
        for xr in entry.get("uniProtKBCrossReferences", []):
            if xr.get("database") != "EMBL":
                continue
            contig = xr.get("id")
            props = {p.get("key"): p.get("value") for p in xr.get("properties", [])}
            pid = props.get("ProteinId")
            if contig and pid and pid != "-":
                cand = (contig.strip(), pid.strip())
                if props.get("MoleculeType") == "Genomic_DNA":
                    best = cand
                    break
                best = best or cand
        if best:
            out[acc] = best
    return out


# ---------- stage 2: EMBL protein -> /coded_by contig:start..end ----------
CODED = re.compile(r'/coded_by="?(?:complement\()?([A-Za-z0-9_.]+):(?:<)?(\d+)\.\.(?:>)?(\d+)')


def coded_by(protein_id):
    url = (f"{EUTILS}/efetch.fcgi?db=protein&id={urllib.parse.quote(protein_id)}"
           f"&rettype=gp&retmode=text")
    if API_KEY:
        url += f"&api_key={API_KEY}"
    txt = get(url)
    m = CODED.search(txt)
    if not m:
        return None
    comp = "complement(" in txt[max(0, m.start() - 30):m.start() + 30]
    return m.group(1), int(m.group(2)), int(m.group(3)), ("-" if comp else "+")


# ---------- stage 3: fetch region, parse CDS ----------
def fetch_region(contig, start, end):
    s = max(1, start - WINDOW)
    e = end + WINDOW
    url = (f"{EUTILS}/efetch.fcgi?db=nuccore&id={urllib.parse.quote(contig)}"
           f"&seq_start={s}&seq_stop={e}&rettype=gbwithparts&retmode=text")
    if API_KEY:
        url += f"&api_key={API_KEY}"
    return get(url), s


def parse_cds(gb_text, offset):
    """Yield dicts for each CDS feature in the GenBank flat file."""
    from io import StringIO
    try:
        from Bio import SeqIO
    except ImportError:
        return []
    recs = []
    try:
        for rec in SeqIO.parse(StringIO(gb_text), "genbank"):
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                q = feat.qualifiers
                recs.append({
                    "protein_id": (q.get("protein_id") or ["-"])[0],
                    "product": (q.get("product") or ["-"])[0].replace("\t", " "),
                    "start": int(feat.location.start) + offset,
                    "end": int(feat.location.end) + offset,
                    "strand": "+" if feat.location.strand == 1 else "-",
                    "translation": (q.get("translation") or [""])[0],
                })
    except Exception:
        return []
    return recs


def done_set():
    seen = set()
    if os.path.exists(OUT_TSV):
        with open(OUT_TSV) as fh:
            for line in fh:
                seen.add(line.split("\t", 1)[0])
    seen.discard("query_acc")
    return seen


COLS = ["query_acc", "group", "contig", "q_start", "q_end", "q_strand",
        "n_protein_id", "n_product", "n_start", "n_end", "n_strand",
        "bp_gap", "gene_index"]


def process(acc, group, contig, protein_id, wlock):
    cb = coded_by(protein_id)
    if not cb:
        log(f"NOCODED {acc} {protein_id}")
        return 0
    ctg, qs, qe, qstrand = cb
    gb, offset = fetch_region(ctg, qs, qe)
    if not gb:
        log(f"NOREGION {acc} {ctg}")
        return 0
    cds = parse_cds(gb, offset)
    if not cds:
        log(f"NOCDS {acc} {ctg}")
        return 0
    cds.sort(key=lambda c: c["start"])
    # locate the query CDS (overlaps qs..qe) to index genes relative to it
    qi = None
    for i, c in enumerate(cds):
        if c["start"] <= qs <= c["end"] or c["start"] <= qe <= c["end"]:
            qi = i
            break
    rows, faa = [], []
    for i, c in enumerate(cds):
        if c["end"] < qs:
            gap = qs - c["end"]
        elif c["start"] > qe:
            gap = c["start"] - qe
        else:
            gap = 0
        gi = (i - qi) if qi is not None else ""
        rows.append([acc, group, ctg, qs, qe, qstrand,
                     c["protein_id"], c["product"], c["start"], c["end"],
                     c["strand"], gap, gi])
        if c["translation"] and c["protein_id"] != "-":
            faa.append(f">{c['protein_id']} {acc}|{group}|gi={gi}|gap={gap}\n{c['translation']}\n")
    with wlock:
        with open(OUT_TSV, "a", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        if faa:
            with open(OUT_FAA, "a") as fh:
                fh.write("".join(faa))
    return len(rows)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    only = sys.argv[1] if len(sys.argv) > 1 else None   # 'bg', 'mqn', or a file path
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    # If arg is a path, treat it as an explicit accession allow-list.
    allow = None
    if only and os.path.exists(only):
        allow = {l.strip() for l in open(only) if l.strip()}
        print(f"[refetch] restricting to {len(allow)} accessions from {only}")
        only = None

    targets = []
    with open(SRC, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = (row.get("group") or "").strip()
            a = (row.get("accession") or "").strip()
            if g in ("mqn", "bg") and a:
                if only and g != only:
                    continue
                if allow is not None and a not in allow:
                    continue
                targets.append((a, g))
    seen = done_set()
    targets = [t for t in targets if t[0] not in seen]
    if limit:
        targets = targets[:limit]

    if not os.path.exists(OUT_TSV):
        with open(OUT_TSV, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerow(COLS)

    print(f"[refetch] group={only or 'all'} todo={len(targets)} window=+/-{WINDOW}bp "
          f"apikey={'yes' if API_KEY else 'no'}", flush=True)

    # stage 1: batched UniProt mapping
    mapping = {}
    B = 60
    for i in range(0, len(targets), B):
        chunk = [a for a, _ in targets[i:i + B]]
        mapping.update(uniprot_batch(chunk))
        if (i // B) % 5 == 0:
            print(f"  [uniprot] {min(i+B, len(targets))}/{len(targets)} mapped={len(mapping)}", flush=True)
    print(f"[refetch] uniprot mapping complete: {len(mapping)}/{len(targets)}", flush=True)

    wlock = threading.Lock()
    ok = rowsum = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {}
        for acc, grp in targets:
            if acc not in mapping:
                continue
            contig, pid = mapping[acc]
            futs[ex.submit(process, acc, grp, contig, pid, wlock)] = acc
        total = len(futs)
        for n, fut in enumerate(as_completed(futs), 1):
            try:
                r = fut.result()
            except Exception as e:
                r = 0
                log(f"ERR {futs[fut]} {e}")
            if r:
                ok += 1
                rowsum += r
            if n % 100 == 0 or n == total:
                print(f"  [region] {n}/{total} ok={ok} rows={rowsum}", flush=True)
    print(f"[refetch] DONE group={only or 'all'} genomes_ok={ok} neighbor_rows={rowsum}", flush=True)


if __name__ == "__main__":
    main()
