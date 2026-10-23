#!/usr/bin/env python3
"""
Audit 2: structural comparison of neighborhood proteins.

Rationale: an HMM finds family members that are recognisable at the SEQUENCE
level. Remote homologs whose sequence has diverged below the Pfam gathering
threshold are invisible to both the text search AND the HMM, but remain
recognisable by FOLD. So we:

  1. take the neighborhood proteins (optionally only the annotation-opaque
     ones, which are the ones most likely to hide a missed DUF4229),
  2. dereplicate them so we do not pay for near-identical sequences,
  3. map GenBank/EMBL protein IDs -> UniProt accessions,
  4. download the corresponding AlphaFold DB models,
  5. Foldseek-search them against reference DUF4229 structures.

A background protein that is structurally DUF4229-like but sequence-divergent
would be a genuine hidden false negative, and would matter for the headline.

Stages are separable so a partial run is still useful.
"""
import csv
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
STRUCT = os.path.join(AUD, "struct")
PDBDIR = os.path.join(STRUCT, "pdb")
TSV = os.path.join(AUD, "neighbors_seq.tsv")
FAA = os.path.join(AUD, "neighbors_seq.faa")

OPAQUE = re.compile(
    r"hypothetical|uncharacteri[sz]ed|unknown function|^duf\d+|"
    r"predicted protein|conserved protein|putative protein", re.I)

RATE = 8.0
_lock = threading.Lock()
_last = [0.0]


def throttle():
    with _lock:
        now = time.time()
        w = _last[0] + 1.0 / RATE - now
        if w > 0:
            time.sleep(w)
        _last[0] = time.time()


def get(url, tries=3, binary=False):
    for i in range(tries):
        throttle()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "radicalSAM-audit/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read() if binary else r.read().decode("utf-8", "replace")
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(1.0 * (i + 1))
    return None


# ---------------- stage 1: select + dereplicate ----------------
def stage_select(opaque_only=True, group=None):
    """Write the candidate protein FASTA subset, then CD-HIT dereplicate."""
    os.makedirs(STRUCT, exist_ok=True)
    keep = {}
    with open(TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if group and r["group"] != group:
                continue
            pid = r["n_protein_id"]
            if pid == "-":
                continue
            if opaque_only and not OPAQUE.search(r["n_product"] or ""):
                continue
            keep[pid] = r["n_product"]
    print(f"[select] candidate proteins: {len(keep)} (opaque_only={opaque_only}, group={group or 'all'})")

    sub = os.path.join(STRUCT, "candidates.faa")
    n = 0
    with open(FAA) as fh, open(sub, "w") as out:
        pid = None
        for line in fh:
            if line.startswith(">"):
                pid = line[1:].split()[0]
            elif pid in keep:
                out.write(f">{pid}\n{line.strip()}\n")
                n += 1
                pid = None
    print(f"[select] wrote {n} sequences -> {sub}")

    derep = os.path.join(STRUCT, "candidates_nr.faa")
    try:
        subprocess.run(["cd-hit", "-i", sub, "-o", derep, "-c", "0.90",
                        "-n", "5", "-M", "2000", "-d", "0"],
                       check=True, capture_output=True)
        m = sum(1 for l in open(derep) if l.startswith(">"))
        print(f"[derep] CD-HIT 90% -> {m} representatives ({derep})")
    except Exception as e:
        print(f"[derep] cd-hit failed ({e}); using undereplicated set")
        derep = sub
    return derep


# ---------------- stage 2: map protein IDs -> UniProt ----------------
def stage_map(fasta):
    ids = [l[1:].split()[0] for l in open(fasta) if l.startswith(">")]
    print(f"[map] mapping {len(ids)} protein IDs to UniProt")
    out = {}
    B = 200
    for i in range(0, len(ids), B):
        chunk = ids[i:i + B]
        q = " OR ".join(f'xref:embl-{c.split(".")[0]}' for c in chunk)
        url = ("https://rest.uniprot.org/uniprotkb/search?query="
               + urllib.parse.quote(q)
               + "&fields=accession,xref_embl&format=json&size=500")
        txt = get(url)
        if not txt:
            continue
        try:
            data = json.loads(txt)
        except Exception:
            continue
        for e in data.get("results", []):
            acc = e.get("primaryAccession")
            for xr in e.get("uniProtKBCrossReferences", []):
                if xr.get("database") != "EMBL":
                    continue
                props = {p.get("key"): p.get("value") for p in xr.get("properties", [])}
                pid = props.get("ProteinId")
                if pid:
                    out[pid] = acc
                    out[pid.split(".")[0]] = acc
        print(f"  [map] {min(i+B, len(ids))}/{len(ids)} mapped={len(out)}", flush=True)
    mp = os.path.join(STRUCT, "protein_to_uniprot.tsv")
    with open(mp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "uniprot"])
        for k, v in out.items():
            w.writerow([k, v])
    print(f"[map] wrote {mp}")
    return out


# ---------------- stage 3: fetch AlphaFold DB models ----------------
def afdb_one(uni):
    dest = os.path.join(PDBDIR, f"{uni}.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return True
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uni}-F1-model_v4.pdb"
    data = get(url, binary=True)
    if not data or not data.startswith(b"ATOM") and b"ATOM" not in data[:2000]:
        return False
    with open(dest, "wb") as fh:
        fh.write(data)
    return True


def stage_fetch(mapping, limit=0):
    os.makedirs(PDBDIR, exist_ok=True)
    unis = sorted(set(mapping.values()))
    if limit:
        unis = unis[:limit]
    print(f"[afdb] fetching {len(unis)} AlphaFold models")
    ok = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(afdb_one, u): u for u in unis}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                if f.result():
                    ok += 1
            except Exception:
                pass
            if n % 100 == 0 or n == len(futs):
                print(f"  [afdb] {n}/{len(futs)} downloaded={ok}", flush=True)
    print(f"[afdb] done: {ok} structures in {PDBDIR}")
    return ok


# ---------------- stage 4: foldseek ----------------
def stage_foldseek(query_pdb):
    if not os.path.exists(query_pdb):
        print(f"[foldseek] query structure missing: {query_pdb}")
        return
    res = os.path.join(STRUCT, "foldseek_hits.tsv")
    tmp = os.path.join(STRUCT, "tmp")
    cmd = ["foldseek", "easy-search", query_pdb, PDBDIR, res, tmp,
           "--format-output", "query,target,fident,alnlen,evalue,bits,alntmscore",
           "-e", "10", "--max-seqs", "2000"]
    print("[foldseek]", " ".join(cmd[:4]), "...")
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print("[foldseek] failed:", e.stderr.decode()[:400])
        return
    rows = [l.split("\t") for l in open(res)] if os.path.exists(res) else []
    print(f"[foldseek] {len(rows)} hits -> {res}")
    for r in rows[:25]:
        print("   ", "\t".join(r).strip())


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("select", "all"):
        derep = stage_select(opaque_only=True)
    else:
        derep = os.path.join(STRUCT, "candidates_nr.faa")
    if stage in ("map", "all"):
        mapping = stage_map(derep)
    else:
        mapping = {}
        p = os.path.join(STRUCT, "protein_to_uniprot.tsv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                mapping[r["protein_id"]] = r["uniprot"]
    if stage in ("fetch", "all"):
        stage_fetch(mapping)
    if stage in ("foldseek", "all"):
        q = os.path.join(ROOT, "ml", "DUF4229_esmfold_atlas.pdb")
        stage_foldseek(q)


if __name__ == "__main__":
    main()
