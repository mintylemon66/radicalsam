#!/usr/bin/env python3
"""
Full structural sweep of the background, done tractably.

Question: could a background protein be a structural DUF4229 mimic that BOTH
the text search and the HMM missed? Remote homologs can retain fold while
losing sequence identity.

Fetching AlphaFold models for all 13,148 dereplicated opaque background
proteins is not feasible. But a structural DUF4229 homolog must share DUF4229's
physical profile: a small protein (88 to 128 aa) with roughly 2 to 3
transmembrane helices. An 800-residue soluble enzyme cannot be a DUF4229
mimic. So we filter on that profile first, which is a motivated biophysical
prefilter rather than an arbitrary cut, then fetch structures only for the
survivors and Foldseek them against the DUF4229 reference.

Stages: filter -> map to UniProt -> fetch AFDB (v6 API) -> foldseek.
"""
import csv
import json
import os
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ST = os.path.join(ROOT, "audit", "struct")
PDB = os.path.join(ST, "sweep_pdb")
SRC = os.path.join(ST, "candidates_nr.faa")
REF = os.path.join(ROOT, "ml", "DUF4229_esmfold_atlas.pdb")

# DUF4229 profile: 88-128 aa, mean ~2.6 TM helices. Allow generous margins.
LEN_MIN, LEN_MAX = 60, 170
MIN_TM = 1

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}

RATE = 8.0
_lk = threading.Lock(); _last = [0.0]


def throttle():
    with _lk:
        w = _last[0] + 1.0 / RATE - time.time()
        if w > 0:
            time.sleep(w)
        _last[0] = time.time()


def get(url, binary=False, tries=3, timeout=60):
    for i in range(tries):
        throttle()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "radicalSAM-audit/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read() if binary else r.read().decode("utf-8", "replace")
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(1.0 * (i + 1))
    return None


def tm_helices(seq, window=19, thresh=1.6):
    """Count Kyte-Doolittle hydrophobic stretches (same params as the paper)."""
    if len(seq) < window:
        return 0
    vals = []
    for i in range(len(seq) - window + 1):
        w = seq[i:i + window]
        vals.append(sum(KD.get(c, 0) for c in w) / window)
    n, inside = 0, False
    for v in vals:
        if v >= thresh and not inside:
            n += 1
            inside = True
        elif v < thresh:
            inside = False
    return n


def stage_filter():
    kept = {}
    pid, seq = None, []
    def flush():
        if pid and seq:
            s = "".join(seq)
            if LEN_MIN <= len(s) <= LEN_MAX and tm_helices(s) >= MIN_TM:
                kept[pid] = s
    for line in open(SRC):
        if line.startswith(">"):
            flush()
            pid, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    flush()
    out = os.path.join(ST, "sweep_filtered.faa")
    with open(out, "w") as fh:
        for k, v in kept.items():
            fh.write(f">{k}\n{v}\n")
    print(f"[filter] {len(kept)} of the dereplicated set match the DUF4229 profile "
          f"({LEN_MIN}-{LEN_MAX} aa, >={MIN_TM} TM helix)")
    return out, kept


def stage_map(ids):
    """Map GenBank/EMBL CDS protein IDs to UniProt accessions.

    Uses UniProt's dedicated ID-mapping service. An earlier version built
    'xref:embl-X OR xref:embl-Y ...' search URLs in batches of 150, which
    silently resolved only ~3% because long OR-queries exceed URL/result
    limits. The ID-mapping endpoint is the correct tool and is reliable at
    scale.
    """
    ids = [i.split(".")[0] for i in ids]
    out = {}
    B = 800
    for i in range(0, len(ids), B):
        chunk = ids[i:i + B]
        try:
            data = urllib.parse.urlencode({
                "from": "EMBL-GenBank-DDBJ_CDS", "to": "UniProtKB",
                "ids": ",".join(chunk)}).encode()
            throttle()
            job = json.load(urllib.request.urlopen(urllib.request.Request(
                "https://rest.uniprot.org/idmapping/run", data=data), timeout=90))["jobId"]
        except Exception as e:
            print(f"  [map] submit failed: {e}")
            continue
        # poll
        for _ in range(60):
            time.sleep(2.0)
            try:
                st = json.load(urllib.request.urlopen(
                    f"https://rest.uniprot.org/idmapping/status/{job}", timeout=60))
            except Exception:
                continue
            if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
                break
        # paginate results
        url = (f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job}"
               f"?format=json&size=500")
        while url:
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "radicalSAM-audit/1.0"})
                throttle()
                with urllib.request.urlopen(req, timeout=90) as r:
                    res = json.loads(r.read().decode())
                    link = r.headers.get("Link", "")
            except Exception:
                break
            for x in res.get("results", []):
                frm = x.get("from")
                to = x.get("to")
                acc = to.get("primaryAccession") if isinstance(to, dict) else to
                if frm and acc:
                    out[frm] = acc
            url = None
            if 'rel="next"' in link:
                url = link.split("<")[1].split(">")[0]
        print(f"  [map] {min(i+B,len(ids))}/{len(ids)} mapped={len(out)}", flush=True)
    print(f"[map] {len(out)}/{len(ids)} protein IDs resolved to UniProt "
          f"({100.0*len(out)/max(1,len(ids)):.1f}%)")
    return out


def afdb(acc, pid):
    dest = os.path.join(PDB, f"{pid}.pdb")
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return True
    meta = get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
    if not meta:
        return False
    try:
        url = json.loads(meta)[0].get("pdbUrl")
    except Exception:
        return False
    if not url:
        return False
    d = get(url, binary=True, timeout=90)
    if not d or b"ATOM" not in d[:5000]:
        return False
    open(dest, "wb").write(d)
    return True


def stage_fetch(mapping):
    os.makedirs(PDB, exist_ok=True)
    items = list(mapping.items())
    print(f"[afdb] fetching up to {len(items)} structures")
    ok = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(afdb, acc, pid): pid for pid, acc in items}
        for n, f in enumerate(as_completed(futs), 1):
            try:
                if f.result():
                    ok += 1
            except Exception:
                pass
            if n % 200 == 0 or n == len(futs):
                print(f"  [afdb] {n}/{len(futs)} got={ok}", flush=True)
    print(f"[afdb] {ok} structures downloaded")
    return ok


def stage_foldseek():
    res = os.path.join(ST, "sweep_foldseek.tsv")
    cmd = ["foldseek", "easy-search", REF, PDB, res, os.path.join(ST, "tmp_sweep"),
           "--format-output", "query,target,fident,alnlen,evalue,bits,alntmscore",
           "-e", "1000", "--max-seqs", "20000"]
    subprocess.run(cmd, check=True, capture_output=True)
    rows = [l.rstrip("\n").split("\t") for l in open(res)] if os.path.exists(res) else []
    rows = [r for r in rows if len(r) >= 7]
    rows.sort(key=lambda r: -float(r[6]))
    print(f"\n[foldseek] {len(rows)} structural hits vs DUF4229 reference")
    print(f"  {'target':<20}{'TMscore':<10}{'evalue':<12}{'bits':<8}{'fident'}")
    for r in rows[:20]:
        print(f"  {r[1][:19]:<20}{float(r[6]):<10.3f}{r[4]:<12}{r[5]:<8}{r[2]}")
    strong = [r for r in rows if float(r[6]) >= 0.5]
    print(f"\n  hits with TM-score >= 0.5 (same-fold threshold): {len(strong)}")
    if not strong:
        print("  -> NO background protein adopts the DUF4229 fold.")
        print("     The background zero survives a structural search, not just sequence.")
    else:
        for r in strong:
            print(f"     {r[1]} TM={r[6]} e={r[4]}  <-- INVESTIGATE")
    return rows


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.makedirs(ST, exist_ok=True)
    filt, kept = stage_filter()
    if stage == "filter":
        return
    mapping = stage_map(kept.keys())
    with open(os.path.join(ST, "sweep_map.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["protein_id", "uniprot"])
        for k, v in mapping.items():
            w.writerow([k, v])
    if stage == "map":
        return
    stage_fetch(mapping)
    stage_foldseek()


if __name__ == "__main__":
    main()
