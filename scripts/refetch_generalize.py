#!/usr/bin/env python3
"""Rebuild BioB / LipA neighbourhoods from sequence, mirroring the MqnC/MqnE refetch.

The generalization results in Results.md 3.5 were produced by the ORIGINAL
product-name text pipeline and were never re-derived from sequence, unlike every
other number in the paper. This closes that inconsistency.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import refetch_neighborhoods_seq as R

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R.SRC     = os.path.join(ROOT, "audit", "generalize_targets.tsv")
R.OUT_TSV = os.path.join(ROOT, "audit", "generalize_neighbors_seq.tsv")
R.OUT_FAA = os.path.join(ROOT, "audit", "generalize_neighbors_seq.faa")
R.LOG     = os.path.join(ROOT, "audit", "generalize_refetch.log")

import csv
_orig = R.main
def main():
    targets = []
    with open(R.SRC, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a = (row.get("accession") or "").strip()
            g = (row.get("group") or "").strip()
            if a and g in ("BioB", "LipA"):
                targets.append((a, g))
    seen = R.done_set()
    targets = [t for t in targets if t[0] not in seen]
    if not os.path.exists(R.OUT_TSV):
        with open(R.OUT_TSV, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerow(R.COLS)
    print(f"[gen] todo={len(targets)} window=+/-{R.WINDOW}bp", flush=True)
    mapping = {}
    B = 60
    for i in range(0, len(targets), B):
        mapping.update(R.uniprot_batch([a for a, _ in targets[i:i+B]]))
        if (i // B) % 10 == 0:
            print(f"  [uniprot] {min(i+B,len(targets))}/{len(targets)} mapped={len(mapping)}", flush=True)
    print(f"[gen] mapped {len(mapping)}/{len(targets)}", flush=True)
    import threading
    from concurrent.futures import ThreadPoolExecutor, as_completed
    wlock = threading.Lock()
    ok = rows = 0
    with ThreadPoolExecutor(max_workers=R.WORKERS) as ex:
        futs = {}
        for acc, grp in targets:
            if acc in mapping:
                contig, pid = mapping[acc]
                futs[ex.submit(R.process, acc, grp, contig, pid, wlock)] = acc
        for n, f in enumerate(as_completed(futs), 1):
            try:
                r = f.result()
            except Exception:
                r = 0
            if r: ok += 1; rows += r
            if n % 200 == 0 or n == len(futs):
                print(f"  [region] {n}/{len(futs)} ok={ok} rows={rows}", flush=True)
    print(f"[gen] DONE genomes_ok={ok} neighbor_rows={rows}", flush=True)

if __name__ == "__main__":
    main()
