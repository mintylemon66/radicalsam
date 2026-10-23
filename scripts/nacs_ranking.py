#!/usr/bin/env python3
"""Novel Accessory Candidate Score (Methods 2.4), computed, plus the EFI-GNT comparison.

    NACS = log2(enrichment) * novelty * log2(recurrence)

Until now NACS was described in Methods but never implemented, so the claim that
DUF4229 "ranks first" was an assertion. This computes it from raw data and ranks
the same candidate pool that EFI-GNT ranked, so the two are directly comparable.

Window: gap <= 8000 bp AND |gene offset| <= 10, anchor removed BY COORDINATES
(see ANCHOR_CONVENTION_2026-07-26.md), matching every other canonical count.

Run: python3 scripts/nacs_ranking.py
Writes audit/nacs_ranking.tsv only with --write.
"""
import collections
import csv
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WINDOW_BP, WINDOW_GENES = 8000, 10

# novelty is a binary indicator: a family whose only name is "DUF####"/"UPF####"
# has no assigned function, so it scores 1. Anything with a functional name is
# already characterized and is down-weighted to NOVELTY_KNOWN. The value is a
# stated modelling choice, not a fitted parameter; the ranking is reported for
# several values in the sensitivity block below.
NOVELTY_KNOWN = 0.25
MIN_GENOMES = 10          # do not rank a family seen in fewer target genomes
UNCHAR = re.compile(r"^(DUF|UPF)\d+", re.I)


def load_anchor_map():
    """(query, neighbour_protein) -> True if that neighbour IS the anchor gene."""
    anch, genus = {}, {}
    with open(os.path.join(ROOT, "audit", "neighbors_seq.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                is_anchor = (abs(int(r["n_start"]) - int(r["q_start"])) <= 1
                             and abs(int(r["n_end"]) - int(r["q_end"])) <= 1)
            except (ValueError, KeyError):
                is_anchor = False
            anch[(r["query_acc"], r["n_protein_id"])] = is_anchor
    return anch


def load_genus():
    g = {}
    with open(os.path.join(ROOT, "data", "metadata_org.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r.get("Sequence_ID", "")
            m = re.match(r"^(?:sp|tr)\|([^|]+)\|", sid)
            acc = m.group(1) if m else sid
            org = (r.get("Organism_Name") or "").strip()
            if org:
                g[acc] = org.split()[0]
    return g


def main():
    anch = load_anchor_map()
    genus_of = load_genus()

    fam_name = {}
    hits = collections.defaultdict(lambda: collections.defaultdict(set))   # pfam -> group -> genomes
    fam_genera = collections.defaultdict(set)
    all_q = collections.defaultdict(set)

    with open(os.path.join(ROOT, "audit", "pfam_wide.domtbl")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 23:
                continue
            acc = p[4].split(".")[0] if p[4] != "-" else p[3]
            fam_name[acc] = p[3]
            d = p[-1].split("|")
            if len(d) < 4:
                continue
            q, grp = d[0], d[1]
            try:
                gi = int(d[2].split("=")[1])
                gap = int(d[3].split("=")[1])
            except (IndexError, ValueError):
                continue
            all_q[grp].add(q)
            if gap > WINDOW_BP or abs(gi) > WINDOW_GENES:
                continue
            if anch.get((q, p[0]), False):     # the anchor enzyme itself
                continue
            hits[acc][grp].add(q)
            if grp == "mqn" and q in genus_of:
                fam_genera[acc].add(genus_of[q])

    n_t, n_b = len(all_q["mqn"]), len(all_q["bg"])
    print(f"target genomes {n_t:,} | background genomes {n_b:,}")

    rows = []
    for acc, grps in hits.items():
        k1, k2 = len(grps.get("mqn", ())), len(grps.get("bg", ()))
        if k1 < MIN_GENOMES:
            continue
        enr = (k1 / n_t) / ((k2 / n_b) if k2 else (0.5 / n_b))
        rec = len(fam_genera.get(acc, ())) or 1
        name = fam_name.get(acc, acc)
        nov = 1.0 if UNCHAR.match(name) else NOVELTY_KNOWN
        nacs = math.log2(max(enr, 1e-9)) * nov * math.log2(max(rec, 2))
        rows.append({"pfam": acc, "name": name, "k1": k1, "n1": n_t, "k2": k2, "n2": n_b,
                     "enrichment": round(enr, 3), "genera": rec,
                     "novelty": nov, "nacs": round(nacs, 3)})

    rows.sort(key=lambda r: -r["nacs"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    print(f"families ranked (>= {MIN_GENOMES} target genomes): {len(rows)}\n")
    print(f"{'rank':>4}  {'pfam':<10}{'name':<22}{'k1/n1':>12}{'enrich':>10}{'genera':>7}{'NACS':>9}")
    for r in rows[:15]:
        mark = "  <== DUF4229" if r["pfam"] == "PF14012" else ""
        print(f"{r['rank']:>4}  {r['pfam']:<10}{r['name'][:21]:<22}"
              f"{str(r['k1'])+'/'+str(r['n1']):>12}{r['enrichment']:>10.1f}"
              f"{r['genera']:>7}{r['nacs']:>9.2f}{mark}")

    duf = next((r for r in rows if r["pfam"] == "PF14012"), None)
    print(f"\nDUF4229 (PF14012): rank {duf['rank']} of {len(rows)}" if duf else "\nDUF4229 not ranked")

    # sensitivity: does the top rank depend on the novelty down-weight?
    print("\nsensitivity of DUF4229's rank to the novelty down-weight:")
    for nk in (0.1, 0.25, 0.5, 0.75, 1.0):
        tmp = []
        for r in rows:
            nov = 1.0 if UNCHAR.match(r["name"]) else nk
            tmp.append((math.log2(max(r["enrichment"], 1e-9)) * nov * math.log2(max(r["genera"], 2)), r["pfam"]))
        tmp.sort(key=lambda t: -t[0])
        rk = [i for i, (_, p) in enumerate(tmp, 1) if p == "PF14012"][0]
        print(f"  novelty(known) = {nk:<5} -> DUF4229 rank {rk}")

    if "--write" in sys.argv:
        out = os.path.join(ROOT, "audit", "nacs_ranking.tsv")
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t",
                               fieldnames=["rank", "pfam", "name", "k1", "n1", "k2", "n2",
                                           "enrichment", "genera", "novelty", "nacs"])
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"\nwrote {os.path.relpath(out, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
