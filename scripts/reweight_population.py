#!/usr/bin/env python3
"""
Correct the stratified-sample rates back to a POPULATION estimate.

The target-side fetch deliberately over-sampled label-positive genomes (all
200 of them) plus a random draw of label-negatives. Raw rates from that
sample are therefore inflated and must NOT be reported directly.

Population estimate:
    rate = [ N_pos * P(HMM+ | label+) + N_neg * P(HMM+ | label-) ] / N_total

where N_pos / N_neg come from the FULL original dataset, and the two
conditional probabilities are measured on whatever subset we actually
fetched. This is standard post-stratification.
"""
import csv
import os
import re
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FULL = os.path.join(ROOT, "data", "mqn_cooccur_FULL.tsv")
AUD = os.path.join(ROOT, "audit")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
DOMT = os.path.join(AUD, "snap", "PF14012_vs_neighbors.domtbl")

ORIG_BP = 8000
ORIG_GENES = 10


def main():
    # --- population strata from the ORIGINAL full table ---
    label_pos_all, label_neg_all = set(), set()
    for r in csv.DictReader(open(FULL), delimiter="\t"):
        if r["group"].strip() != "mqn":
            continue
        a = r["accession"].strip()
        if re.search(r"DUF4229", r["products"] or "", re.I):
            label_pos_all.add(a)
        else:
            label_neg_all.add(a)
    N_pos, N_neg = len(label_pos_all), len(label_neg_all)
    N_tot = N_pos + N_neg
    print(f"POPULATION (original full table): mqn total={N_tot}  "
          f"label+={N_pos}  label-={N_neg}")

    # --- HMM hits ---
    hmm = set()
    for line in open(DOMT):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) >= 8:
            hmm.add(f[0])

    # --- per-genome status in the fetched sample ---
    wide = defaultdict(bool)
    orig = defaultdict(bool)
    grp = {}
    for r in csv.DictReader(open(TSV), delimiter="\t"):
        if r["group"] != "mqn":
            continue
        acc = r["query_acc"]
        grp[acc] = "mqn"
        if r["n_protein_id"] in hmm:
            wide[acc] = True
            try:
                gap = int(r["bp_gap"])
            except ValueError:
                gap = 10 ** 9
            gi = r["gene_index"]
            gi = int(gi) if gi not in ("", None) else 0
            if gap <= ORIG_BP and abs(gi) <= ORIG_GENES:
                orig[acc] = True

    fetched = set(grp)
    f_pos = fetched & label_pos_all
    f_neg = fetched & label_neg_all
    print(f"FETCHED: {len(fetched)} mqn genomes  (label+ {len(f_pos)}, label- {len(f_neg)})")

    def cond(subset, table):
        if not subset:
            return None
        return sum(1 for a in subset if table.get(a)) / len(subset)

    for name, table in (("original window (<=8kb)", orig), ("wide window (<=25kb)", wide)):
        p_pos = cond(f_pos, table)
        p_neg = cond(f_neg, table)
        print(f"\n--- {name} ---")
        print(f"  P(HMM+ | label+) = {p_pos:.4f}" if p_pos is not None else "  no label+ fetched")
        print(f"  P(HMM+ | label-) = {p_neg:.4f}" if p_neg is not None else "  no label- fetched")
        if p_pos is None or p_neg is None:
            continue
        est = (N_pos * p_pos + N_neg * p_neg) / N_tot
        # binomial SE on the label- stratum dominates the uncertainty
        import math
        se_neg = math.sqrt(max(p_neg * (1 - p_neg), 1e-12) / len(f_neg))
        se = (N_neg / N_tot) * se_neg
        lo, hi = max(0.0, est - 1.96 * se), est + 1.96 * se
        print(f"  POPULATION ESTIMATE = {100*est:.2f}%  (95% CI {100*lo:.2f}% to {100*hi:.2f}%)")
        print(f"  implied genome count = {est*N_tot:.0f} / {N_tot}")

    print("\nORIGINAL PUBLISHED FIGURE: 200/3647 = 5.48% (text label, <=8kb window)")


if __name__ == "__main__":
    main()
