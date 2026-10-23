#!/usr/bin/env python3
"""
Audit 1: how much of the DUF4229 co-occurrence result depends on GenBank
product-name TEXT rather than on sequence evidence?

Reads data/mqn_cooccur_FULL.tsv (accession, group, products) and reports:
  (a) group sizes and the label-based DUF4229 counts (reproducing the headline)
  (b) the opaque-annotation pool ("hypothetical protein" etc.), which is the
      false-negative risk: a real DUF4229 protein annotated as "hypothetical"
      is invisible to a text search
  (c) how the headline numbers would move under plausible miscall rates

No network. Read-only.
"""
import csv
import re
import sys
from collections import Counter

SRC = "data/mqn_cooccur_FULL.tsv"

OPAQUE = re.compile(
    r"hypothetical|uncharacteri[sz]ed|unknown function|^duf\d+|"
    r"predicted protein|conserved protein|putative protein",
    re.I,
)


def main():
    groups = Counter()
    duf4229_hits = Counter()
    neighbor_total = 0
    opaque_total = 0
    opaque_by_group = Counter()
    neigh_by_group = Counter()
    duf_families = Counter()
    per_row_counts = []

    with open(SRC, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            grp = (row.get("group") or "").strip()
            prods = (row.get("products") or "")
            groups[grp] += 1

            items = [p.strip() for p in prods.split(";") if p.strip()]
            neighbor_total += len(items)
            neigh_by_group[grp] += len(items)

            n_duf = 0
            n_opaque = 0
            for p in items:
                if OPAQUE.search(p):
                    opaque_total += 1
                    opaque_by_group[grp] += 1
                for m in re.findall(r"DUF\d+", p, re.I):
                    duf_families[m.upper()] += 1
                if re.search(r"DUF4229", p, re.I):
                    n_duf += 1
            if n_duf:
                duf4229_hits[grp] += 1
            per_row_counts.append((grp, len(items), n_duf))

    print("=" * 68)
    print("AUDIT 1: label-based detection and its false-negative exposure")
    print("=" * 68)

    print("\n[a] Neighborhood sets")
    for g, n in groups.most_common():
        print(f"    group={g:<12} neighborhoods={n:>6}")
    print(f"    TOTAL neighborhoods = {sum(groups.values())}")

    print("\n[b] DUF4229 detected by PRODUCT-NAME TEXT MATCH")
    for g in groups:
        hit = duf4229_hits.get(g, 0)
        tot = groups[g]
        pct = 100.0 * hit / tot if tot else 0.0
        print(f"    group={g:<12} {hit:>5} / {tot:<6} = {pct:5.2f}%")

    print("\n[c] Neighbor gene entries and OPAQUE annotations (risk pool)")
    for g in groups:
        nb = neigh_by_group[g]
        op = opaque_by_group[g]
        pct = 100.0 * op / nb if nb else 0.0
        print(f"    group={g:<12} neighbors={nb:>8}  opaque={op:>8} ({pct:5.1f}%)")
    pct_all = 100.0 * opaque_total / neighbor_total if neighbor_total else 0
    print(f"    TOTAL           neighbors={neighbor_total:>8}  opaque={opaque_total:>8} ({pct_all:5.1f}%)")
    print("\n    -> Any real DUF4229 protein annotated as 'hypothetical protein'")
    print("       is INVISIBLE to the text search. This pool bounds the")
    print("       false-negative rate in BOTH the target and background sets.")

    print("\n[d] All DUF families seen by text match (top 20)")
    for fam, n in duf_families.most_common(20):
        print(f"    {fam:<12} {n}")

    print("\n[e] Mean neighborhood size by group (sanity: comparable windows?)")
    for g in groups:
        nb = neigh_by_group[g]
        tot = groups[g]
        print(f"    group={g:<12} mean genes/neighborhood = {nb/tot:5.1f}")

    print("\n[f] WHY THIS MATTERS FOR THE HEADLINE")
    mqn_key = "mqn" if "mqn" in groups else None
    bg_keys = [g for g in groups if g != mqn_key]
    if mqn_key and bg_keys:
        bg = bg_keys[0]
        t_hit, t_tot = duf4229_hits.get(mqn_key, 0), groups[mqn_key]
        b_hit, b_tot = duf4229_hits.get(bg, 0), groups[bg]
        print(f"    label-based target : {t_hit}/{t_tot}")
        print(f"    label-based bkgrnd : {b_hit}/{b_tot}")
        print("    The background ZERO is the load-bearing number (it drives")
        print("    the odds ratio and the p-value). If even a handful of the")
        print("    background 'hypothetical' neighbors are truly DUF4229, the")
        print("    zero becomes nonzero and the effect size drops sharply.")
        print("    -> This is why the sequence-based (HMM) recount is required.")


if __name__ == "__main__":
    sys.exit(main())
