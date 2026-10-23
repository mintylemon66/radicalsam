#!/usr/bin/env python3
"""
Re-derive the accessory-category enrichment table and its Benjamini-Hochberg
FDR from SEQUENCE, replacing the original product-name text groupings.

Why: the original category table (menaquinone pathway, shikimate, ferredoxin,
Fe-S assembly, ...) was built by matching free-text GenBank product names, the
same method shown to be only 56.9% sensitive for DUF4229. Its q-values are
therefore not trustworthy. Here each category is defined by an explicit set of
Pfam domains, every neighbor protein is scanned with HMMER at the Pfam
gathering threshold, and enrichment + Fisher + BH-FDR are recomputed.

Categories deliberately EXCLUDED: the original "reductase/oxidoreductase"
grouping has no clean domain definition (it spans dozens of unrelated
families), so re-deriving it would be arbitrary. It is dropped rather than
faked. Pfams PF00378 (ECH_1) and PF00501 (AMP-binding) were downloaded but are
excluded from the menaquinone category because both are large generic
metabolic families that would inflate it.
"""
import csv
import os
import subprocess
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
FAA = os.path.join(AUD, "snap", "neighbors_seq.faa")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
HMMDIR = os.path.join(ROOT, "hmm", "cat")
OUT = os.path.join(AUD, "category_fdr_sequence_based.tsv")

ORIG_BP, ORIG_GENES = 8000, 10

CATEGORIES = {
    "menaquinone / futalosine pathway (validation)": ["PF02621", "PF01048"],
    "prenyltransferase (UbiA)":                      ["PF01040"],
    "shikimate / chorismate precursor":              ["PF01202", "PF01488", "PF00275",
                                                      "PF01761", "PF01487"],
    "thioredoxin / glutaredoxin (redox)":            ["PF00085", "PF00462", "PF13098"],
    "Fe-S cluster assembly":                         ["PF01592", "PF01106", "PF01521",
                                                      "PF00266"],
    "ferredoxin (classic electron donor)":           ["PF00037", "PF00111", "PF12838",
                                                      "PF13187"],
    "transporter / membrane":                        ["PF00005", "PF00664", "PF07690"],
    "DUF4229 (discovery)":                           ["PF14012"],
}


def hmm_path(pf):
    for p in (os.path.join(HMMDIR, f"{pf}.hmm"), os.path.join(ROOT, "hmm", f"{pf}.hmm")):
        if os.path.exists(p):
            return p
    return None


def scan(pf):
    """Return set of protein IDs matching this Pfam at the gathering threshold."""
    dom = os.path.join(AUD, "cat_domtbl", f"{pf}.domtbl")
    os.makedirs(os.path.dirname(dom), exist_ok=True)
    hp = hmm_path(pf)
    if hp is None:
        print(f"  [warn] no HMM for {pf}")
        return set()
    if not (os.path.exists(dom) and os.path.getmtime(dom) > os.path.getmtime(FAA)):
        subprocess.run(["hmmsearch", "--cut_ga", "--noali", "-o", os.devnull,
                        "--domtblout", dom, hp, FAA], check=True)
    hits = set()
    for line in open(dom):
        if not line.startswith("#"):
            f = line.split()
            if len(f) >= 8:
                hits.add(f[0])
    return hits


def main():
    print("[scan] running HMMER for each category Pfam (gathering threshold)")
    pf_hits = {}
    for pfs in CATEGORIES.values():
        for pf in pfs:
            if pf not in pf_hits:
                pf_hits[pf] = scan(pf)
                print(f"    {pf}: {len(pf_hits[pf])} protein hits")

    rows = []
    grp = {}
    with open(TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grp[r["query_acc"]] = r["group"]
            rows.append(r)

    def is_anchor(r):
        """True if this neighbour row IS the anchor gene itself.

        Identify the anchor by COORDINATES, not by gene_index == 0. The index is a
        rank among coordinate-sorted CDS, and it slips whenever a neighbouring gene
        overlaps the anchor: measured across the whole table, 979 real neighbours sit
        at index 0 without being the anchor, and 981 true anchors sit at a nonzero
        index. The old index-based rule therefore deleted real neighbours and kept
        real anchors. See ANCHOR_CONVENTION_2026-07-26.md.
        """
        try:
            return (abs(int(r["n_start"]) - int(r["q_start"])) <= 1
                    and abs(int(r["n_end"]) - int(r["q_end"])) <= 1)
        except (TypeError, ValueError, KeyError):
            return False

    def in_window(r):
        try:
            gap = int(r["bp_gap"])
        except ValueError:
            return False
        gi = r["gene_index"]
        gi = int(gi) if gi not in ("", None) else 0
        return gap <= ORIG_BP and abs(gi) <= ORIG_GENES and not is_anchor(r)

    present = {c: defaultdict(bool) for c in CATEGORIES}
    for r in rows:
        if not in_window(r):
            continue
        pid = r["n_protein_id"]
        for cat, pfs in CATEGORIES.items():
            if any(pid in pf_hits[pf] for pf in pfs):
                present[cat][r["query_acc"]] = True

    n_t = sum(1 for g in grp.values() if g == "mqn")
    n_b = sum(1 for g in grp.values() if g == "bg")

    from scipy.stats import fisher_exact
    results = []
    for cat in CATEGORIES:
        kt = sum(1 for a, g in grp.items() if g == "mqn" and present[cat][a])
        kb = sum(1 for a, g in grp.items() if g == "bg" and present[cat][a])
        _, p = fisher_exact([[kt, n_t - kt], [kb, n_b - kb]], alternative="greater")
        pt, pb = 100.0 * kt / n_t, 100.0 * kb / n_b
        enr = (pt / pb) if pb > 0 else float("inf")
        results.append([cat, kt, n_t, pt, kb, n_b, pb, enr, p])

    # Benjamini-Hochberg
    order = sorted(range(len(results)), key=lambda i: results[i][8])
    m = len(results)
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        idx = m - rank + 1
        val = min(prev, results[i][8] * m / idx)
        q[i] = val
        prev = val

    print("\n" + "=" * 108)
    print("CATEGORY ENRICHMENT, SEQUENCE-BASED (Pfam HMM at gathering threshold)")
    print("=" * 108)
    print(f"{'category':<46}{'target':<16}{'background':<16}{'enrich':<10}{'p':<12}{'q (BH)':<12}")
    for i in sorted(range(m), key=lambda i: results[i][8]):
        cat, kt, nt, pt, kb, nb, pb, enr, p = results[i]
        e = "inf" if enr == float("inf") else f"{enr:.1f}x"
        sig = "" if q[i] < 0.05 else "   (NS)"
        print(f"{cat:<46}{f'{kt}/{nt}={pt:.1f}%':<16}{f'{kb}/{nb}={pb:.1f}%':<16}"
              f"{e:<10}{p:<12.2e}{q[i]:<12.2e}{sig}")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "pfams", "target_hits", "target_n", "target_pct",
                    "bg_hits", "bg_n", "bg_pct", "enrichment", "p", "q_BH"])
        for i in range(m):
            cat, kt, nt, pt, kb, nb, pb, enr, p = results[i]
            w.writerow([cat, ";".join(CATEGORIES[cat]), kt, nt, f"{pt:.2f}",
                        kb, nb, f"{pb:.2f}",
                        ("inf" if enr == float("inf") else f"{enr:.2f}"),
                        f"{p:.3e}", f"{q[i]:.3e}"])
    print(f"\n[out] {OUT}")
    print("\nNOTE: 'reductase/oxidoreductase' from the original table is DROPPED:")
    print("      it has no clean Pfam definition and re-deriving it would be arbitrary.")


if __name__ == "__main__":
    main()
