#!/usr/bin/env python3
"""
THE key missing control: is the background zero futalosine-specific, or is it
merely TAXONOMIC?

DUF4229 is nearly exclusive to Actinomycetota (499/500 sequences). If the
background set contains few actinobacteria, it would score zero for reasons
having nothing to do with the futalosine pathway, and the specificity claim
would be confounded.

This script:
  1. assigns a phylum to every target and background genome (local metadata),
  2. reports the phylum composition of both sets (the confound, quantified),
  3. recomputes DUF4229 co-occurrence restricted to ACTINOBACTERIAL genomes
     only, which is the taxonomy-matched comparison,
  4. reruns Fisher on that matched subset.

If DUF4229 is still absent from actinobacterial background genomes, the
specificity claim survives and becomes much stronger.
"""
import csv
import os
import re
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
META = os.path.join(ROOT, "data", "metadata_org.tsv")
AUD = os.path.join(ROOT, "audit")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
DOMT = os.path.join(AUD, "snap", "PF14012_vs_neighbors.domtbl")

ORIG_BP = 8000
ORIG_GENES = 10
ACTINO = re.compile(r"actinomycetota|actinobacteri", re.I)


def load_taxonomy():
    """accession -> (organism, phylum)"""
    tax = {}
    with open(META, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r.get("Sequence_ID", "")
            parts = sid.split("|")
            acc = parts[1] if len(parts) >= 2 else sid
            tax[acc] = (r.get("Organism_Name", "?"), r.get("Phylum", "?"))
    return tax


def main():
    tax = load_taxonomy()
    print(f"[tax] taxonomy loaded for {len(tax)} sequences")

    hmm = set()
    for line in open(DOMT):
        if not line.startswith("#"):
            f = line.split()
            if len(f) >= 8:
                hmm.add(f[0])

    grp = {}
    pos_orig = defaultdict(bool)
    pos_wide = defaultdict(bool)
    for r in csv.DictReader(open(TSV), delimiter="\t"):
        acc = r["query_acc"]
        grp[acc] = r["group"]
        if r["n_protein_id"] in hmm:
            pos_wide[acc] = True
            try:
                gap = int(r["bp_gap"])
            except ValueError:
                gap = 10 ** 9
            gi = r["gene_index"]
            gi = int(gi) if gi not in ("", None) else 0
            if gap <= ORIG_BP and abs(gi) <= ORIG_GENES:
                pos_orig[acc] = True

    # ---- phylum composition ----
    comp = {"mqn": Counter(), "bg": Counter()}
    unknown = Counter()
    for acc, g in grp.items():
        if g not in comp:
            continue
        t = tax.get(acc)
        if not t:
            unknown[g] += 1
            comp[g]["<no taxonomy>"] += 1
        else:
            comp[g][t[1] or "?"] += 1

    print("\n" + "=" * 72)
    print("PHYLUM COMPOSITION of the two sets (this is the confound)")
    print("=" * 72)
    for g in ("mqn", "bg"):
        tot = sum(comp[g].values())
        print(f"\n  {g}  (n={tot})")
        for ph, n in comp[g].most_common(10):
            flag = "  <-- ACTINO" if ACTINO.search(ph or "") else ""
            print(f"    {ph[:44]:<46}{n:>6}  ({100.0*n/tot:5.1f}%){flag}")

    def is_actino(acc):
        t = tax.get(acc)
        return bool(t and ACTINO.search(t[1] or ""))

    # ---- taxonomy-matched contingency ----
    print("\n" + "=" * 72)
    print("TAXONOMY-MATCHED COMPARISON (actinobacterial genomes only)")
    print("=" * 72)

    for label, table in (("original window (<=8kb)", pos_orig),
                         ("wide window (<=25kb)", pos_wide)):
        rows = {}
        for setname, want_actino in (("ALL genomes", None), ("ACTINO only", True),
                                     ("NON-actino only", False)):
            t = {"mqn": [0, 0], "bg": [0, 0]}
            for acc, g in grp.items():
                if g not in t:
                    continue
                if want_actino is not None and is_actino(acc) != want_actino:
                    continue
                t[g][1] += 1
                if table.get(acc):
                    t[g][0] += 1
            rows[setname] = t

        print(f"\n--- {label} ---")
        print(f"  {'subset':<18}{'target':<22}{'background':<22}{'Fisher p':<12}{'OR':<10}")
        for setname, t in rows.items():
            m, b = t["mqn"], t["bg"]
            ms = f"{m[0]}/{m[1]}" + (f" = {100.0*m[0]/m[1]:.2f}%" if m[1] else "")
            bs = f"{b[0]}/{b[1]}" + (f" = {100.0*b[0]/b[1]:.2f}%" if b[1] else "")
            p = orr = "n/a"
            if m[1] and b[1]:
                try:
                    from scipy.stats import fisher_exact
                    _, pv = fisher_exact([[m[0], m[1] - m[0]],
                                          [b[0], b[1] - b[0]]], alternative="greater")
                    p = f"{pv:.2e}"
                    a, bb, c, d = m[0]+.5, m[1]-m[0]+.5, b[0]+.5, b[1]-b[0]+.5
                    orr = f"{(a/bb)/(c/d):.1f}x"
                except ImportError:
                    pass
            print(f"  {setname:<18}{ms:<22}{bs:<22}{p:<12}{orr:<10}")

    # ---- the decisive number ----
    actino_bg = [a for a, g in grp.items() if g == "bg" and is_actino(a)]
    actino_bg_pos = [a for a in actino_bg if pos_wide.get(a)]
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    print(f"  actinobacterial BACKGROUND genomes available : {len(actino_bg)}")
    print(f"  of those, carrying DUF4229 (wide window)     : {len(actino_bg_pos)}")
    if len(actino_bg) < 30:
        print("  -> TOO FEW actinobacterial background genomes for a powered test.")
        print("     The background is taxonomically mismatched, and a purpose-built")
        print("     actinobacterial background must be constructed. THIS IS A REAL")
        print("     LIMITATION and must be stated in the paper.")
    elif not actino_bg_pos:
        print("  -> DUF4229 is absent even from taxonomy-matched background genomes.")
        print("     The specificity claim is NOT a taxonomic artifact. Claim survives.")
    else:
        print("  -> DUF4229 DOES occur in actinobacterial background genomes.")
        print("     The background zero was partly taxonomic. Effect size must be")
        print("     recomputed on the matched subset and the claim softened.")


if __name__ == "__main__":
    main()
