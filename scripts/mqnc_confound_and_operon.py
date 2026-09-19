#!/usr/bin/env python3
"""
Check whether the DUF4229 signal near MqnC comes from an MqnE in the same operon, then
measure how often DUF4229 and the futalosine gene sit on the same strand across every pair.
"""
import csv
import os
import re
from collections import defaultdict, Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
DOMT = os.path.join(AUD, "snap", "PF14012_vs_neighbors.domtbl")
FASTA = os.path.join(ROOT, "data", "mqn_targets.fasta")
META = os.path.join(ROOT, "data", "metadata_org.tsv")

ORIG_BP, ORIG_GENES = 8000, 10
ACTINO = re.compile(r"actinomycetota|actinobacteri", re.I)
# product-name signatures of futalosine radical SAM enzymes
FUT = re.compile(r"aminofutalosine|dehypoxanthin|futalosine|mqnE|mqnC|cofG|cofH", re.I)


def classify_targets():
    cls = {}
    for line in open(FASTA):
        if not line.startswith(">"):
            continue
        p = line.split("|")
        if len(p) < 2:
            continue
        acc, d = p[1], line.lower()
        gn = re.search(r"gn=(\w+)", d)
        gn = gn.group(1) if gn else ""
        if "aminofutalosine" in d or gn in ("mqne", "cofg"):
            cls[acc] = "MqnE"
        elif "dehypoxanthin" in d or gn in ("mqnc", "cofh"):
            cls[acc] = "MqnC"
        else:
            cls[acc] = "other"
    return cls


def main():
    cls = classify_targets()
    tax = {}
    for r in csv.DictReader(open(META), delimiter="\t"):
        q = r["Sequence_ID"].split("|")
        tax[q[1] if len(q) >= 2 else r["Sequence_ID"]] = r.get("Phylum", "?")

    hmm = set()
    for line in open(DOMT):
        if not line.startswith("#"):
            f = line.split()
            if len(f) >= 8:
                hmm.add(f[0])

    rows_by_genome = defaultdict(list)
    grp = {}
    for r in csv.DictReader(open(TSV), delimiter="\t"):
        grp[r["query_acc"]] = r["group"]
        rows_by_genome[r["query_acc"]].append(r)

    def in_window(r):
        try:
            gap = int(r["bp_gap"])
        except ValueError:
            return False
        gi = r["gene_index"]
        gi = int(gi) if gi not in ("", None) else 0
        return gap <= ORIG_BP and abs(gi) <= ORIG_GENES

    # MqnC confound
    print("=" * 72)
    print("MqnC confound: is the MqnC signal just a nearby MqnE?")
    print("=" * 72)
    mqnc_pos, mqnc_pos_with_fut = [], []
    for acc, rows in rows_by_genome.items():
        if grp.get(acc) != "mqn" or cls.get(acc) != "MqnC":
            continue
        duf = [r for r in rows if r["n_protein_id"] in hmm and in_window(r)]
        if not duf:
            continue
        mqnc_pos.append(acc)
        # does this MqnC neighborhood also contain a futalosine enzyme (i.e. MqnE)?
        has_fut = any(FUT.search(r["n_product"] or "") and in_window(r)
                      and int(r["gene_index"] or 0) != 0 for r in rows)
        if has_fut:
            mqnc_pos_with_fut.append(acc)
    n, k = len(mqnc_pos), len(mqnc_pos_with_fut)
    print(f"  DUF4229-positive MqnC genomes                : {n}")
    print(f"  ...that also contain a futalosine gene nearby: {k}"
          + (f"  ({100.0*k/n:.1f}%)" if n else ""))
    if n and k / n > 0.8:
        print("  -> The MqnC signal is mostly explained by a nearby MqnE.")
        print("     DUF4229 goes with MqnE, and MqnC picks it up only by proximity.")
    elif n:
        print("  -> A substantial share of MqnC-positive genomes lack a nearby")
        print("     futalosine gene, so the MqnC signal is not purely proximity.")

    # same-strand rate across all DUF4229 pairs
    print("\n" + "=" * 72)
    print("Operon strand across all DUF4229 pairs")
    print("=" * 72)
    same = diff = 0
    gaps = []
    same_act = diff_act = 0
    for acc, rows in rows_by_genome.items():
        if grp.get(acc) != "mqn":
            continue
        is_act = bool(ACTINO.search(tax.get(acc, "") or ""))
        for r in rows:
            if r["n_protein_id"] not in hmm or not in_window(r):
                continue
            q, nstr = r["q_strand"], r["n_strand"]
            if not q or not nstr:
                continue
            if q == nstr:
                same += 1
                if is_act:
                    same_act += 1
            else:
                diff += 1
                if is_act:
                    diff_act += 1
            try:
                gaps.append(int(r["bp_gap"]))
            except ValueError:
                pass
    tot = same + diff
    print(f"  DUF4229-radicalSAM pairs measured : {tot}")
    if tot:
        print(f"  same strand                       : {same} ({100.0*same/tot:.1f}%)")
        print(f"  opposite strand                   : {diff} ({100.0*diff/tot:.1f}%)")
    ta = same_act + diff_act
    if ta:
        print(f"  (actinobacteria only: {same_act}/{ta} same-strand = {100.0*same_act/ta:.1f}%)")
    if gaps:
        gaps.sort()
        print(f"  intergenic gap: median={gaps[len(gaps)//2]} bp, "
              f"p25={gaps[len(gaps)//4]}, p75={gaps[3*len(gaps)//4]}")
        close = sum(1 for g in gaps if g <= 200)
        print(f"  gap <= 200 bp (typical operon spacing): {close}/{len(gaps)} "
              f"({100.0*close/len(gaps):.1f}%)")
    # binomial test against a 50% null
    if tot:
        try:
            from scipy.stats import binomtest
            p = binomtest(same, tot, 0.5, alternative="greater").pvalue
            print(f"  binomial test vs 50% null: p = {p:.3e}")
        except Exception:
            pass


if __name__ == "__main__":
    main()
