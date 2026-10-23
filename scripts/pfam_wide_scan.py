#!/usr/bin/env python3
"""
Genome-wide Pfam enrichment scan: is MqnP/UbiA genuinely DUF4229's top
partner, or merely the top of the 24-domain panel tested so far?

Everything reported about the DUF4229 gene cluster to date rests on a hand-
picked set of 24 Pfam domains. That is a real limitation: a family never
scanned for cannot be found. This scans the neighbour proteins against the
FULL Pfam-A library and ranks every family by enrichment in DUF4229-positive
versus DUF4229-negative futalosine genomes.

Two comparisons are made, and they answer different questions:
  (A) target vs background            -> what is enriched near futalosine enzymes
  (B) DUF4229-positive vs DUF4229-negative targets, both futalosine
                                      -> what is specific to the DUF4229 cluster

(B) is the one that matters for the mechanism, because it controls for the
pathway itself.

Benjamini-Hochberg FDR is applied across all tested families.
"""
import csv
import os
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
FAA = os.path.join(AUD, "snap", "neighbors_seq.faa")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
PFAM = os.path.join(ROOT, "hmm", "pfamA", "Pfam-A.hmm")
DOM = os.path.join(AUD, "pfam_wide.domtbl")
ORIG_BP, ORIG_GENES = 8000, 10
MIN_GENOMES = 15          # ignore families too rare to test meaningfully


def ensure_pfam():
    gz = PFAM + ".gz"
    if not os.path.exists(PFAM):
        if not os.path.exists(gz):
            sys.exit("Pfam-A.hmm.gz not downloaded yet")
        print("[pfam] decompressing ...", flush=True)
        subprocess.run(["gunzip", "-k", gz], check=True)
    if not os.path.exists(PFAM + ".h3i"):
        print("[pfam] hmmpress ...", flush=True)
        subprocess.run(["hmmpress", PFAM], check=True, capture_output=True)


def run_scan(threads=6):
    if os.path.exists(DOM) and os.path.getmtime(DOM) > os.path.getmtime(FAA):
        print(f"[scan] reusing {DOM}")
        return
    print("[scan] hmmsearch Pfam-A vs neighbour proteins (this is the long step)",
          flush=True)
    subprocess.run(["hmmsearch", "--cut_ga", "--noali", "-o", os.devnull,
                    "--domtblout", DOM, "--cpu", str(threads), PFAM, FAA],
                   check=True)


def main():
    ensure_pfam()
    run_scan()

    fam_of = defaultdict(set)      # pfam accession -> set of protein ids
    fam_name = {}
    with open(DOM) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            prot, fname, facc = f[0], f[3], f[4]
            acc = facc.split(".")[0] if facc != "-" else fname
            fam_of[acc].add(prot)
            fam_name[acc] = fname
    print(f"[scan] families with >=1 hit: {len(fam_of)}")

    rows = defaultdict(list)
    grp = {}
    with open(TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grp[r["query_acc"]] = r["group"]
            rows[r["query_acc"]].append(r)

    def gi(r):
        try:
            return int(r["gene_index"])
        except (TypeError, ValueError):
            return None

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

    def inwin(r):
        try:
            gap = int(r["bp_gap"])
        except ValueError:
            return False
        g = gi(r)
        return (gap <= ORIG_BP and g is not None and abs(g) <= ORIG_GENES
                and not is_anchor(r))

    duf = fam_of.get("PF14012", set())
    tgt_pos, tgt_neg, bg = [], [], []
    for acc, rr in rows.items():
        g = grp.get(acc)
        if g == "bg":
            bg.append(acc)
        elif g == "mqn":
            (tgt_pos if any(r["n_protein_id"] in duf and inwin(r) for r in rr)
             else tgt_neg).append(acc)
    print(f"[sets] target DUF4229+ {len(tgt_pos)}, target DUF4229- {len(tgt_neg)}, "
          f"background {len(bg)}")

    def present(acc_list, members):
        return sum(1 for a in acc_list
                   if any(r["n_protein_id"] in members and inwin(r) for r in rows[a]))

    from scipy.stats import fisher_exact
    resA, resB = [], []
    for acc, members in fam_of.items():
        kt = present(tgt_pos, members) + present(tgt_neg, members)
        nt = len(tgt_pos) + len(tgt_neg)
        kb = present(bg, members)
        nb = len(bg)
        if kt + kb >= MIN_GENOMES:
            _, p = fisher_exact([[kt, nt - kt], [kb, nb - kb]], alternative="greater")
            enr = ((kt / nt) / (kb / nb)) if kb else float("inf")
            resA.append((acc, fam_name[acc], kt, nt, kb, nb, enr, p))
        kp = present(tgt_pos, members)
        kn = present(tgt_neg, members)
        if kp + kn >= MIN_GENOMES:
            _, p = fisher_exact([[kp, len(tgt_pos) - kp],
                                 [kn, len(tgt_neg) - kn]], alternative="greater")
            enr = ((kp / len(tgt_pos)) / (kn / len(tgt_neg))) if kn else float("inf")
            resB.append((acc, fam_name[acc], kp, len(tgt_pos), kn, len(tgt_neg), enr, p))

    def bh(res):
        order = sorted(range(len(res)), key=lambda i: res[i][7])
        m = len(res)
        q = [0.0] * m
        prev = 1.0
        for rank, i in enumerate(reversed(order), 1):
            idx = m - rank + 1
            prev = min(prev, res[i][7] * m / idx)
            q[i] = prev
        return q

    for title, res, lab in (
        ("(A) TARGET vs BACKGROUND: enriched near futalosine enzymes", resA, "bg"),
        ("(B) DUF4229-POSITIVE vs DUF4229-NEGATIVE targets: the DUF4229 cluster", resB, "DUF-"),
    ):
        q = bh(res)
        idx = sorted(range(len(res)), key=lambda i: res[i][7])
        print("\n" + "=" * 100)
        print(title)
        print(f"families tested: {len(res)}")
        print("=" * 100)
        print(f"  {'Pfam':<10}{'name':<24}{'target':<16}{lab:<16}{'enrich':<10}{'p':<12}{'q(BH)'}")
        for i in idx[:25]:
            acc, nm, k1, n1, k2, n2, e, p = res[i]
            es = "inf" if e == float("inf") else f"{e:.1f}x"
            print(f"  {acc:<10}{nm[:23]:<24}{f'{k1}/{n1}':<16}{f'{k2}/{n2}':<16}"
                  f"{es:<10}{p:<12.1e}{q[i]:.1e}")
        out = os.path.join(AUD, f"pfam_wide_{'A' if lab=='bg' else 'B'}.tsv")
        with open(out, "w") as fh:
            fh.write("pfam\tname\tk1\tn1\tk2\tn2\tenrichment\tp\tq_BH\n")
            for i in idx:
                acc, nm, k1, n1, k2, n2, e, p = res[i]
                fh.write(f"{acc}\t{nm}\t{k1}\t{n1}\t{k2}\t{n2}\t"
                         f"{'inf' if e==float('inf') else f'{e:.3f}'}\t{p:.3e}\t{q[i]:.3e}\n")
        print(f"  [out] {out}")


if __name__ == "__main__":
    main()
