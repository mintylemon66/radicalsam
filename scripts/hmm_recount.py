#!/usr/bin/env python3
"""
Count DUF4229 in the neighborhoods from sequence, using HMMER against PF14012 at the
Pfam gathering threshold, and compare those calls with the GenBank product names. The
refetched neighborhoods span 25 kb on each side, so this also measures how far DUF4229
sits from the query gene and whether the wider window changes any genome's call.
"""
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
FAA = os.path.join(AUD, os.environ.get("SNAP",""), "neighbors_seq.faa")
TSV = os.path.join(AUD, os.environ.get("SNAP",""), "neighbors_seq.tsv")
HMM = os.path.join(ROOT, "hmm", "PF14012.hmm")
DOMT = os.path.join(AUD, os.environ.get("SNAP",""), "PF14012_vs_neighbors.domtbl")

ORIG_BP = 8000      # original pipeline window
ORIG_GENES = 10     # original pipeline gene window


def run_hmmsearch():
    if os.path.exists(DOMT) and os.path.getmtime(DOMT) > os.path.getmtime(FAA):
        print(f"[hmm] reusing {os.path.basename(DOMT)}")
        return
    cmd = ["hmmsearch", "--cut_ga", "--domtblout", DOMT, "--noali", "-o", os.devnull,
           HMM, FAA]
    print("[hmm] running:", " ".join(os.path.basename(c) for c in cmd))
    subprocess.run(cmd, check=True)


def parse_domtbl():
    hits = {}
    if not os.path.exists(DOMT):
        return hits
    with open(DOMT) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 14:
                continue
            pid, full_ev, full_sc = f[0], float(f[6]), float(f[7])
            prev = hits.get(pid)
            if prev is None or full_sc > prev[1]:
                hits[pid] = (full_ev, full_sc)
    return hits


def main():
    if not os.path.exists(FAA):
        sys.exit("no neighbor FASTA yet; run refetch_neighborhoods_seq.py first")

    run_hmmsearch()
    hmm_hits = parse_domtbl()
    print(f"[hmm] PF14012 sequence hits (GA threshold): {len(hmm_hits)} proteins\n")

    # load neighbor table
    rows = []
    with open(TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)

    genomes = {}
    for r in rows:
        genomes.setdefault(r["query_acc"], r["group"])

    label_pos = defaultdict(bool)  # DUF4229 by product name
    hmm_pos_wide = defaultdict(bool)  # DUF4229 by HMM within 25 kb
    hmm_pos_orig = defaultdict(bool)  # DUF4229 by HMM within the original window
    prot_label = set()
    prot_hmm = set()
    dists = []                         # (group, bp_gap, gene_index, by_label)

    for r in rows:
        g, acc = r["group"], r["query_acc"]
        pid = r["n_protein_id"]
        prod = r["n_product"]
        try:
            gap = int(r["bp_gap"])
        except ValueError:
            gap = None
        gi = r["gene_index"]
        gi = int(gi) if gi not in ("", None) else None

        is_label = bool(re.search(r"DUF4229", prod, re.I))
        is_hmm = pid in hmm_hits

        if is_label:
            label_pos[acc] = True
            prot_label.add(pid)
        if is_hmm:
            prot_hmm.add(pid)
            hmm_pos_wide[acc] = True
            in_orig = (gap is not None and gap <= ORIG_BP) and \
                      (gi is None or abs(gi) <= ORIG_GENES)
            if in_orig:
                hmm_pos_orig[acc] = True
            dists.append((g, gap, gi, is_label))

    # protein-level agreement between product name and HMM
    both = prot_label & prot_hmm
    only_label = prot_label - prot_hmm
    only_hmm = prot_hmm - prot_label
    print("=" * 70)
    print("Protein level: product name vs HMM calls")
    print("=" * 70)
    print(f"  DUF4229 by both name and HMM        : {len(both)}")
    print(f"  DUF4229 by name only                : {len(only_label)}   possible false positives")
    print(f"  DUF4229 by HMM only                 : {len(only_hmm)}   missed by the name search")
    if prot_hmm:
        rec = 100.0 * len(both) / len(prot_hmm)
        print(f"  name-search recall against HMM       : {rec:.1f}%")

    # genome-level counts
    def tally(fn):
        t = {"mqn": [0, 0], "bg": [0, 0]}   # [positive, total]
        for acc, grp in genomes.items():
            if grp not in t:
                continue
            t[grp][1] += 1
            if fn.get(acc):
                t[grp][0] += 1
        return t

    print("\n" + "=" * 70)
    print("Genome level: DUF4229 presence by method")
    print("=" * 70)
    print(f"{'method':<34}{'target (mqn)':<20}{'background (bg)':<20}")
    for name, fn in (("text label (any window)", label_pos),
                     (f"HMM, original window (<={ORIG_BP}bp)", hmm_pos_orig),
                     ("HMM, wide window (<=25kb)", hmm_pos_wide)):
        t = tally(fn)
        m, b = t["mqn"], t["bg"]
        ms = f"{m[0]}/{m[1]} = {100.0*m[0]/m[1]:.2f}%" if m[1] else "n/a"
        bs = f"{b[0]}/{b[1]} = {100.0*b[0]/b[1]:.2f}%" if b[1] else "n/a"
        print(f"{name:<34}{ms:<20}{bs:<20}")

    # Fisher test on the HMM calls
    try:
        from scipy.stats import fisher_exact
        print("\n" + "=" * 70)
        print("Fisher's exact test (one-sided) on the HMM calls")
        print("=" * 70)
        for name, fn in ((f"HMM original window", hmm_pos_orig),
                         ("HMM wide window", hmm_pos_wide)):
            t = tally(fn)
            m, b = t["mqn"], t["bg"]
            if not (m[1] and b[1]):
                continue
            table = [[m[0], m[1] - m[0]], [b[0], b[1] - b[0]]]
            odds, p = fisher_exact(table, alternative="greater")
            # Haldane-Anscombe corrected OR
            a, bb, c, d = m[0] + .5, m[1] - m[0] + .5, b[0] + .5, b[1] - b[0] + .5
            hor = (a / bb) / (c / d)
            print(f"  {name:<24} p = {p:.3g}   OR(HA-corrected) = {hor:.1f}x")
    except ImportError:
        print("(scipy unavailable; skipping Fisher)")

    # how far DUF4229 sits from the query gene
    print("\n" + "=" * 70)
    print("Distance of DUF4229 from the query radical SAM gene (HMM hits)")
    print("=" * 70)
    if dists:
        bybp = sorted(d[1] for d in dists if d[1] is not None)
        print(f"  n DUF4229 neighbor genes = {len(dists)}")
        if bybp:
            def pct(p):
                return bybp[min(len(bybp) - 1, int(p * len(bybp)))]
            print(f"  bp gap  min={bybp[0]}  p25={pct(.25)}  median={pct(.5)}  "
                  f"p75={pct(.75)}  max={bybp[-1]}")
            beyond = sum(1 for x in bybp if x > ORIG_BP)
            print(f"  beyond the original {ORIG_BP} bp window: {beyond}/{len(bybp)} "
                  f"({100.0*beyond/len(bybp):.1f}%)  <- invisible to the original run")
        gis = sorted(d[2] for d in dists if d[2] is not None)
        if gis:
            from collections import Counter
            print("\n  gene-offset histogram (negative = upstream):")
            c = Counter(gis)
            for k in sorted(c):
                if abs(k) <= 25:
                    print(f"    {k:>+4} | {'#' * min(60, c[k])} {c[k]}")

    # write the summary table
    out = os.path.join(AUD, "audit1_summary.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["method", "group", "positive", "total"])
        for name, fn in (("text_label", label_pos),
                         ("hmm_orig_window", hmm_pos_orig),
                         ("hmm_wide_window", hmm_pos_wide)):
            t = tally(fn)
            for grp in ("mqn", "bg"):
                w.writerow([name, grp, t[grp][0], t[grp][1]])
    print(f"\n[out] {out}")


if __name__ == "__main__":
    main()
