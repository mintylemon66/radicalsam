#!/usr/bin/env python3
"""
Gene-cluster architecture around DUF4229 (the "check gene clusters" item).

Distance and strand tell us DUF4229 sits ~2 genes from the radical SAM gene on
the same strand. This asks the next question: WHAT ELSE is in that cluster, and
how conserved is its composition?

Three outputs:
  (a) positional profile: which Pfam families occupy each gene offset relative
      to the radical SAM gene, in DUF4229-positive vs DUF4229-negative genomes;
  (b) cluster composition: which families co-occur with DUF4229 specifically,
      i.e. enriched in DUF4229-positive target genomes over DUF4229-negative
      ones (this controls for "futalosine pathway" and isolates the DUF4229
      sub-architecture);
  (c) a consensus operon diagram for the most common architecture.

All local; no network.
"""
import csv
import os
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
CATDIR = os.path.join(AUD, "cat_domtbl")
DUF_DOM = os.path.join(AUD, "snap", "PF14012_vs_neighbors.domtbl")
RS_DOM = os.path.join(AUD, "PF04055_vs_all_neighbors.domtbl")

WIN_GENES = 10


def load_dom(path):
    s = set()
    if not os.path.exists(path):
        return s
    for line in open(path):
        if not line.startswith("#"):
            f = line.split()
            if len(f) >= 8:
                s.add(f[0])
    return s


def main():
    duf = load_dom(DUF_DOM)
    rsam = load_dom(RS_DOM)
    pf = {}
    names = {}
    if os.path.isdir(CATDIR):
        for fn in os.listdir(CATDIR):
            if not fn.endswith(".domtbl"):
                continue
            acc = fn[:-len(".domtbl")]
            pf[acc] = load_dom(os.path.join(CATDIR, fn))
            for cand in (os.path.join(ROOT, "hmm", "cat", acc + ".hmm"),
                         os.path.join(ROOT, "hmm", acc + ".hmm")):
                if os.path.exists(cand):
                    for line in open(cand, errors="replace"):
                        if line.startswith("NAME"):
                            names[acc] = line.split()[1]
                            break
                    break
    pf["PF14012"] = duf
    names["PF14012"] = "DUF4229"
    pf["PF04055"] = rsam
    names["PF04055"] = "Radical_SAM"

    rows = defaultdict(list)
    grp = {}
    with open(TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grp[r["query_acc"]] = r["group"]
            rows[r["query_acc"]].append(r)

    def gi(r):
        v = r["gene_index"]
        try:
            return int(v)
        except (TypeError, ValueError):
            return None

    # classify target genomes by DUF4229 presence within the window
    pos, neg = [], []
    for acc, rr in rows.items():
        if grp.get(acc) != "mqn":
            continue
        has = any(r["n_protein_id"] in duf and gi(r) is not None
                  and abs(gi(r)) <= WIN_GENES for r in rr)
        (pos if has else neg).append(acc)
    print(f"target genomes: DUF4229-positive={len(pos)}  DUF4229-negative={len(neg)}")

    # ---------- (a) positional profile ----------
    print("\n" + "=" * 76)
    print("(a) POSITIONAL PROFILE in DUF4229-positive genomes")
    print("    most common Pfam at each offset from the radical SAM gene")
    print("=" * 76)
    at = defaultdict(Counter)
    posset = set(pos)
    for acc in pos:
        for r in rows[acc]:
            g = gi(r)
            if g is None or abs(g) > 6:
                continue
            for acc_pf, s in pf.items():
                if r["n_protein_id"] in s:
                    at[g][acc_pf] += 1
    for off in range(-6, 7):
        c = at.get(off)
        if not c:
            continue
        top = "; ".join(f"{names.get(a,a)}={n}" for a, n in c.most_common(3))
        star = "  <== DUF4229 peak" if off in (2, -2) else ""
        print(f"  offset {off:>+3} | {top}{star}")

    # ---------- (b) what co-occurs WITH DUF4229 specifically ----------
    print("\n" + "=" * 76)
    print("(b) FAMILIES ENRICHED IN DUF4229-POSITIVE vs DUF4229-NEGATIVE TARGETS")
    print("    (both groups are futalosine genomes, so this isolates the")
    print("     DUF4229 sub-architecture rather than the pathway itself)")
    print("=" * 76)
    from scipy.stats import fisher_exact
    res = []
    for acc_pf, s in pf.items():
        if acc_pf == "PF14012":
            continue
        kp = sum(1 for a in pos if any(r["n_protein_id"] in s and gi(r) is not None
                                       and abs(gi(r)) <= WIN_GENES for r in rows[a]))
        kn = sum(1 for a in neg if any(r["n_protein_id"] in s and gi(r) is not None
                                       and abs(gi(r)) <= WIN_GENES for r in rows[a]))
        if not (pos and neg):
            continue
        _, p = fisher_exact([[kp, len(pos) - kp], [kn, len(neg) - kn]],
                            alternative="greater")
        pp, pn = 100.0 * kp / len(pos), 100.0 * kn / len(neg)
        res.append((names.get(acc_pf, acc_pf), acc_pf, kp, pp, kn, pn,
                    (pp / pn if pn else float("inf")), p))
    res.sort(key=lambda x: x[7])
    print(f"  {'family':<18}{'Pfam':<10}{'DUF+':<16}{'DUF-':<16}{'enrich':<9}{'p'}")
    for nm, a, kp, pp, kn, pn, e, p in res[:12]:
        es = "inf" if e == float("inf") else f"{e:.1f}x"
        print(f"  {nm:<18}{a:<10}{f'{kp} ({pp:.0f}%)':<16}{f'{kn} ({pn:.0f}%)':<16}{es:<9}{p:.2e}")

    # ---------- (c) consensus architecture ----------
    print("\n" + "=" * 76)
    print("(c) CONSENSUS CLUSTER ARCHITECTURE (DUF4229-positive genomes)")
    print("=" * 76)
    arch = Counter()
    for acc in pos:
        slot = {}
        for r in rows[acc]:
            g = gi(r)
            if g is None or abs(g) > 3:
                continue
            lbl = None
            if r["n_protein_id"] in duf:
                lbl = "DUF4229"
            elif r["n_protein_id"] in rsam:
                lbl = "radSAM"
            if lbl:
                slot[g] = lbl
        key = tuple(sorted(slot.items()))
        if key:
            arch[key] += 1
    for key, n in arch.most_common(6):
        line = "  " + "  ".join(f"[{o:+d}:{l}]" for o, l in key)
        print(f"{line}   x{n}")
    print("\n  (offset 0 is the query radical SAM gene itself)")


if __name__ == "__main__":
    main()
