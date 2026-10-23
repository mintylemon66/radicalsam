#!/usr/bin/env python3
"""
AF3 AUDIT: verify every transcribed ipTM value, compute LIS, and build the
iPTM-vs-LIS scatter Philipp asked for.

Why this exists. Until now the AlphaFold results existed only as five ipTM
numbers hand-copied into NOTEBOOK.md, with no output files saved. Nothing was
reproducible, and the MqnC-MqnE complex (which the wet lab is to be built on)
had never been verified. This script recomputes everything from the raw
AlphaFold Server downloads.

LIS (Local Interaction Score), following Kim et al.: the inter-chain block of
the PAE matrix is inverted and scaled to [0,1] with a 12 A cutoff,
    s = (cutoff - PAE) / cutoff   for PAE < cutoff, else 0
and averaged over the interacting positions. Computed in both chain directions
and averaged. LIA (Local Interaction Area) is the count of residue pairs under
the cutoff. LIS is more sensitive than ipTM for small or transient interfaces,
which is why ipTM alone is considered weak.

Thresholds used (Philipp): iPTM > 0.6 AND LIS > 0.4 = possible real interaction.
"""
import json
import os
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(ROOT, "data", "af3_outputs", "extracted")
OUT = os.path.join(ROOT, "audit", "af3")
PAE_CUTOFF = 12.0

# Values transcribed into NOTEBOOK.md, to be verified against the raw files.
CLAIMED = {
    "mqne_x_mqnc_metabolon": 0.87,
    "verify_mqne_mqnc_actinobacteria": 0.87,
    "control_mqne_x_biob_unrelated_radicalsam": 0.18,
    "posctrl_barnase_barstar": 0.92,
    "negctrl_mqne_gfp_streptomyces_yangpuensis": 0.13,
    "mqne_duf4229_streptomyces_yangpuensis": 0.12,
    "mqne_duf4229_micromonospora_thermarum": 0.16,
    "mqne_duf4229_kitasatospora_papulosa": 0.19,
    "mqne_duf4229sol_streptomyces_yangpuensis": 0.19,
    "mqne_duf4229sol_micromonospora_thermarum": 0.20,
    "mqne_duf4229sol_kitasatospora_papulosa": 0.19,
}


def compute_lis(full_data_path):
    """Return (LIS, LIA, n_tokens) from an AF3 full_data json."""
    with open(full_data_path) as fh:
        d = json.load(fh)
    pae = d.get("pae")
    chains = d.get("token_chain_ids")
    if not pae or not chains:
        return None, None, None
    n = len(chains)
    uniq = []
    for c in chains:
        if c not in uniq:
            uniq.append(c)
    if len(uniq) < 2:
        return None, None, n
    a, b = uniq[0], uniq[1]
    idx_a = [i for i, c in enumerate(chains) if c == a]
    idx_b = [i for i, c in enumerate(chains) if c == b]

    tot, cnt = 0.0, 0
    for i in idx_a:                       # direction A -> B
        row = pae[i]
        for j in idx_b:
            v = row[j]
            if v < PAE_CUTOFF:
                tot += (PAE_CUTOFF - v) / PAE_CUTOFF
                cnt += 1
    for i in idx_b:                       # direction B -> A
        row = pae[i]
        for j in idx_a:
            v = row[j]
            if v < PAE_CUTOFF:
                tot += (PAE_CUTOFF - v) / PAE_CUTOFF
                cnt += 1
    lis = (tot / cnt) if cnt else 0.0
    return lis, cnt, n


def main():
    os.makedirs(OUT, exist_ok=True)
    if not os.path.isdir(EXT):
        sys.exit(f"no extracted AF3 outputs at {EXT}")

    rows = []
    for job in sorted(os.listdir(EXT)):
        jdir = os.path.join(EXT, job)
        if not os.path.isdir(jdir):
            continue
        key = job.replace("fold_", "")
        iptms, ptms, lis_vals, lia_vals = [], [], [], []
        for f in sorted(os.listdir(jdir)):
            if "summary_confidences" in f:
                with open(os.path.join(jdir, f)) as fh:
                    s = json.load(fh)
                if s.get("iptm") is not None:
                    iptms.append(float(s["iptm"]))
                if s.get("ptm") is not None:
                    ptms.append(float(s["ptm"]))
            elif "full_data" in f:
                lis, lia, _ = compute_lis(os.path.join(jdir, f))
                if lis is not None:
                    lis_vals.append(lis)
                    lia_vals.append(lia)
        if not iptms:
            continue
        rows.append({
            "job": key,
            "iptm_best": max(iptms), "iptm_mean": sum(iptms) / len(iptms),
            "ptm_best": max(ptms) if ptms else None,
            "lis_best": max(lis_vals) if lis_vals else None,
            "lis_mean": (sum(lis_vals) / len(lis_vals)) if lis_vals else None,
            "lia_best": max(lia_vals) if lia_vals else None,
            "n_seeds": len(iptms),
        })

    # ---------- verification against the notebook ----------
    print("=" * 92)
    print("PART 1: VERIFY THE HAND-TRANSCRIBED ipTM VALUES")
    print("=" * 92)
    print(f"  {'job':<46}{'claimed':<10}{'best':<8}{'mean':<8}{'verdict'}")
    mism = 0
    for r in rows:
        claim = None
        for k, v in CLAIMED.items():
            if k in r["job"]:
                claim = v
                break
        if claim is None:
            print(f"  {r['job'][:45]:<46}{'-':<10}{r['iptm_best']:<8.2f}{r['iptm_mean']:<8.2f}(not in notebook)")
            continue
        ok = abs(r["iptm_best"] - claim) <= 0.02
        if not ok:
            mism += 1
        print(f"  {r['job'][:45]:<46}{claim:<10.2f}{r['iptm_best']:<8.2f}{r['iptm_mean']:<8.2f}"
              f"{'MATCH' if ok else '*** MISMATCH ***'}")
    print(f"\n  mismatches: {mism}")

    # ---------- LIS ----------
    print("\n" + "=" * 92)
    print("PART 2: LIS (Local Interaction Score), PAE cutoff 12 A")
    print("Philipp's criterion: iPTM > 0.6 AND LIS > 0.4")
    print("=" * 92)
    print(f"  {'job':<46}{'iPTM':<8}{'LIS':<8}{'LIA':<8}{'call'}")
    for r in sorted(rows, key=lambda x: -(x["lis_best"] or 0)):
        call = "PASS" if (r["iptm_best"] > 0.6 and (r["lis_best"] or 0) > 0.4) else "fail"
        print(f"  {r['job'][:45]:<46}{r['iptm_best']:<8.2f}"
              f"{(r['lis_best'] or 0):<8.3f}{(r['lia_best'] or 0):<8}{call}")

    # ---------- null distribution ----------
    print("\n" + "=" * 92)
    print("PART 3: NULL DISTRIBUTION (non-partner pairs)")
    print("=" * 92)
    null_keys = ("negctrl", "control_mqne_x_biob", "mqne_x_reductase",
                 "mqne_x_other_duf", "duf4229")
    nulls = [r for r in rows if any(k in r["job"] for k in null_keys)]
    if nulls:
        ni = [r["iptm_best"] for r in nulls]
        nl = [(r["lis_best"] or 0) for r in nulls]
        print(f"  n non-partner pairs = {len(nulls)}")
        print(f"  iPTM range {min(ni):.2f} to {max(ni):.2f}")
        print(f"  LIS  range {min(nl):.3f} to {max(nl):.3f}")
        met = [r for r in rows if "mqnc" in r["job"]]
        for r in met:
            print(f"\n  {r['job']}: iPTM {r['iptm_best']:.2f}, LIS {(r['lis_best'] or 0):.3f}")
            print(f"    exceeds every non-partner on iPTM: {r['iptm_best'] > max(ni)}")
            print(f"    exceeds every non-partner on LIS : {(r['lis_best'] or 0) > max(nl)}")

    with open(os.path.join(OUT, "af3_audit.tsv"), "w") as fh:
        fh.write("job\tiptm_best\tiptm_mean\tptm_best\tlis_best\tlis_mean\tlia_best\tn_seeds\n")
        for r in rows:
            fh.write(f"{r['job']}\t{r['iptm_best']:.3f}\t{r['iptm_mean']:.3f}\t"
                     f"{r['ptm_best']}\t{r['lis_best']:.4f}\t{r['lis_mean']:.4f}\t"
                     f"{r['lia_best']}\t{r['n_seeds']}\n")
    print(f"\n[out] {os.path.join(OUT, 'af3_audit.tsv')}")
    return rows


if __name__ == "__main__":
    main()
