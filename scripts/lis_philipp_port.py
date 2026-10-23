#!/usr/bin/env python3
"""
LIS computed with Philipp's implementation, ported from his AF3-pipeline.

SOURCE / ATTRIBUTION
    https://github.com/FilyCode/AF3-pipeline
    04_post_processing_results.ipynb
    Functions `transform_pae_matrix` and `calculate_mean_lis` are reproduced
    here essentially verbatim (numpy logic unchanged). The contact-map step
    (`cLIS`) follows his `calculate_contact_map`: CB atoms, CA for glycine,
    8 A threshold. The underlying metric is LIS from Kim et al.

    Local additions are only: reading AlphaFold *Server* output layout
    (his pipeline reads local AF3 install layout), and the comparison table
    against our independent reimplementation.

WHY THIS EXISTS
    Our first LIS numbers came from an independent reimplementation written
    from the published method description, not from his code. This re-runs the
    calculation with his actual implementation so the reported values are the
    established ones, and our version becomes a reproducibility check rather
    than the primary result.

THREE DIFFERENCES FOUND vs our reimplementation
    1. Seed aggregation: he takes the MEAN of LIS across the 5 seeds
       (`group['LIS'].mean()`); we reported the MAX.
    2. He also computes cLIS, a contact-map-filtered LIS, which we did not.
    3. His notebook's own filter is iPTM > 0.4 and LIS > 0.25, whereas the
       thresholds he quoted to us were iPTM > 0.6 and LIS > 0.4. We report
       against both.
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(ROOT, "data", "af3_outputs", "extracted")
OUT = os.path.join(ROOT, "audit", "af3")
PAE_CUTOFF = 12
CONTACT_THRESHOLD = 8.0


# ---------------------------------------------------------------- his code --
def transform_pae_matrix(pae_matrix, pae_cutoff):
    transformed_pae = np.zeros_like(pae_matrix)
    within_cutoff = pae_matrix < pae_cutoff
    transformed_pae[within_cutoff] = 1 - (pae_matrix[within_cutoff] / pae_cutoff)
    return transformed_pae


def calculate_mean_lis(transformed_pae, subunit_number):
    cum_lengths = np.cumsum(subunit_number)
    start_indices = np.concatenate(([0], cum_lengths[:-1]))
    mean_lis_matrix = np.zeros((len(subunit_number), len(subunit_number)))
    for i in range(len(subunit_number)):
        for j in range(len(subunit_number)):
            start_i, end_i = start_indices[i], cum_lengths[i]
            start_j, end_j = start_indices[j], cum_lengths[j]
            submatrix = transformed_pae[start_i:end_i, start_j:end_j]
            relevant_elements = submatrix[submatrix > 0]
            mean_lis_matrix[i, j] = relevant_elements.mean() if relevant_elements.size else 0.0
    return mean_lis_matrix
# ------------------------------------------------------------ end his code --


def contact_map_from_cif(cif_path, n_tokens):
    """CB atoms (CA for GLY), 8 A threshold. Follows his calculate_contact_map."""
    coords, chains = [], []
    try:
        with open(cif_path) as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                p = line.split()
                if len(p) < 16:      # AF3 Server CIF has 18 fields
                    continue
                atom, res = p[3], p[5]
                if not (atom == "CB" or (res == "GLY" and atom == "CA")):
                    continue
                try:
                    coords.append([float(p[10]), float(p[11]), float(p[12])])
                except (ValueError, IndexError):
                    continue
                chains.append(p[6])
    except OSError:
        return None, None
    if len(coords) < 2:
        return None, None
    c = np.asarray(coords)
    d = np.linalg.norm(c[:, None, :] - c[None, :, :], axis=-1)
    return (d < CONTACT_THRESHOLD).astype(int), chains


def subunit_lengths(chain_ids):
    lens, cur, last = [], 0, None
    for c in chain_ids:
        if last is None or c == last:
            cur += 1
        else:
            lens.append(cur)
            cur = 1
        last = c
    lens.append(cur)
    return lens


def main():
    if not os.path.isdir(EXT):
        sys.exit("no extracted AF3 outputs")
    os.makedirs(OUT, exist_ok=True)

    rows = []
    for job in sorted(os.listdir(EXT)):
        jdir = os.path.join(EXT, job)
        if not os.path.isdir(jdir):
            continue
        per_seed_lis, per_seed_clis, iptms = [], [], []
        for f in sorted(os.listdir(jdir)):
            path = os.path.join(jdir, f)
            if "summary_confidences" in f:
                s = json.load(open(path))
                if s.get("iptm") is not None:
                    iptms.append(float(s["iptm"]))
            elif "full_data" in f:
                d = json.load(open(path))
                pae = np.asarray(d["pae"], dtype=float)
                chains = d["token_chain_ids"]
                subs = subunit_lengths(chains)
                if len(subs) < 2:
                    continue
                t = transform_pae_matrix(pae, PAE_CUTOFF)
                m = calculate_mean_lis(t, subs)
                # interchain LIS = mean of the two off-diagonal blocks
                per_seed_lis.append((m[0, 1] + m[1, 0]) / 2.0)

                cif = path.replace("full_data", "model").replace(".json", ".cif")
                cmap, cch = contact_map_from_cif(cif, len(chains))
                if cmap is not None and cmap.shape[0] == pae.shape[0]:
                    tc = t * cmap
                    mc = calculate_mean_lis(tc, subs)
                    per_seed_clis.append((mc[0, 1] + mc[1, 0]) / 2.0)
        if not iptms:
            continue
        rows.append({
            "job": job.replace("fold_", ""),
            "iptm_mean": float(np.mean(iptms)),
            "iptm_max": float(np.max(iptms)),
            "lis_mean": float(np.mean(per_seed_lis)) if per_seed_lis else 0.0,
            "lis_max": float(np.max(per_seed_lis)) if per_seed_lis else 0.0,
            "clis_mean": float(np.mean(per_seed_clis)) if per_seed_clis else None,
            "n_seeds": len(iptms),
        })

    print("=" * 104)
    print("LIS via Philipp's implementation (AF3-pipeline 04_post_processing_results.ipynb)")
    print("His aggregation: MEAN across seeds.  cLIS = contact-map filtered (CB, 8 A)")
    print("=" * 104)
    print(f"  {'job':<44}{'iPTM(mn)':<10}{'LIS(mn)':<10}{'LIS(mx)':<10}{'cLIS':<10}"
          f"{'his 0.4/0.25':<14}{'quoted 0.6/0.4'}")
    for r in sorted(rows, key=lambda x: -x["lis_mean"]):
        c = f"{r['clis_mean']:.3f}" if r["clis_mean"] is not None else "n/a"
        his = "PASS" if (r["iptm_mean"] > 0.4 and r["lis_mean"] > 0.25) else "fail"
        quoted = "PASS" if (r["iptm_mean"] > 0.6 and r["lis_mean"] > 0.4) else "fail"
        print(f"  {r['job'][:43]:<44}{r['iptm_mean']:<10.3f}{r['lis_mean']:<10.3f}"
              f"{r['lis_max']:<10.3f}{c:<10}{his:<14}{quoted}")

    out = os.path.join(OUT, "lis_philipp_implementation.tsv")
    with open(out, "w") as fh:
        fh.write("job\tiptm_mean\tiptm_max\tlis_mean\tlis_max\tclis_mean\tn_seeds\n")
        for r in rows:
            c = "" if r["clis_mean"] is None else f"{r['clis_mean']:.4f}"
            fh.write(f"{r['job']}\t{r['iptm_mean']:.4f}\t{r['iptm_max']:.4f}\t"
                     f"{r['lis_mean']:.4f}\t{r['lis_max']:.4f}\t{c}\t{r['n_seeds']}\n")
    print(f"\n[out] {out}")

    # ---- reproducibility check against our independent reimplementation ----
    ours_path = os.path.join(OUT, "af3_audit.tsv")
    if os.path.exists(ours_path):
        import csv
        ours = {r["job"]: float(r["lis_best"])
                for r in csv.DictReader(open(ours_path), delimiter="\t")}
        print("\n" + "=" * 104)
        print("REPRODUCIBILITY CHECK: his implementation vs our independent reimplementation")
        print("(compare like with like: his LIS(max) against our lis_best)")
        print("=" * 104)
        print(f"  {'job':<44}{'his LIS(mx)':<14}{'ours':<14}{'abs diff'}")
        diffs = []
        for r in sorted(rows, key=lambda x: -x["lis_max"]):
            o = ours.get(r["job"])
            if o is None:
                continue
            d = abs(r["lis_max"] - o)
            diffs.append(d)
            print(f"  {r['job'][:43]:<44}{r['lis_max']:<14.4f}{o:<14.4f}{d:.4f}")
        if diffs:
            print(f"\n  max absolute difference: {max(diffs):.4f}")
            print(f"  mean absolute difference: {sum(diffs)/len(diffs):.4f}")


if __name__ == "__main__":
    main()
