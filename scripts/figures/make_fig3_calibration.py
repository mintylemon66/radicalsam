#!/usr/bin/env python3
"""Figure 3. Cross-validation protocol calibration on signal-free data.
Input : results/cv_protocol_calibration.tsv
Output: latex/figures/fig_calibration.png
"""
import csv, os
import numpy as np
import matplotlib.pyplot as plt

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
rows = list(csv.DictReader(open("results/cv_protocol_calibration.tsv"), delimiter="\t"))
lab = [" ".join(r["protocol"].split()) for r in rows]
nm = [float(r["null_median"]) for r in rows]
fp = [100.0 * float(r["fpr_at_0.55"]) for r in rows]
x = np.arange(len(lab))

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

axes[0].bar(x, nm, color="blue", edgecolor="black")
axes[0].axhline(0.5, color="r", linestyle="--")
axes[0].set_xticks(x)
axes[0].set_xticklabels(lab, rotation=20, ha="right")
axes[0].set_ylabel("Null median AUC")
axes[0].set_ylim(0.4, 0.8)
axes[0].set_title("Null Median AUC on Signal-Free Data")

axes[1].bar(x, fp, color="green", edgecolor="black")
axes[1].set_xticks(x)
axes[1].set_xticklabels(lab, rotation=20, ha="right")
axes[1].set_ylabel("Replicates exceeding AUC 0.55 (%)")
axes[1].set_title("False Positive Rate at AUC 0.55")

fig.savefig("latex/figures/fig_calibration.png", dpi=300, bbox_inches="tight")
print("wrote latex/figures/fig_calibration.png")
