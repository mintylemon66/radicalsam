#!/usr/bin/env python3
"""Figure. AlphaFold 3 interface scores for every co-folded pair.
Input : results/af3_interface_scores.tsv
Output: latex/figures/fig_af3.png
"""
import csv, os
import matplotlib.pyplot as plt

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
rows = list(csv.DictReader(open("results/af3_interface_scores.tsv"), delimiter="\t"))

def group(r):
    if "positive" in r["role"]: return "positive control"
    if "candidate complex" in r["role"]: return "MqnE + MqnC"
    return "all other pairs"

colors = {"positive control": "green", "MqnE + MqnC": "blue", "all other pairs": "grey"}
fig, ax = plt.subplots(figsize=(6, 4.5))
for g in colors:
    xs = [float(r["iptm"]) for r in rows if group(r) == g]
    ys = [float(r["lis"])  for r in rows if group(r) == g]
    ax.scatter(xs, ys, c=colors[g], edgecolors="black", s=70, label=g)
ax.axvline(0.6, color="r", linestyle="--")
ax.axhline(0.4, color="r", linestyle="--")
ax.set_xlabel("Interface predicted TM-score (ipTM)")
ax.set_ylabel("Local Interaction Score (LIS)")
ax.set_title("AlphaFold 3 Interface Scores")
ax.legend(fontsize=9)
fig.savefig("latex/figures/fig_af3.png", dpi=300, bbox_inches="tight")
print("wrote latex/figures/fig_af3.png (%d pairs)" % len(rows))
