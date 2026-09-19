#!/usr/bin/env python3
import csv, os
import numpy as np
import matplotlib.pyplot as plt

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
rows = list(csv.DictReader(open("results/coevo_panel.tsv"), delimiter="\t"))
lab = [r["label"].replace(" (control)", "").replace(" (candidate)", "") for r in rows]
raw = [float(r["raw_r"]) for r in rows]
par = [float(r["partial_r"]) for r in rows]
x = np.arange(len(lab))
w = 0.38

fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(x - w / 2, raw, w, color="orange", edgecolor="black", alpha=0.5, label="Raw")
ax.bar(x + w / 2, par, w, color="cyan", edgecolor="black", alpha=0.5,
       label="Ancestry removed")
ax.set_xticks(x)
ax.set_xticklabels(lab, rotation=20, ha="right")
ax.set_ylabel("Mantel correlation")
ax.set_title("Co-evolution: Candidate vs Control Families")
ax.legend()

fig.savefig("latex/figures/fig_coevolution.png", dpi=300, bbox_inches="tight")
print("wrote latex/figures/fig_coevolution.png")
