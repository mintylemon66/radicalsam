#!/usr/bin/env python3
"""Figure 1. DUF4229 co-occurrence across radical SAM families, Actinomycetota only.
Input : results/duf4229_cross_family_specificity.tsv
Output: latex/figures/fig_specificity.png
"""
import csv, os
import matplotlib.pyplot as plt

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
rows = [r for r in csv.DictReader(
    (l for l in open("results/duf4229_cross_family_specificity.tsv") if not l.startswith("#")),
    delimiter="\t")]

names = {"MqnE": "MqnE", "MqnC": "MqnC", "background_radicalSAM": "background\nradical SAM",
         "BioB_biotin_synthase": "BioB", "LipA_lipoyl_synthase": "LipA"}
lab = [names.get(r["family"], r["family"]) for r in rows]
pos = [int(r["actino_DUF4229pos"]) for r in rows]
den = [int(r["actino_denom"]) for r in rows]
rate = [100.0 * p / d for p, d in zip(pos, den)]

fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(lab, rate, color="blue", edgecolor="black")
for i, (p, d, r) in enumerate(zip(pos, den, rate)):
    ax.text(i, r + 2, "%d/%d" % (p, d), ha="center")
ax.set_ylabel("Genomes with DUF4229 in window (%)")
ax.set_ylim(0, 100)
ax.set_title("DUF4229 Co-occurrence by Anchor Enzyme")
fig.savefig("latex/figures/fig_specificity.png", dpi=300, bbox_inches="tight")
print("wrote latex/figures/fig_specificity.png")
