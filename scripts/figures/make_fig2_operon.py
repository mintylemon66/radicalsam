#!/usr/bin/env python3
import json, os, collections
import numpy as np
import matplotlib.pyplot as plt

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
d = json.load(open("data/operon_pairs_566.json"))
gaps = [abs(g) for g in d["gaps"]]
offs = [o for o in d["offs"] if abs(o) <= 8]
same, tot = d["same"], d["same"] + d["diff"]
med = int(np.median(gaps))

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

c = collections.Counter(offs)
ks = sorted(c)
axes[0].bar(ks, [c[k] for k in ks], color="blue", edgecolor="black")
axes[0].set_xlabel("Gene offset from futalosine enzyme")
axes[0].set_ylabel("Frequency")
axes[0].set_title("Gene Offset Distribution")

axes[1].hist(gaps, bins=np.arange(0, 4001, 200), color="green", edgecolor="black")
axes[1].axvline(med, color="r", linestyle="--")
axes[1].set_xlabel("Intergenic gap (bp)")
axes[1].set_ylabel("Frequency")
axes[1].set_title("Intergenic Gap Distribution (median %d bp)" % med)

fig.savefig("latex/figures/fig_operon.png", dpi=300, bbox_inches="tight")
print("wrote latex/figures/fig_operon.png  (%d pairs, %d same-strand, median %d bp)"
      % (tot, same, med))
