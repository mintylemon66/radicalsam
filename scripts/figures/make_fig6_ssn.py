#!/usr/bin/env python3
"""Figure. Sequence similarity network of futalosine radical SAM enzymes.
Inputs : results/ssn_nodes.tsv, results/ssn_edges.tsv
Output : latex/figures/fig_ssn.png
Node and edge tables are the delivered network at a bit-score threshold of 380.
"""
import csv, os
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
nodes = list(csv.DictReader(open("results/ssn_nodes.tsv"), delimiter="\t"))
edges = list(csv.DictReader(open("results/ssn_edges.tsv"), delimiter="\t"))

G = nx.Graph()
for r in nodes:
    G.add_node(r["id"], sub=r["subfamily"], n=int(r["n_members"]))
for e in edges:
    G.add_edge(e["source"], e["target"], weight=float(e["bitscore"]))

colors = {"Actino MqnE": "navy", "Actino MqnC": "lightblue",
          "Bacillota MqnE": "orange", "Bacillota MqnC": "green"}
partner = {"Actino MqnE": "DUF4229", "Bacillota MqnE": "DUF1450",
           "Bacillota MqnC": "YtxC", "Actino MqnC": "DUF4229 via MqnE"}

pos = nx.spring_layout(G, seed=5, k=0.5, iterations=400, weight="weight")
fig, ax = plt.subplots(figsize=(7, 6))
nx.draw_networkx_edges(G, pos, alpha=0.15, width=0.5, edge_color="grey", ax=ax)
nx.draw_networkx_nodes(
    G, pos,
    node_color=[colors.get(G.nodes[a]["sub"], "grey") for a in G],
    node_size=[min(180, 40 + 3 * G.nodes[a]["n"]) for a in G],
    edgecolors="black", linewidths=0.4, ax=ax)
leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
              markeredgecolor="black", markersize=9, label=f"{k} ({partner[k]})")
       for k, c in colors.items()]
ax.legend(handles=leg, loc="upper left", fontsize=9)
ax.set_title("Sequence similarity network of futalosine radical SAM enzymes")
ax.axis("off")
fig.savefig("latex/figures/fig_ssn.png", dpi=300, bbox_inches="tight")
print("wrote latex/figures/fig_ssn.png (%d nodes, %d edges)" % (G.number_of_nodes(), G.number_of_edges()))
