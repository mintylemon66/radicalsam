#!/usr/bin/env python3
"""Figure 3: where DUF4229 sits relative to its futalosine anchor enzyme.

    mlenv/bin/python _scripts/figures/make_fig2_operon.py

Reads the reconstructed neighborhoods (neighbors_seq.tsv), the PF14012 hits
from hmmsearch at the Pfam gathering threshold, and the MqnE/MqnC call that the
two-profile HMM classifier made for every anchor. The neighbor table lists 224
anchor and DUF4229 pairs twice, so each pair is kept once (342 pairs, not 566
rows). gene_index counts genes left to right along the contig, so offsets are
flipped for anchors on the minus strand; a positive offset always means
downstream of the anchor.

Writes latex/figures/fig_operon.png: gene offset on the left and the distance
between the two genes on the right, both split by anchor enzyme, and prints
the counts the paper quotes so they can be checked against the text.
"""

import argparse
import csv
import os
from collections import Counter
from statistics import median

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(os.path.dirname(HERE))
NEIGHBORS_TSV = os.path.join(PROJECT, "_audit", "snap", "neighbors_seq.tsv")
DUF4229_HITS = os.path.join(PROJECT, "_audit", "snap", "PF14012_vs_neighbors.domtbl")
ENZYME_CALLS = os.path.join(PROJECT, "_audit", "mqn_seq_classification.tsv")
OUT_PNG = os.path.join(PROJECT, "latex", "figures", "fig_operon.png")

MAX_DISTANCE = 8000  # bp, the operonic window used throughout the paper
MAX_OFFSET = 10      # genes
ENZYMES = ("MqnE", "MqnC")
COLORS = {"MqnE": "blue", "MqnC": "orange"}


def log(msg):
    print(msg, flush=True)


# protein ids that hit PF14012 (the search was run with --cut_ga)
def read_duf4229_hits(path):
    hits = set()
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                hits.add(line.split()[0])
    return hits


# anchor accession -> MqnE or MqnC
def read_enzyme_calls(path):
    with open(path) as handle:
        return {row["accession"]: row["seq_assignment"] for row in csv.DictReader(handle, delimiter="\t")}


def collect_pairs(neighbors_tsv, hits, calls):
    pairs = {}
    rows = 0
    with open(neighbors_tsv) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["group"] != "mqn" or row["n_protein_id"] not in hits:
                continue
            if not row["q_strand"] or not row["n_strand"]:
                continue
            try:
                distance = int(row["bp_gap"])
            except ValueError:
                continue
            offset = int(row["gene_index"] or 0)
            if distance > MAX_DISTANCE or abs(offset) > MAX_OFFSET:
                continue
            rows += 1

            # positive offset = downstream of the anchor, whichever strand it is on
            if row["q_strand"] == "-":
                offset = -offset

            # setdefault keeps the first copy of a pair and ignores the duplicate row
            pairs.setdefault((row["query_acc"], row["n_protein_id"]), {
                "enzyme": calls[row["query_acc"]],
                "offset": offset,
                "distance": distance,
                "same_strand": row["q_strand"] == row["n_strand"],
            })
    return list(pairs.values()), rows


def summarize(pairs, rows):
    log(f"{rows} rows in the window, {len(pairs)} distinct anchor and DUF4229 pairs")
    groups = [("all", pairs)] + [(enzyme, [p for p in pairs if p["enzyme"] == enzyme]) for enzyme in ENZYMES]
    for name, group in groups:
        offsets = Counter(p["offset"] for p in group)
        downstream = sum(count for offset, count in offsets.items() if offset > 0)
        same = sum(p["same_strand"] for p in group)
        log(f"  {name}: {len(group)} pairs, {same} same strand, {offsets[2]} two genes downstream, "
            f"{downstream} downstream, median distance {median(p['distance'] for p in group):.0f} bp")


def plot(pairs, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # left: one bar per gene offset, MqnE and MqnC anchors stacked
    offsets = np.arange(-MAX_OFFSET, MAX_OFFSET + 1)
    bottom = np.zeros(len(offsets))
    for enzyme in ENZYMES:
        counts = Counter(p["offset"] for p in pairs if p["enzyme"] == enzyme)
        heights = np.array([counts[offset] for offset in offsets])
        axes[0].bar(offsets, heights, bottom=bottom, color=COLORS[enzyme], edgecolor="black",
                    label=f"{enzyme} anchor")
        bottom += heights
    axes[0].set_xticks(offsets[::2])
    axes[0].set_ylim(0, bottom.max() * 1.1)  # headroom, or the tallest bar touches the frame
    axes[0].set(title="Gene Offset Distribution", xlabel="Gene offset from anchor (positive = downstream)",
                ylabel="Gene pairs")
    axes[0].legend(frameon=False, loc="upper left")

    # right: distance between the two genes, stacked the same way
    distances = [[p["distance"] for p in pairs if p["enzyme"] == enzyme] for enzyme in ENZYMES]
    axes[1].hist(distances, bins=np.arange(0, MAX_DISTANCE + 1, 250), stacked=True,
                 color=[COLORS[enzyme] for enzyme in ENZYMES], edgecolor="black",
                 label=[f"{enzyme} anchor" for enzyme in ENZYMES])
    middle = median(p["distance"] for p in pairs)
    axes[1].axvline(middle, color="black", linestyle="--", label=f"median {middle:.0f} bp")
    axes[1].set(title="Distance Between the Two Genes", xlabel="Distance (bp)", ylabel="Gene pairs")
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--neighbors", default=NEIGHBORS_TSV)
    parser.add_argument("--hits", default=DUF4229_HITS)
    parser.add_argument("--enzyme-calls", default=ENZYME_CALLS)
    parser.add_argument("--out", default=OUT_PNG)
    args = parser.parse_args()

    hits = read_duf4229_hits(args.hits)
    calls = read_enzyme_calls(args.enzyme_calls)
    log(f"{len(hits)} PF14012 hits, {len(calls):,} classified anchors")
    pairs, rows = collect_pairs(args.neighbors, hits, calls)
    summarize(pairs, rows)
    plot(pairs, args.out)
    log(f"wrote {args.out}")


if __name__ == "__main__":
    main()
