#!/usr/bin/env python3
"""EvoWeaver labels for the fair test: DUF4229 against the genes that also sit near MqnE.

    python3 _scripts/ensemble/build_neighborhood_labels.py      # from the project folder

The benchmark scored DUF4229 against random Pfam families, which any gene picked for its
closeness to MqnE would beat. Here the candidates are the 145 Pfam families that recur within
10 genes of MqnE (neighbor_candidates.tsv), so every candidate shares that closeness. The two
known futalosine genes, MqnA/MqnD and MqnB, are added as reference points. Labels are built
exactly as for the benchmark, on the MqnE contigs of the futalosine system.

Writes _data/ensemble/labels_futalosine_neighborhood.tsv and neighborhood_candidates.tsv.
"""

import csv
import sys

sys.path.insert(0, "_scripts/ensemble")
sys.path.insert(0, "_scripts/evoweaver/genomewide")
from build_training_pairs import build  # noqa: E402

OUT = "_data/ensemble/"
NEIGHBORS = "_data/evoweaver/genomewide/neighbor_candidates.tsv"
KNOWN = {"fam_PF02621": "MqnA/MqnD", "fam_PF01048": "MqnB"}
DUF4229 = "fam_PF14012"


def log(msg):
    print(msg, flush=True)


def neighbor_families(path):
    with open(path) as handle:
        return {row["acc"]: row["name"] for row in csv.DictReader(handle, delimiter="\t")}


def main():
    neighbors = neighbor_families(NEIGHBORS)
    candidates = list(KNOWN) + [DUF4229] + sorted(neighbors)
    labels = build("futalosine", candidates)
    assert labels.get("ANCHOR"), "MqnE has no labels, nothing can be scored against it"

    with open(OUT + "labels_futalosine_neighborhood.tsv", "w") as handle:
        handle.write("family\tlabel\n")
        for family, labs in labels.items():
            for lab in labs:
                handle.write(f"{family}\t{lab}\n")

    placed = [c for c in candidates if labels.get(c)]
    with open(OUT + "neighborhood_candidates.tsv", "w") as handle:
        handle.write("candidate\tname\tkind\n")
        for c in placed:
            kind = "known" if c in KNOWN else "DUF4229" if c == DUF4229 else "neighborhood"
            handle.write(f"{c}\t{KNOWN.get(c) or neighbors.get(c, 'DUF4229')}\t{kind}\n")
    log(f"MqnE anchor: {len(labels['ANCHOR'])} labels")
    log(f"{len(placed)} of {len(candidates)} candidates placed "
        f"({sum(c in neighbors for c in placed)} neighborhood families)")


if __name__ == "__main__":
    main()
