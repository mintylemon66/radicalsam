#!/usr/bin/env python3
"""Labels and the pair table for the neighborhood benchmark.

    python3 _scripts/ensemble/build_neighbor_pairs.py      # after pick_system_neighbors.py

Positives are each pathway's genes, exactly as in build_training_pairs.py. Negatives are that
system's neighborhood families from system_neighbors.tsv. Labels are built the same way, on the
contigs that carry each anchor.

Writes _data/ensemble/labels_nb_<system>.tsv and pair_labels_nb.tsv.
"""

import csv
import os
import sys

sys.path.insert(0, "_scripts/ensemble")
sys.path.insert(0, "_scripts/evoweaver/genomewide")
from build_training_pairs import build, pathway  # noqa: E402
from systems import TRAIN  # noqa: E402

OUT = "_data/ensemble/"
GENES = "_data/evoweaver/genomewide/genes/"


def log(msg):
    print(msg, flush=True)


def main():
    have = {f[:-4] for f in os.listdir(GENES) if f.endswith(".tsv")}
    neighbors = {}
    with open(OUT + "system_neighbors.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            neighbors.setdefault(row["system"], []).append(row["acc"])

    rows = []
    for system in TRAIN:
        positives = sorted({f"fam_{pf}" for pf in pathway.get(system, {}).values()} & have)
        taken = {p[4:] for p in positives}
        negatives = [acc for acc in neighbors.get(system, []) if acc in have and acc not in taken]
        labels = build(system, positives + negatives)
        with open(f"{OUT}labels_nb_{system}.tsv", "w") as handle:
            handle.write("family\tlabel\n")
            for family, labs in labels.items():
                for lab in labs:
                    handle.write(f"{family}\t{lab}\n")
        placed = [c for c in positives + negatives if labels.get(c)]
        rows += [(system, c, int(c in positives)) for c in placed]
        log(f"{system:14s} {sum(c in positives for c in placed)} positives, "
            f"{sum(c not in positives for c in placed)} neighborhood negatives")

    with open(OUT + "pair_labels_nb.tsv", "w") as handle:
        handle.write("system\tcandidate\tis_partner\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    log(f"\n{sum(r[2] for r in rows)} positives and {sum(1 - r[2] for r in rows)} negatives across {len(TRAIN)} systems")


if __name__ == "__main__":
    main()
