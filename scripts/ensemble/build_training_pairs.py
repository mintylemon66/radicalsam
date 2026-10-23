"""Write EvoWeaver labels for every system, and the table saying which candidates are partners.

For each system the anchor enzyme is scored against a set of candidate families. A candidate is a
positive when it is one of that pathway's genes and a negative when it comes from the random Pfam pool.
Labels are genome_contig_strand_index, the form EvoWeaver needs for its gene-organization methods as
well as its profiling ones.

Only genomes carrying the anchor contribute positions, since a family in a genome with no anchor has no
distance to measure. Copies are ordered nearest-the-anchor first, because GeneDistance reads only the
first label of each genome and OrientationMI the first of each genome and contig.
"""

import gzip
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, "_scripts/ensemble")
sys.path.insert(0, "_scripts/evoweaver/genomewide")
from fetch_families import read_genes  # noqa: E402
from systems import HELDOUT, TRAIN  # noqa: E402

OUT = "_data/ensemble/"
kept = set(open("_data/evoweaver/genomewide/proteomes_kept.txt").read().split())
pathway = json.load(open(OUT + "pathway_pfams.json"))
contig_id = {}


def number(name):
    return contig_id.setdefault(name, len(contig_id) + 1)


def positions(system):
    """{(contig, lowercased locus tag): (index, strand)} for the contigs carrying this anchor."""
    where = {}
    for line in gzip.open(f"{OUT}contigs/{system}.tsv.gz", "rt"):
        contig, idx, start, end, strand, tags, pid = line.rstrip("\n").split("\t")
        for t in tags.split(","):
            if t:
                where[(contig, t.casefold())] = (int(idx), strand)
    return where


def place(family, where, anchor_contigs, anchor_accs, is_anchor=False):
    """genome -> list of (contig id, strand, index) for every copy of this family.

    A candidate that is itself one of the anchor proteins is dropped, since it would sit at distance
    zero from itself. That exclusion is skipped for the anchor family, whose every protein is in that set.
    """
    out = defaultdict(list)
    for r in read_genes(family):
        up = r["proteome"]
        if up not in kept or up not in anchor_contigs:
            continue
        if not is_anchor and r["acc"] in anchor_accs.get(up, ()):
            continue
        hits = {(number(c), int(where[(c, t.casefold())][1]), where[(c, t.casefold())][0])
                for c in r["contigs"] if c in anchor_contigs[up]
                for t in r["tags"] if (c, t.casefold()) in where}
        # a copy on none of the anchor's contigs gets a contig of its own, so it is never a neighbor
        out[up] += sorted(hits) or [(number("off:" + r["acc"]), 0, 1)]
    return out


def build(system, candidates):
    where = positions(system)
    anchor_contigs, anchor_accs = defaultdict(set), defaultdict(set)
    for r in read_genes("anchor_" + system):
        if r["proteome"] in kept:
            anchor_contigs[r["proteome"]].update(r["contigs"])
            anchor_accs[r["proteome"]].add(r["acc"])

    placed = {"ANCHOR": place("anchor_" + system, where, anchor_contigs, anchor_accs, is_anchor=True)}
    for fam in candidates:
        placed[fam] = place(fam, where, anchor_contigs, anchor_accs)

    first = {up: v[0] for up, v in placed["ANCHOR"].items() if v}
    labels = defaultdict(list)
    for name, genomes in placed.items():
        for up, copies in genomes.items():
            a = first.get(up)
            if a:
                copies = sorted(set(copies), key=lambda x: (x[0] != a[0], abs(x[2] - a[2])))
            for contig, strand, idx in copies:
                labels[name].append(f"{up}_{contig}_{strand}_{idx}")
    return labels


if __name__ == "__main__":
    have = {f[:-4] for f in os.listdir("_data/evoweaver/genomewide/genes") if f.endswith(".tsv")}
    negatives = sorted(n for n in have if n.startswith("fam_")
                       and n[4:] not in {pf for g in pathway.values() for pf in g.values()})
    rows = []
    for system in list(TRAIN) + list(HELDOUT):
        positives = sorted({f"fam_{pf}" for pf in pathway.get(system, {}).values()} & have)
        # a system's named partner can sit outside its pathway Pfam set (PflB is TIGR01255 in
        # systems.py and PF02901 in the pathway map), and it must never be scored as a negative
        own = "fam_" + {**TRAIN, **HELDOUT}[system]["partner"].split("-")[1]
        candidates = positives + [n for n in negatives if n not in positives and n != own]
        labels = build(system, candidates)
        with open(f"{OUT}labels_{system}.tsv", "w") as f:
            f.write("family\tlabel\n")
            for fam, labs in labels.items():
                for l in labs:
                    f.write(f"{fam}\t{l}\n")
        assert labels.get("ANCHOR"), f"{system}: anchor has no labels, nothing can be scored against it"
        for c in candidates:
            if labels.get(c):
                rows.append((system, c, int(c in positives)))
        print(f"{system:14s} anchor {len(labels['ANCHOR']):6d} labels | "
              f"{sum(1 for c in positives if labels.get(c))} positives, "
              f"{sum(1 for c in candidates if c not in positives and labels.get(c))} negatives", flush=True)
    with open(OUT + "pair_labels.tsv", "w") as f:
        f.write("system\tcandidate\tis_partner\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    train = [r for r in rows if r[0] != "futalosine"]
    print(f"\ntraining set: {sum(r[2] for r in train)} positives, {sum(1 - r[2] for r in train)} negatives")
