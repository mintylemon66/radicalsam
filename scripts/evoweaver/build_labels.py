"""Turn the downloaded families into the labels EvoWeaver reads.

    python3 _scripts/evoweaver/genomewide/build_labels.py      # after fetch_contig_tables.py

Step 7. EvoWeaver takes one label per gene copy, genome_contig_strand_index, and reads both
co-occurrence and gene organization out of them. Three choices here decide what the scores mean:

A protein that is itself an MqnE is dropped from every other family, since it would sit at distance
zero from MqnE and manufacture a perfect organization score. A copy on none of the genome's MqnE
contigs is given a contig of its own, so it counts for presence but never as a neighbor. Labels are
written best-first, because GeneDistance reads only the first label of each genome and OrientationMI
the first of each genome and contig.

The random families are filtered to those within 0.5-2x DUF4229's genome count, so prevalence alone
cannot separate it from the null. The neighborhood families are kept as they are.

Writes labels.tsv.
"""

import gzip
from collections import defaultdict

from fetch_families import proteomes_with, read_genes

OUT = "_data/evoweaver/genomewide/"
TEST = ["MqnE", "MqnC", "DUF4229", "MqnAD", "Fer4"]

kept = set(open(OUT + "proteomes_kept.txt").read().split())
target = len(proteomes_with("DUF4229") & kept)


def comparison_set(path, prefix, matched):
    # families from a candidates file; with matched=True only those whose genome count in the kept
    # proteomes is within 0.5-2x DUF4229's
    out = {}
    for line in open(OUT + path).read().splitlines()[1:]:
        acc, name = line.split("\t")[:2]
        if not matched or 0.5 * target <= len(proteomes_with(acc) & kept) <= 2 * target:
            out[acc] = f"{prefix}_{acc}_{name}"
    return out


nulls = comparison_set("null_candidates.tsv", "null", matched=True)       # random Pfam families
nbrs = comparison_set("neighbor_candidates.tsv", "nbr", matched=False)    # families that recur within 10 genes of MqnE
families = {f: f for f in TEST} | nulls | nbrs

# every proteome's MqnE proteins and contigs, and where each locus tag sits on those contigs
mqne_contigs, mqne_accs = defaultdict(set), defaultdict(set)
for r in read_genes("MqnE"):
    if r["proteome"] in kept:
        mqne_contigs[r["proteome"]].update(r["contigs"])
        mqne_accs[r["proteome"]].add(r["acc"])
where = {}
for line in gzip.open(OUT + "mqne_contig_cds.tsv.gz", "rt"):
    contig, idx, start, end, strand, tags, pid = line.rstrip("\n").split("\t")
    for t in tags.split(","):
        if t:
            where[(contig, t.casefold())] = (int(idx), strand)    # tag case differs between UniProt and GenBank

# EvoWeaver's MoransI reads the contig field as a number, so every contig gets an integer id
contig_id = {}


def number(name):
    return contig_id.setdefault(name, len(contig_id) + 1)


found = defaultdict(lambda: defaultdict(list))
for fam, name in families.items():
    for r in read_genes(fam):
        up = r["proteome"]
        if up not in kept or (fam != "MqnE" and r["acc"] in mqne_accs[up]):
            continue          # a protein that is itself an MqnE would sit at distance 0 from MqnE
        hits = {(number(c), int(where[(c, t.casefold())][1]), where[(c, t.casefold())][0])
                for c in r["contigs"] if c in mqne_contigs.get(up, ())
                for t in r["tags"] if (c, t.casefold()) in where}
        # a copy off every MqnE contig gets a contig of its own, so it never counts as an MqnE neighbor
        found[name][up] += sorted(hits) or [(number("protein:" + r["acc"]), 0, 1)]

# GeneDistance reads the first label of each genome, OrientationMI and MoransI the first of each
# genome and contig, so labels go best-first: MqnE's first copy, then each family's copy closest to it
with open(OUT + "labels.tsv", "w") as out:
    out.write("family\tlabel\n")
    for name, genomes in found.items():
        for up, copies in genomes.items():
            first = found["MqnE"][up][0] if up in found["MqnE"] else None
            if first:
                copies = sorted(set(copies), key=lambda x: (x[0] != first[0], abs(x[2] - first[2])))
            for contig, strand, idx in copies:
                out.write(f"{name}\t{up}_{contig}_{strand}_{idx}\n")

placed = lambda up, copies: any(c == found["MqnE"][up][0][0] for c, _, _ in copies) if up in found["MqnE"] else False
print(f"{len(kept)} proteomes; MqnE in {len(found['MqnE'])}, DUF4229 in {target}; "
      f"{len(nulls)} random families in the 0.5-2x window, {len(nbrs)} MqnE-neighborhood families")
for fam in TEST:
    g = found[fam]
    both = [up for up in g if up in found["MqnE"]]
    print(f"  {fam}: {len(g)} genomes, {len(both)} with MqnE, {sum(placed(up, g[up]) for up in both)} on MqnE's contig")
