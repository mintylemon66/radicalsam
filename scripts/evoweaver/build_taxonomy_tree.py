"""The species tree EvoWeaver needs, built from UniProt lineage strings.

    python3 _scripts/evoweaver/genomewide/build_taxonomy_tree.py      # after fetch_proteomes.py

Step 2. EvoWeaver's phylogenetic-structure and profiling methods need a tree over the genomes being
compared. Rather than infer one, each proteome is hung off its UniProt lineage, one unit of branch
length per taxonomic rank. A taxonomy tree has no branch lengths of its own, so sibling clades are
split into a binary tree with near-zero inner edges, which ReadDendrogram requires.

Proteomes below BUSCO completeness 95% are dropped, since a fragmentary assembly looks like a
genuine absence to a presence/absence method and would count as evidence against co-occurrence.

Writes taxonomy_tree.nwk and proteomes_kept.txt.
"""

import csv
import re
import sys

OUT = "_data/evoweaver/genomewide/"
MIN_BUSCO = 95.0
sys.setrecursionlimit(100000)

keep = {}
for r in csv.DictReader(open(OUT + "proteomes.tsv"), delimiter="\t"):
    m = re.match(r"C:([\d.]+)%", r["BUSCO"] or "")
    if m and float(m.group(1)) >= MIN_BUSCO:
        keep[r["Proteome Id"]] = r["Taxonomic lineage"].split(", ")

# taxonomy as a nested dict: each lineage name is a level, proteomes hang off their last rank
root = {}
for up, lineage in keep.items():
    node = root
    for name in lineage:
        node = node.setdefault(name, {})
    node.setdefault("__leaves__", []).append(up)


TAXON_EDGE = 1.0     # one unit per taxonomy rank
SPLIT_EDGE = 0.001   # inner edges made only to turn a multifurcation into a binary tree


def render(node):
    # returns (newick without its own edge, length of the edge above it)
    parts = [(up, TAXON_EDGE) for up in node.get("__leaves__", [])]
    for name in sorted(k for k in node if k != "__leaves__"):
        sub, length = render(node[name])
        parts.append((sub, length + TAXON_EDGE))
    if len(parts) == 1:
        return parts[0]
    return binary(parts)[0], 0.0


def binary(parts):
    if len(parts) == 1:
        return parts[0]
    mid = len(parts) // 2
    (a, la), (b, lb) = binary(parts[:mid]), binary(parts[mid:])
    return f"({a}:{la:g},{b}:{lb:g})", SPLIT_EDGE


tree = render(root)[0]
with open(OUT + "taxonomy_tree.nwk", "w") as f:
    f.write(tree + ";\n")   # no edge length on the root: ReadDendrogram would turn one into a one-child root
with open(OUT + "proteomes_kept.txt", "w") as f:
    f.write("\n".join(sorted(keep)) + "\n")
print(f"{len(keep)} proteomes with BUSCO >= {MIN_BUSCO}% in the tree")
