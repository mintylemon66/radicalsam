"""The harder comparison set: Pfam families that recur within 10 genes of MqnE.

    python3 _scripts/evoweaver/genomewide/pick_neighbor_families.py      # after pick_null_families.py

Step 5. Random families are an easy comparison, because DUF4229 was found by its proximity to MqnE
and they have none. These families were picked by that same proximity, so they share the property
that selected the candidate and the comparison tests whether anything beyond proximity separates it.

The rule is fixed here before any score is seen: inside 10 genes of MqnE in at least 2% of its
neighborhoods, from the sequence-based Pfam scan rather than product names, and at most 50,000
UniProt proteins, which leaves out only the giant families that would take hours to download. MqnE's
own domain, DUF4229 and the two controls are excluded.

Writes neighbor_candidates.tsv and the per-family gene tables.
"""

import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

from fetch_families import fetch_genes, proteomes_with
from pick_null_families import uniprot_count

OUT = "_data/evoweaver/genomewide/"
WINDOW_GENES = 10
MIN_SHARE = 0.02      # in at least 2% of MqnE neighborhoods
MAX_ENTRIES = 50000   # leaves out only the giant families (ABC transporters, regulators) that would take hours to download
# DUF4229 itself, the two controls already in the run, and Radical_SAM (MqnE's own domain)
SKIP = {"PF14012", "PF02621", "PF00037", "PF04055"}

# Pfam families that recur within 10 genes of MqnE: the comparison set for gene organization,
# since every one of them was picked by the same proximity to MqnE that picked DUF4229
kind = {r["accession"]: r["seq_assignment"] for r in csv.DictReader(open("_audit/mqn_seq_classification.tsv"), delimiter="\t")}
seen, windows, names = defaultdict(set), set(), {}
for line in open("_audit/pfam_wide.domtbl"):
    if line.startswith("#"):
        continue
    p = line.split()
    if len(p) < 23:
        continue
    q, group, gi = p[-1].split("|")[:3]
    if kind.get(q) != "MqnE":
        continue
    windows.add(q)
    gi = int(gi.split("=")[1])
    acc = p[4].split(".")[0]
    if gi != 0 and abs(gi) <= WINDOW_GENES and acc not in SKIP:
        seen[acc].add(q)
        names[acc] = p[3]

common = sorted(a for a in seen if len(seen[a]) >= MIN_SHARE * len(windows))
with ThreadPoolExecutor(6) as pool:
    entries = dict(zip(common, pool.map(uniprot_count, common)))
keep = [a for a in common if entries[a] <= MAX_ENTRIES]
with open(OUT + "neighbor_candidates.tsv", "w") as out:
    out.write("acc\tname\tmqne_windows\tuniprot_entries\n")
    for a in keep:
        out.write(f"{a}\t{names[a]}\t{len(seen[a])}\t{entries[a]}\n")
print(f"{len(common)} families in >= {MIN_SHARE:.0%} of {len(windows)} MqnE neighborhoods; {len(keep)} with at most {MAX_ENTRIES:,} UniProt proteins")

with ThreadPoolExecutor(4) as pool:
    for done, acc in enumerate(pool.map(lambda a: fetch_genes(a, f"xref:pfam-{a}") and a, keep), 1):
        print(f"{acc}: {len(proteomes_with(acc))} proteomes ({done}/{len(keep)})", flush=True)
