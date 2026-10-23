"""Fetch every protein family the ensemble needs, one TSV per family.

For each of the nine training systems we need the anchor enzyme and its verified partner. We also need a
pool of unrelated families to act as negatives, drawn at random from Pfam with a fixed seed. The held-out
futalosine system is fetched too, since it has to be scored once the model is fit.

Each row carries the proteome, the locus tags and the EMBL contigs, which is what build_training_pairs.py
needs to place a gene on a contig. Families already on disk are skipped, so this can be rerun.
"""

import json
import os
import random
import re
import sys

from systems import HELDOUT, TRAIN

sys.path.insert(0, "_scripts/evoweaver/genomewide")
from fetch_families import GENES, fetch_genes, proteomes_with  # noqa: E402

OUT = "_data/ensemble/"
N_NEGATIVE_POOL = 60
SEED = 11
PFAM = "_hmm/pfamA/Pfam-A.hmm"
# anything in a pathway we train on, so a random draw cannot accidentally be a true partner
SKIP = re.compile(r"biotin|lipoyl|lipoate|pyrroloquinoline|mycofactocin|hydrogenase|molybdopterin|"
                  r"thiamin|formate|queuosine|menaquinone|futalosine|radical.SAM", re.I)


def wanted():
    """Every family to fetch, as {name: uniprot xref}.

    Anchors keep their system name. Every candidate family, whether a pathway gene or a random negative,
    is fetched once under fam_<accession>, so a family shared by two systems is downloaded only once.
    """
    out = {}
    for sys_name, s in list(TRAIN.items()) + list(HELDOUT.items()):
        out[f"anchor_{sys_name}"] = s["anchor"]
        out[f"fam_{s['partner'].split('-')[1]}"] = s["partner"]
        if "candidate_of_interest" in s:
            out[f"fam_{s['candidate_of_interest'].split('-')[1]}"] = s["candidate_of_interest"]
    if os.path.exists(OUT + "pathway_pfams.json"):          # written by map_pathway_genes.py
        for genes in json.load(open(OUT + "pathway_pfams.json")).values():
            for pf in genes.values():
                out[f"fam_{pf}"] = f"pfam-{pf}"
    return out


def negative_pool():
    """Random Pfam families, excluding anything named after a pathway we train on."""
    fams, acc = [], None
    for line in open(PFAM):
        if line.startswith("ACC"):
            acc = line.split()[1].split(".")[0]
        elif line.startswith("DESC") and acc and not SKIP.search(line):
            fams.append(acc)
    random.Random(SEED).shuffle(fams)
    return {f"fam_{a}": f"pfam-{a}" for a in fams[:N_NEGATIVE_POOL]}


if __name__ == "__main__":
    targets = wanted() | negative_pool()
    print(f"{len(targets)} families to fetch")
    for name, xref in targets.items():
        fetch_genes(name, f"xref:{xref}")
        print(f"  {name:28s} {xref:22s} {len(proteomes_with(name)):6d} proteomes", flush=True)
    print(f"\ndone, TSVs in {GENES}")
