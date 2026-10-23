#!/usr/bin/env python3
"""Record which kept reference proteomes encode each classical menaquinone (men) enzyme.

    python3 _scripts/men_neighborhood/fetch_genome_men.py        # from the project root

The windows only show what sits next to DUF4229. Whether a genome runs the
classical (Men) pathway at all is a whole-genome question, answered here from
UniProtKB: for every NCBIfam model in GENOME_MODELS the script downloads the
UniProtKB entries of bacterial and archaeal reference proteomes that InterPro
matches to that model (InterProScan runs each NCBIfam HMM at the model's own
cutoff), and writes one row per kept proteome with the number of proteins per
gene to GENOME_TSV. A genome counts as running the classical pathway when it
has at least one MenD and one MenB.

The raw UniProt tables are saved in UNIPROT_DIR, so a rerun skips them. The
query is the one in _scripts/evoweaver/genomewide/fetch_families.py.
"""

import argparse
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

GENOMEWIDE_DIR = "_data/evoweaver/genomewide/"
KEPT_TXT = GENOMEWIDE_DIR + "proteomes_kept.txt"

OUT_DIR = "_data/men_neighborhood/"
UNIPROT_DIR = OUT_DIR + "uniprot/"
GENOME_TSV = OUT_DIR + "genome_men.tsv"

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
SCOPE = "(keyword:KW-1185) AND ((taxonomy_id:2) OR (taxonomy_id:2157))"   # reference proteomes
FIELDS = "accession,xref_proteomes,length"
PAGE_SIZE = 500
RETRIES = 5

# gene -> NCBIfam models; the models with Actinomycetota in their range are listed because
# TIGR01923 (MenE) matches only about 135 actinobacterial proteins in UniProtKB
GENOME_MODELS = {
    "MenD": ["TIGR00173"],
    "MenB": ["TIGR01929", "NF006186", "NF005637"],
    "MenA": ["TIGR00751"],
    "MenE": ["TIGR01923", "NF005877", "NF002966"],
    "MqnP": ["NF041584"],
}


def log(msg):
    print(msg, flush=True)


def fetch_model(model):
    """UniProtKB entries matched to one NCBIfam model, saved as TSV; returns the path."""
    path = UNIPROT_DIR + model + ".tsv"
    if os.path.exists(path):
        return path
    params = {"query": f"(xref:ncbifam-{model}) AND {SCOPE}", "fields": FIELDS,
              "format": "tsv", "size": PAGE_SIZE}
    url = SEARCH + "?" + urllib.parse.urlencode(params)
    pages = []
    while url:                                   # UniProt pages results; the Link header points to the next page
        for attempt in range(RETRIES):
            try:
                response = urllib.request.urlopen(url, timeout=300)
                body = response.read().decode()
                break
            except Exception as e:
                log(f"  {model} retry {attempt + 1}: {e}")
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError(f"UniProt failed for {model}")
        pages.append(body if not pages else body.split("\n", 1)[1])
        found = re.search(r'<([^>]+)>; rel="next"', response.headers.get("Link") or "")
        url = found.group(1) if found else None
    os.makedirs(UNIPROT_DIR, exist_ok=True)
    with open(path + ".part", "w") as out:
        out.write("".join(pages))
    os.rename(path + ".part", path)
    return path


def read_model(path):
    """proteome -> set of UniProt accessions."""
    by_proteome = defaultdict(set)
    for line in open(path).read().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        for p in parts[1].split(";"):
            p = p.strip()
            if p.startswith("UP"):
                by_proteome[p.split(":")[0]].add(parts[0])
    return by_proteome


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", default=GENOME_TSV)
    args = parser.parse_args()

    kept = open(KEPT_TXT).read().split()
    proteins = {gene: defaultdict(set) for gene in GENOME_MODELS}
    for gene, models in GENOME_MODELS.items():
        for model in models:
            per_proteome = read_model(fetch_model(model))
            log(f"{gene} {model}: {sum(len(v) for v in per_proteome.values()):,} proteins "
                f"in {len(per_proteome):,} reference proteomes")
            for up, accs in per_proteome.items():
                proteins[gene][up] |= accs

    genes = list(GENOME_MODELS)
    with open(args.output, "w") as out:
        out.write("proteome\t" + "\t".join(f"n_{g}" for g in genes) + "\tclassical\t"
                  + "\t".join(f"acc_{g}" for g in genes) + "\n")
        for up in kept:
            counts = [len(proteins[g].get(up, ())) for g in genes]
            classical = int(bool(proteins["MenD"].get(up)) and bool(proteins["MenB"].get(up)))
            accs = [",".join(sorted(proteins[g].get(up, ()))) for g in genes]
            out.write(f"{up}\t" + "\t".join(map(str, counts)) + f"\t{classical}\t" + "\t".join(accs) + "\n")
    n_classical = sum(1 for up in kept if proteins["MenD"].get(up) and proteins["MenB"].get(up))
    log(f"wrote {args.output}: {len(kept):,} kept proteomes, {n_classical:,} with MenD and MenB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
