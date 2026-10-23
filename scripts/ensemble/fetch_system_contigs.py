"""Fetch NCBI feature tables for every contig that carries one of the anchor enzymes.

EvoWeaver needs each gene's position on its contig, and UniProt gives only the contig accession and the
locus tag. The feature table supplies the gene order those tags sit in. One gzip per system, appended as
it goes, so an interrupted run resumes instead of starting over.
"""

import gzip
import os
import sys

sys.path.insert(0, "_scripts/ensemble")
sys.path.insert(0, "_scripts/evoweaver/genomewide")
from fetch_contig_tables import BATCH, fetch_ft, parse_ft  # noqa: E402
from fetch_families import read_genes  # noqa: E402
from systems import HELDOUT, TRAIN  # noqa: E402
from concurrent.futures import ThreadPoolExecutor  # noqa: E402

OUT = "_data/ensemble/contigs/"
kept = set(open("_data/evoweaver/genomewide/proteomes_kept.txt").read().split())


def fetch_for(system):
    path = OUT + system + ".tsv.gz"
    wanted = sorted({c for r in read_genes("anchor_" + system) if r["proteome"] in kept
                     for c in r["contigs"]})
    done = {l.split("\t")[0] for l in gzip.open(path, "rt")} if os.path.exists(path) else set()
    todo = [c for c in wanted if c not in done]
    print(f"{system}: {len(wanted)} contigs carry the anchor, {len(todo)} still to fetch", flush=True)
    if not todo:
        return
    batches = [",".join(todo[i:i + BATCH]) for i in range(0, len(todo), BATCH)]
    with gzip.open(path, "at") as out, ThreadPoolExecutor(3) as pool:
        for n, text in enumerate(pool.map(fetch_ft, batches), 1):
            for contig, cds in parse_ft(text).items():
                cds.sort()
                for idx, (start, end, strand, tags, pid) in enumerate(cds, 1):
                    out.write(f"{contig}\t{idx}\t{start}\t{end}\t{strand}\t{','.join(tags)}\t{pid}\n")
            if n % 10 == 0 or n == len(batches):
                print(f"  {system}: {n}/{len(batches)} batches", flush=True)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for system in list(TRAIN) + list(HELDOUT):
        fetch_for(system)
    print("all systems done")
