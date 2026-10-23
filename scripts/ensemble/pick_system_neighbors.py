#!/usr/bin/env python3
"""Neighborhood negatives for the retrained combiner: the families that recur near each anchor.

    python3 _scripts/ensemble/pick_system_neighbors.py      # from the project folder

The first combiner learned partners against random Pfam families. Near MqnE every candidate is
already close by, so it leaned on co-occurrence and ranked DUF4229 32nd of 146. Here each of the
eight training systems gets the negatives it meets in use: Pfam families that recur within 10 genes
of its anchor. The rule is the one that defined MqnE's 145 neighborhood families
(pick_neighbor_families.py): within 10 genes in at least 2% of anchor loci, at most 50,000 UniProt
proteins, and no pathway family. Every setting is fixed here before any result is seen, and nothing
from the futalosine system is used.

Up to 400 anchor loci per system (one per genome, seed 4229) are read from the contig feature
tables. The proteins within 10 genes are mapped to UniProt for their Pfam domains, and the families
that pass are downloaded genome-wide for EvoWeaver.

Writes _data/ensemble/system_neighbors.tsv and neighbor_pfam_map.tsv; gene lists go to genes/.
"""

import csv
import gzip
import json
import os
import random
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "_scripts/ensemble")
sys.path.insert(0, "_scripts/evoweaver/genomewide")
from fetch_families import fetch_genes, read_genes  # noqa: E402
from pick_null_families import uniprot_count  # noqa: E402
from systems import TRAIN  # noqa: E402

OUT = "_data/ensemble/"
KEPT = "_data/evoweaver/genomewide/proteomes_kept.txt"
PATHWAY = OUT + "pathway_pfams.json"
LOCI_PER_SYSTEM = 400
SEED = 4229
WINDOW_GENES = 10
MIN_SHARE = 0.02
MAX_ENTRIES = 50000
# the radical SAM domain every anchor carries, and DUF4229 itself, which must stay out of training
SKIP = {"PF04055", "PF14012"}
IDMAP = "https://rest.uniprot.org/idmapping/"
CHUNK = 90000


def log(msg):
    print(msg, flush=True)


def sample_loci(system, kept, rng):
    first = {}
    for row in read_genes("anchor_" + system):
        if row["proteome"] in kept and row["tags"] and row["contigs"]:
            first.setdefault(row["proteome"], row)
    genomes = sorted(first)
    rng.shuffle(genomes)
    return [first[g] for g in genomes[:LOCI_PER_SYSTEM]]


def neighbor_proteins(system, loci):
    """One list of protein IDs per anchor locus: the CDS within 10 genes on the anchor's contig."""
    wanted = {c for row in loci for c in row["contigs"]}
    genes = defaultdict(list)
    with gzip.open(f"{OUT}contigs/{system}.tsv.gz", "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 7 and parts[0] in wanted:
                genes[parts[0]].append((int(parts[1]), {t.casefold() for t in parts[5].split(",") if t}, parts[6]))
    windows = []
    for row in loci:
        tags = {t.casefold() for t in row["tags"]}
        for contig in row["contigs"]:
            at = next((idx for idx, t, _ in genes.get(contig, []) if t & tags), None)
            if at is not None:
                windows.append([pid for idx, _, pid in genes[contig] if 0 < abs(idx - at) <= WINDOW_GENES and pid])
                break
    return windows


def map_to_pfam(pids):
    """protein ID -> set of Pfam accessions, through UniProt ID mapping."""
    pfam = {}
    for start in range(0, len(pids), CHUNK):
        chunk = pids[start:start + CHUNK]
        data = urllib.parse.urlencode({"from": "EMBL-GenBank-DDBJ_CDS", "to": "UniProtKB", "ids": ",".join(chunk)}).encode()
        job = json.load(urllib.request.urlopen(IDMAP + "run", data, timeout=300))["jobId"]
        while True:
            time.sleep(10)
            status = json.load(urllib.request.urlopen(IDMAP + "status/" + job, timeout=300))
            if status.get("jobStatus") not in ("NEW", "RUNNING"):
                break
        if status.get("jobStatus") == "ERROR":
            raise RuntimeError(f"UniProt ID mapping failed: {status}")
        url = IDMAP + "uniprotkb/results/stream/" + job + "?" + urllib.parse.urlencode({"format": "tsv", "fields": "accession,xref_pfam"})
        text = urllib.request.urlopen(url, timeout=1200).read().decode()
        for line in text.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 3:
                pfam.setdefault(parts[0], set()).update(p for p in parts[2].split(";") if p)
        log(f"  mapped {start + len(chunk):,} of {len(pids):,} protein IDs")
    return pfam


def recurring(windows, pfam, skip):
    counts = defaultdict(int)
    for window in windows:
        for family in {f for pid in window for f in pfam.get(pid, ())} - skip:
            counts[family] += 1
    return {f: n for f, n in counts.items() if n >= MIN_SHARE * len(windows)}


def main():
    kept = set(open(KEPT).read().split())
    pathway = {pf for genes in json.load(open(PATHWAY)).values() for pf in genes.values()}
    skip = SKIP | pathway | {s["partner"].split("-")[1] for s in TRAIN.values()}
    rng = random.Random(SEED)

    windows = {}
    for system in TRAIN:
        loci = sample_loci(system, kept, rng)
        windows[system] = neighbor_proteins(system, loci)
        log(f"{system:14s} {len(loci)} anchor loci sampled, {len(windows[system])} placed on their contig")

    map_path = OUT + "neighbor_pfam_map.tsv"
    if os.path.exists(map_path):
        with open(map_path) as handle:
            pfam = {row[0]: set(filter(None, row[1].split(";"))) for row in csv.reader(handle, delimiter="\t") if row[2] == "1"}
    else:
        pids = sorted({pid for ws in windows.values() for w in ws for pid in w})
        log(f"{len(pids):,} neighbor proteins to map to UniProt")
        pfam = map_to_pfam(pids)
        with open(map_path, "w") as handle:
            for pid in pids:
                handle.write(f"{pid}\t{';'.join(sorted(pfam.get(pid, ())))}\t{int(pid in pfam)}\n")
    log(f"{len(pfam):,} neighbor proteins found in UniProtKB")

    picked = {s: recurring(windows[s], pfam, skip) for s in TRAIN}
    families = sorted({f for fams in picked.values() for f in fams})
    with ThreadPoolExecutor(6) as pool:
        entries = dict(zip(families, pool.map(uniprot_count, families)))
    failed = [f for f in families if entries[f] < 0]
    if failed:
        log(f"  UniProt count failed for {len(failed)} families, left out: {' '.join(failed)}")
    keep = {s: {f: n for f, n in fams.items() if 0 < entries[f] <= MAX_ENTRIES} for s, fams in picked.items()}

    with open(OUT + "system_neighbors.tsv", "w") as handle:
        handle.write("system\tacc\twindows_with\twindows\tuniprot_entries\n")
        for system, fams in keep.items():
            for f, n in sorted(fams.items()):
                handle.write(f"{system}\t{f}\t{n}\t{len(windows[system])}\t{entries[f]}\n")
    for system in TRAIN:
        log(f"{system:14s} {len(picked[system])} families in >= {MIN_SHARE:.0%} of windows, "
            f"{len(keep[system])} with at most {MAX_ENTRIES:,} UniProt proteins")

    todo = sorted({f for fams in keep.values() for f in fams})
    log(f"downloading gene lists for {len(todo)} families (ones already on disk are skipped)")
    with ThreadPoolExecutor(4) as pool:
        for done, _ in enumerate(pool.map(lambda a: fetch_genes(a, f"xref:pfam-{a}"), todo), 1):
            if done % 20 == 0 or done == len(todo):
                log(f"  gene lists: {done} of {len(todo)}")


if __name__ == "__main__":
    main()
