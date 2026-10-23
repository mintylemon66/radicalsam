#!/usr/bin/env python3
"""Look for mqnP just beyond the 10-gene window at the futalosine DUF4229 loci that have none inside it.

    python3 _scripts/men_neighborhood/check_mqnp_beyond_window.py   # after analyze_futalosine.py

analyze_futalosine.py counts mqnP (NCBIfam NF041584) within WINDOW genes of
DUF4229. For every locus where no mqnP call falls inside, this script takes the
CDS from WINDOW + 1 to REACH genes on either side (duf_contig_cds.tsv.gz holds
every CDS of the contig), downloads those proteins from NCBI, and runs
hmmsearch with the MqnP model at its gathering threshold. It writes the nearest
mqnP at each locus to BEYOND_TSV and prints the counts by order. This step has
no control windows, so it only describes where mqnP sits; the enrichment test
is the 10-gene one in analyze_futalosine.py.
"""

import argparse
import gzip
import os
import shutil
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

from call_men_genes import read_tblout, run_hmmsearch
from fetch_window_proteins import bare_headers, fasta_ids, fetch_fasta

FUTALOSINE_DIR = "_data/men_neighborhood/futalosine/"
CDS_TSV = FUTALOSINE_DIR + "duf_contig_cds.tsv.gz"
WINDOWS_TSV = FUTALOSINE_DIR + "windows.tsv"
LOCI_TSV = FUTALOSINE_DIR + "duf4229_loci.tsv"
MQNP_HMM = FUTALOSINE_DIR + "hmm/NF041584.1.HMM"     # downloaded by call_men_genes.py --futalosine
BEYOND_FAA = FUTALOSINE_DIR + "beyond_window_proteins.faa"
BEYOND_TBL = FUTALOSINE_DIR + "hmmsearch_mqnp_beyond_window.tbl"
BEYOND_TSV = FUTALOSINE_DIR + "mqnp_beyond_window.tsv"

WINDOW = 10            # as in build_windows.py
REACH = 20
BATCH_SIZE = 500
THREADS = 3


def log(msg=""):
    print(msg, flush=True)


def read_tsv(path):
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def centers(path):
    """(proteome, duf_acc) -> (contig, idx, strand) of each DUF4229 gene."""
    found = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row["kind"] == "duf" and row["offset"] == "0":
                found[(row["proteome"], row["duf_acc"])] = (row["contig"], int(row["center_idx"]), row["center_strand"])
    return found


def ring_members(cds_path, loci, center):
    """locus -> CDS between WINDOW + 1 and REACH genes from DUF4229, with offset and strand relative to it."""
    wanted = defaultdict(list)
    for key in loci:
        wanted[center[key][0]].append(key)
    members = defaultdict(list)
    with gzip.open(cds_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row["contig"] not in wanted or not row["protein_id"]:
                continue
            for key in wanted[row["contig"]]:
                contig, idx, strand = center[key]
                step = int(row["idx"]) - idx
                if WINDOW < abs(step) <= REACH:
                    members[key].append({"offset": step if strand == "+" else -step,
                                         "rel_strand": "same" if row["strand"] == strand else "opposite",
                                         "protein_id": row["protein_id"], "product": row["product"]})
    return members


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--threads", type=int, default=THREADS)
    parser.add_argument("--cpus", type=int, default=os.cpu_count() or 4)
    args = parser.parse_args()
    if shutil.which("hmmsearch") is None:
        sys.exit("hmmsearch not on PATH; install HMMER 3")

    rows = {(r["proteome"], r["duf_acc"]): r for r in read_tsv(LOCI_TSV)}
    center = centers(WINDOWS_TSV)
    lacking = sorted(key for key, r in rows.items() if not r["MqnP_offset"])
    members = ring_members(CDS_TSV, lacking, center)
    wanted = sorted({m["protein_id"] for ms in members.values() for m in ms})
    have = fasta_ids(BEYOND_FAA)
    todo = [p for p in wanted if p not in have]
    log(f"{len(lacking):,} of {len(rows):,} DUF4229 loci have no mqnP within {WINDOW} genes; "
        f"{len(wanted):,} proteins {WINDOW + 1} to {REACH} genes away, {len(todo):,} still to fetch")

    batches = [todo[i:i + BATCH_SIZE] for i in range(0, len(todo), BATCH_SIZE)]
    with open(BEYOND_FAA, "a") as out, ThreadPoolExecutor(args.threads) as pool:
        for text in pool.map(fetch_fasta, batches):
            out.write(bare_headers(text))
    missing = [p for p in wanted if p not in fasta_ids(BEYOND_FAA)]
    log(f"  {len(wanted) - len(missing):,} sequences, {len(missing)} not returned")

    run_hmmsearch(MQNP_HMM, BEYOND_FAA, BEYOND_TBL, ["--cut_ga"], args.cpus)
    mqnp = {target for target, *_ in read_tblout(BEYOND_TBL)}

    nearest = {}
    for key in lacking:
        hits = [m for m in members.get(key, []) if m["protein_id"] in mqnp]
        if hits:
            nearest[key] = min(hits, key=lambda m: (abs(m["offset"]), m["offset"]))
    with open(BEYOND_TSV, "w") as out:
        out.write("proteome\tduf_acc\tgenus\torder\tMqnE_offset\tmqnP_offset\trel_strand\tprotein_id\tproduct\n")
        for key in lacking:
            r, m = rows[key], nearest.get(key)
            cells = [str(m["offset"]), m["rel_strand"], m["protein_id"], m["product"]] if m else ["", "", "", ""]
            out.write("\t".join([key[0], key[1], r["genus"], r["order"], r["MqnE_offset"]] + cells) + "\n")

    log()
    log(f"mqnP {WINDOW + 1} to {REACH} genes from DUF4229 at the loci with none within {WINDOW}, by order")
    by_order = defaultdict(list)
    for key in lacking:
        by_order[rows[key]["order"]].append(key)
    for order, keys in sorted(by_order.items(), key=lambda x: -len(x[1])):
        hits = [nearest[k] for k in keys if k in nearest]
        offsets = Counter(f"{m['offset']:+d}" for m in hits)
        same = sum(1 for m in hits if m["rel_strand"] == "same")
        log(f"  {order:<22} {len(hits):>4} of {len(keys):>4} loci; same strand {same}; "
            f"offsets {', '.join(f'{o}: {n}' for o, n in offsets.most_common(5))}")
    with_e = [k for k in lacking if rows[k]["MqnE_offset"]]
    log(f"  of these loci, {sum(1 for k in with_e if k in nearest)} of {len(with_e)} with mqnE inside the window "
        f"and {sum(1 for k in lacking if k in nearest and not rows[k]['MqnE_offset'])} of {len(lacking) - len(with_e)} "
        f"without it have mqnP {WINDOW + 1} to {REACH} genes away")
    within = sum(1 for r in rows.values() if r["MqnP_offset"]) + len(nearest)
    log(f"  mqnP within {REACH} genes of DUF4229 in all: {within:,} of {len(rows):,} loci")
    return 0


if __name__ == "__main__":
    sys.exit(main())
