#!/usr/bin/env python3
"""Download the protein sequence of every gene in the DUF4229 and control windows from NCBI.

    python3 _scripts/men_neighborhood/fetch_window_proteins.py   # after build_windows.py
    python3 _scripts/men_neighborhood/fetch_window_proteins.py --futalosine   # genomes WITH MqnE

Reads the INSDC protein_id of every CDS in windows.tsv and fetches the
sequences with efetch (db=protein, rettype=fasta), BATCH_SIZE accessions per
request and 3 requests at a time. Each FASTA header is cut to the bare
protein_id so it matches windows.tsv. A rerun only fetches accessions that are
not in PROTEINS_FAA yet; accessions NCBI does not return are listed at the end.
With --futalosine both files are the ones in FUTALOSINE_DIR.
"""

import argparse
import os
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

OUT_DIR = "_data/men_neighborhood/"
WINDOWS_TSV = OUT_DIR + "windows.tsv"
PROTEINS_FAA = OUT_DIR + "window_proteins.faa"
FUTALOSINE_DIR = OUT_DIR + "futalosine/"      # --futalosine: the same files for the DUF4229 genomes with MqnE

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 500
THREADS = 3
RETRIES = 5


def log(msg):
    print(msg, flush=True)


def in_folder(path, futalosine):
    """path itself, or the file of the same name in FUTALOSINE_DIR for the --futalosine run."""
    return FUTALOSINE_DIR + path[len(OUT_DIR):] if futalosine else path


def wanted_ids(path):
    ids = set()
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col = header.index("protein_id")
        for line in handle:
            pid = line.rstrip("\n").split("\t")[col]
            if pid:
                ids.add(pid)
    return sorted(ids)


def fasta_ids(path):
    if not os.path.exists(path):
        return set()
    return {line[1:].split()[0] for line in open(path) if line.startswith(">")}


def fetch_fasta(ids):
    params = urllib.parse.urlencode({"db": "protein", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"})
    for attempt in range(RETRIES):
        try:
            time.sleep(1)
            return urllib.request.urlopen(EFETCH, params.encode(), timeout=600).read().decode()
        except Exception as e:
            log(f"  retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))
    log(f"  gave up on a batch starting {ids[0]}")
    return ""


def bare_headers(text):
    """'>CAF19163.1 hypothetical protein [...]' -> '>CAF19163.1'; drops blank lines."""
    lines = []
    for line in text.splitlines():
        if line.startswith(">"):
            lines.append(line.split()[0])
        elif line.strip():
            lines.append(line.strip())
    return "\n".join(lines) + "\n" if lines else ""


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--futalosine", action="store_true",
                        help=f"the DUF4229 genomes that carry MqnE; files are read from and written to {FUTALOSINE_DIR}")
    parser.add_argument("--windows", help=f"default {WINDOWS_TSV}")
    parser.add_argument("--output", help=f"default {PROTEINS_FAA}")
    parser.add_argument("--threads", type=int, default=THREADS)
    args = parser.parse_args()
    proteins_faa = args.output or in_folder(PROTEINS_FAA, args.futalosine)

    wanted = wanted_ids(args.windows or in_folder(WINDOWS_TSV, args.futalosine))
    have = fasta_ids(proteins_faa)
    todo = [p for p in wanted if p not in have]
    log(f"{len(wanted):,} proteins in the windows, {len(todo):,} still to fetch")

    batches = [todo[i:i + BATCH_SIZE] for i in range(0, len(todo), BATCH_SIZE)]
    with open(proteins_faa, "a") as out, ThreadPoolExecutor(args.threads) as pool:
        for n, text in enumerate(pool.map(fetch_fasta, batches), 1):
            out.write(bare_headers(text))
            if n % 50 == 0 or n == len(batches):
                log(f"  {n}/{len(batches)} batches")

    have = fasta_ids(proteins_faa)
    missing = [p for p in wanted if p not in have]
    log(f"done; {len(have):,} sequences, {len(missing)} not returned" +
        (f": {', '.join(missing[:20])}" if missing else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
