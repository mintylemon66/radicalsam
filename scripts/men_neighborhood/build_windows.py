#!/usr/bin/env python3
"""List the genes within WINDOW genes of each DUF4229 gene, and of matched control genes.

    python3 _scripts/men_neighborhood/build_windows.py           # after fetch_duf_contigs.py
    python3 _scripts/men_neighborhood/build_windows.py --futalosine   # genomes WITH MqnE, in FUTALOSINE_DIR

For every DUF4229 protein in targets.tsv the script finds its CDS in
duf_contig_cds.tsv.gz (by UniProtKB accession, else by locus tag) and writes
the CDS within WINDOW genes on either side to WINDOWS_TSV. Offsets are counted
in the direction DUF4229 is transcribed, so +1 is the next gene downstream and
-1 the next gene upstream; rel_strand says whether a gene is on DUF4229's strand.

Control windows are the same windows around CONTROLS_PER_LOCUS random genes
drawn from the same contig. A control gene encodes a protein whose length is
within LENGTH_TOLERANCE of the DUF4229 protein's, and it sits more than
2 * WINDOW genes from any DUF4229 gene, so no control window overlaps a DUF4229
window. Because the controls come from the same contigs, each genome's own men
genes are equally available to both kinds of window.

With --futalosine every file is read from and written to FUTALOSINE_DIR (the
DUF4229 genomes that carry MqnE, from fetch_duf_contigs.py --futalosine); the
rules and the seed are the same.
"""

import argparse
import gzip
import random
import sys
from collections import defaultdict

OUT_DIR = "_data/men_neighborhood/"
TARGETS_TSV = OUT_DIR + "targets.tsv"
CDS_TSV = OUT_DIR + "duf_contig_cds.tsv.gz"
WINDOWS_TSV = OUT_DIR + "windows.tsv"
UNPLACED_TSV = OUT_DIR + "unplaced_duf4229.tsv"
FUTALOSINE_DIR = OUT_DIR + "futalosine/"      # --futalosine: the same files for the DUF4229 genomes with MqnE

WINDOW = 10
CONTROLS_PER_LOCUS = 5
LENGTH_TOLERANCE = 0.25        # control protein length within 25% of the DUF4229 protein's
SEED = 4229

WINDOW_COLUMNS = ["window_id", "kind", "proteome", "duf_acc", "contig", "center_idx", "center_strand",
                  "center_protein", "center_length", "member_idx", "offset", "rel_strand",
                  "protein_id", "locus_tags", "gene", "uniprot", "pseudo", "length_aa", "product"]


def log(msg):
    print(msg, flush=True)


def in_folder(path, futalosine):
    """path itself, or the file of the same name in FUTALOSINE_DIR for the --futalosine run."""
    return FUTALOSINE_DIR + path[len(OUT_DIR):] if futalosine else path


def read_targets(path):
    rows = []
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    for line in lines[1:]:
        rows.append(dict(zip(header, line.split("\t"))))
    return rows


def read_cds(path):
    """contig -> list of CDS dicts ordered by idx (idx 1 is element 0)."""
    contigs = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            row["idx"] = int(row["idx"])
            row["start"], row["end"] = int(row["start"]), int(row["end"])
            row["length_aa"] = (row["end"] - row["start"] + 1) // 3 - 1
            contigs[row["contig"]].append(row)
    for rows in contigs.values():
        rows.sort(key=lambda r: r["idx"])
    return contigs


def place(target, contigs, by_uniprot, by_tag, by_gene, by_protein):
    """(contig, idx, how) for a DUF4229 protein, or None when its CDS is not found.

    Tries the UniProtKB cross-reference, then the locus tag, then the gene
    qualifier, since some GenBank records (BA000036, C. glutamicum) keep the
    ordered locus name UniProt uses (Cgl0447) in the gene qualifier, and last
    the INSDC protein ID that fetch_duf_contigs.py looked up for entries with
    no locus tag.
    """
    for contig in target["contigs"].split(","):
        for c, idx in by_uniprot.get(target["duf_acc"], ()):
            if c == contig:
                return contig, idx, "uniprot"
    for index, how in ((by_tag, "locus_tag"), (by_gene, "gene_qualifier")):
        for contig in target["contigs"].split(","):
            for tag in target["locus_tags"].split(","):
                idx = index.get((contig, tag.casefold()))
                if tag and idx:
                    return contig, idx, how
    for pid in target.get("protein_ids", "").split(","):
        if pid in by_protein and by_protein[pid][0] in target["contigs"].split(","):
            return by_protein[pid][0], by_protein[pid][1], "protein_id"
    return None


def window_rows(window_id, kind, target, contig, center_idx, cds):
    """One row per CDS within WINDOW genes of center_idx, the center itself at offset 0."""
    center = cds[center_idx - 1]
    rows = []
    for member in cds[max(0, center_idx - 1 - WINDOW):center_idx + WINDOW]:
        step = member["idx"] - center_idx
        offset = step if center["strand"] == "+" else -step
        rows.append([window_id, kind, target["proteome"], target["duf_acc"], contig, center_idx,
                     center["strand"], center["protein_id"], center["length_aa"], member["idx"], offset,
                     "same" if member["strand"] == center["strand"] else "opposite",
                     member["protein_id"], member["locus_tags"], member["gene"], member["uniprot"],
                     member["pseudo"], member["length_aa"], member["product"]])
    return rows


def pick_controls(rng, cds, duf_positions, center_length):
    """Up to CONTROLS_PER_LOCUS random length-matched genes away from every DUF4229 window."""
    low, high = center_length * (1 - LENGTH_TOLERANCE), center_length * (1 + LENGTH_TOLERANCE)
    candidates = [c["idx"] for c in cds
                  if c["pseudo"] == "0" and low <= c["length_aa"] <= high
                  and all(abs(c["idx"] - d) > 2 * WINDOW for d in duf_positions)]
    return sorted(rng.sample(candidates, min(CONTROLS_PER_LOCUS, len(candidates))))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--futalosine", action="store_true",
                        help=f"the DUF4229 genomes that carry MqnE; files are read from and written to {FUTALOSINE_DIR}")
    parser.add_argument("--targets", help=f"default {TARGETS_TSV}")
    parser.add_argument("--cds", help=f"default {CDS_TSV}")
    parser.add_argument("--output", help=f"default {WINDOWS_TSV}")
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    windows_tsv = args.output or in_folder(WINDOWS_TSV, args.futalosine)
    unplaced_tsv = in_folder(UNPLACED_TSV, args.futalosine)

    targets = read_targets(args.targets or in_folder(TARGETS_TSV, args.futalosine))
    contigs = read_cds(args.cds or in_folder(CDS_TSV, args.futalosine))
    log(f"{len(targets):,} DUF4229 proteins; {len(contigs):,} contigs with "
        f"{sum(len(v) for v in contigs.values()):,} CDS")

    by_uniprot, by_tag, by_gene, by_protein = defaultdict(list), {}, {}, {}
    for contig, cds in contigs.items():
        for c in cds:
            for acc in c["uniprot"].split(","):
                if acc:
                    by_uniprot[acc].append((contig, c["idx"]))
            for tag in c["locus_tags"].split(","):
                if tag:
                    by_tag[(contig, tag.casefold())] = c["idx"]     # tag case differs between UniProt and GenBank
            if c["gene"]:
                by_gene[(contig, c["gene"].casefold())] = c["idx"]
            if c["protein_id"]:
                by_protein[c["protein_id"]] = (contig, c["idx"])

    placed, unplaced = [], []
    for t in targets:
        spot = place(t, contigs, by_uniprot, by_tag, by_gene, by_protein)
        (placed if spot else unplaced).append((t, spot))
    how = defaultdict(int)
    for t, spot in placed:
        how[spot[2]] += 1
    log(f"placed {len(placed):,} DUF4229 genes ({dict(how)}); {len(unplaced)} not found on their contigs")

    duf_positions = defaultdict(list)
    for t, (contig, idx, _) in placed:
        duf_positions[contig].append(idx)

    rng = random.Random(args.seed)
    n_windows = defaultdict(int)
    short_of_controls = 0
    with open(windows_tsv, "w") as out:
        out.write("\t".join(WINDOW_COLUMNS) + "\n")
        for t, (contig, idx, _) in sorted(placed, key=lambda x: (x[0]["proteome"], x[0]["duf_acc"])):
            cds = contigs[contig]
            locus = f"{t['proteome']}|{t['duf_acc']}"
            for row in window_rows(locus + "|duf", "duf", t, contig, idx, cds):
                out.write("\t".join(map(str, row)) + "\n")
            n_windows["duf"] += 1
            controls = pick_controls(rng, cds, duf_positions[contig], cds[idx - 1]["length_aa"])
            short_of_controls += len(controls) < CONTROLS_PER_LOCUS
            for k, center in enumerate(controls, 1):
                for row in window_rows(f"{locus}|control{k}", "control", t, contig, center, cds):
                    out.write("\t".join(map(str, row)) + "\n")
                n_windows["control"] += 1

    with open(unplaced_tsv, "w") as out:
        out.write("proteome\tduf_acc\tlocus_tags\tcontigs\tcontig_fetched\n")
        for t, _ in unplaced:
            fetched = any(c in contigs for c in t["contigs"].split(","))
            out.write(f"{t['proteome']}\t{t['duf_acc']}\t{t['locus_tags']}\t{t['contigs']}\t{int(fetched)}\n")

    log(f"wrote {windows_tsv}: {n_windows['duf']:,} DUF4229 windows, {n_windows['control']:,} control windows; "
        f"{short_of_controls:,} loci had fewer than {CONTROLS_PER_LOCUS} eligible control genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
