#!/usr/bin/env python3
"""Fetch the NCBI feature table of every contig that carries DUF4229 in a genome without MqnE.

    python3 _scripts/men_neighborhood/fetch_duf_contigs.py               # from the project root
    python3 _scripts/men_neighborhood/fetch_duf_contigs.py --futalosine  # the genomes WITH MqnE

DUF4229 (Pfam PF14012) is in 3,136 of the 16,025 kept reference proteomes, and
2,205 of those have no MqnE (NCBIfam TIGR03700). This script lists those
genomes and their DUF4229 proteins in TARGETS_TSV, then downloads the NCBI
feature table of every EMBL contig a DUF4229 protein sits on and keeps one row
per CDS in CDS_TSV. A rerun only fetches the contigs that are not in CDS_TSV yet.

With --futalosine the targets are the other 931 genomes, the ones that do carry
MqnE, and both files go to FUTALOSINE_DIR, so the files of the main run stay as
they are. build_windows.py, fetch_window_proteins.py and call_men_genes.py take
the same flag, and analyze_futalosine.py counts the result.

The fetch reuses the approach of _scripts/evoweaver/genomewide/fetch_contig_tables.py
(efetch rettype=ft, 20 contigs per request, 3 requests at a time, which is NCBI's
limit without an API key). This version also keeps each CDS's gene name,
product name and UniProtKB cross-reference, which the later scripts use.
"""

import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

GENOMEWIDE_DIR = "_data/evoweaver/genomewide/"
GENES_DIR = GENOMEWIDE_DIR + "genes/"
KEPT_TXT = GENOMEWIDE_DIR + "proteomes_kept.txt"
PROTEOMES_TSV = GENOMEWIDE_DIR + "proteomes.tsv"

OUT_DIR = "_data/men_neighborhood/"
TARGETS_TSV = OUT_DIR + "targets.tsv"
CDS_TSV = OUT_DIR + "duf_contig_cds.tsv.gz"
FUTALOSINE_DIR = OUT_DIR + "futalosine/"      # --futalosine: the same files for the DUF4229 genomes with MqnE

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{}.json"
BATCH_SIZE = 20
THREADS = 3
RETRIES = 5
MAX_CDS_SPAN = 100_000     # a CDS longer than this is a join across the origin of a circular contig

CDS_COLUMNS = ["contig", "idx", "start", "end", "strand", "locus_tags", "gene",
               "protein_id", "uniprot", "pseudo", "product"]


def log(msg):
    print(msg, flush=True)


def in_folder(path, futalosine):
    """path itself, or the file of the same name in FUTALOSINE_DIR for the --futalosine run."""
    return FUTALOSINE_DIR + path[len(OUT_DIR):] if futalosine else path


def read_genes(name):
    """One row per (protein, proteome) from a UniProt family table in GENES_DIR.

    Same reader as _scripts/evoweaver/genomewide/fetch_families.py.
    """
    rows = []
    for line in open(GENES_DIR + name + ".tsv").read().splitlines()[1:]:
        if line.count("\t") != 4:      # UniProt's TSV has stray lines holding only a comma
            continue
        acc, oln, orf, prot, embl = line.split("\t")
        tags = (oln + " " + orf).split()
        proteomes = [p.strip().split(":")[0] for p in prot.split(";") if p.strip().startswith("UP")]
        contigs = [c for c in embl.split(";") if c]
        for up in proteomes:
            rows.append({"acc": acc, "proteome": up, "tags": tags, "contigs": contigs})
    return rows


def split_lineage(lineage):
    """Phylum, order, family, genus and species from a UniProt lineage string."""
    names = [x.strip() for x in lineage.split(",")]
    phylum = names[3] if len(names) > 3 else ""
    order = next((x for x in names if x.endswith("ales")), "")
    family = next((x for x in names if x.endswith("aceae")), "")
    species = names[-1]
    genus = ""
    if family and names.index(family) + 1 < len(names):
        after = names[names.index(family) + 1]
        if " " not in after:
            genus = after
    if not genus:
        words = [w for w in species.split() if w != "Candidatus"]
        genus = words[0] if words else species
    return phylum, order, family, genus, species


def read_lineages(path):
    lineages = {}
    for line in open(path).read().splitlines()[1:]:
        parts = line.split("\t")
        lineages[parts[0]] = (parts[1], parts[2])
    return lineages


def embl_protein_ids(acc):
    """INSDC protein IDs of a UniProtKB entry, for placing a protein that has no locus tag."""
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(UNIPROT_ENTRY.format(acc), timeout=120) as response:
                entry = json.load(response)
            break
        except Exception as e:
            log(f"  {acc} retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))
    else:
        return []
    ids = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref["database"] == "EMBL":
            ids += [p["value"] for p in xref.get("properties", []) if p["key"] == "ProteinId" and p["value"] != "-"]
    return ids


def write_targets(path, futalosine=False):
    """Kept proteomes with DUF4229 and without MqnE (with MqnE if futalosine), one row per DUF4229 protein."""
    kept = set(open(KEPT_TXT).read().split())
    with_mqne = {r["proteome"] for r in read_genes("MqnE") if r["proteome"] in kept}
    with_duf = [r for r in read_genes("DUF4229") if r["proteome"] in kept]
    rows = [r for r in with_duf if (r["proteome"] in with_mqne) == futalosine]
    lineages = read_lineages(PROTEOMES_TSV)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as out:
        out.write("proteome\torganism_id\tphylum\torder\tfamily\tgenus\tspecies\tduf_acc\tlocus_tags\tcontigs\t"
                  "protein_ids\n")
        for r in sorted(rows, key=lambda r: (r["proteome"], r["acc"])):
            organism_id, lineage = lineages[r["proteome"]]
            phylum, order, family, genus, species = split_lineage(lineage)
            protein_ids = [] if r["tags"] else embl_protein_ids(r["acc"])     # the few entries with no locus tag
            out.write(f"{r['proteome']}\t{organism_id}\t{phylum}\t{order}\t{family}\t{genus}\t{species}\t"
                      f"{r['acc']}\t{','.join(r['tags'])}\t{','.join(r['contigs'])}\t{','.join(protein_ids)}\n")
    log(f"DUF4229 in {len({r['proteome'] for r in with_duf}):,} kept proteomes; "
        f"{len({r['proteome'] for r in rows}):,} have {'MqnE' if futalosine else 'no MqnE'} "
        f"({len(rows):,} DUF4229 proteins)")
    return rows


def to_int(field):
    return int(re.sub(r"\D", "", field))          # drops the < > ^ marks of partial features


def parse_feature_table(text):
    """NCBI feature table text -> {contig: [cds dict, ...]}, one dict per CDS, sorted by start."""
    features = []
    contig = None
    for line in text.splitlines():
        if line.startswith(">Feature"):
            accession = line.split("|")[1] if "|" in line else line.split()[1]
            contig = accession.split(".")[0]
            features.append({"contig": contig, "key": None})        # marks a contig with no features
            continue
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] and parts[2]:              # a new feature: start, end, key
            a, b = to_int(parts[0]), to_int(parts[1])
            features.append({"contig": contig, "key": parts[2], "spans": [(min(a, b), max(a, b))],
                             "strand": "-" if a > b else "+", "quals": []})
        elif len(parts) == 2 and parts[0] and features and features[-1]["key"]:   # more parts of a join
            a, b = to_int(parts[0]), to_int(parts[1])
            features[-1]["spans"].append((min(a, b), max(a, b)))
        elif len(parts) >= 4 and not parts[0] and features and features[-1]["key"]:
            features[-1]["quals"].append((parts[3], parts[4] if len(parts) > 4 else ""))

    contigs = {}
    gene = None
    for f in features:
        contigs.setdefault(f["contig"], [])
        if f["key"] is None:
            gene = None
            continue
        spans = f["spans"]
        lo, hi = min(a for a, b in spans), max(b for a, b in spans)
        if hi - lo > MAX_CDS_SPAN:
            lo, hi = spans[0]
        quals = f["quals"]
        tags = [v for k, v in quals if k in ("locus_tag", "old_locus_tag")]
        names = [v for k, v in quals if k == "gene"]
        if f["key"] == "gene":
            gene = {"lo": lo, "hi": hi, "tags": tags, "names": names}
            continue
        if f["key"] != "CDS":
            continue
        # the CDS takes the locus tag of the gene feature that contains it
        inside = gene is not None and gene["lo"] <= lo and hi <= gene["hi"]
        if not tags and inside:
            tags = gene["tags"]
        if not names and inside:
            names = gene["names"]
        protein_ids = [v.split("|")[1] if "|" in v else v for k, v in quals if k == "protein_id"]
        uniprot = [v.split(":", 1)[1] for k, v in quals if k == "db_xref" and v.startswith("UniProtKB/")]
        products = [v for k, v in quals if k == "product"]
        pseudo = any(k == "pseudo" for k, v in quals) or not protein_ids
        contigs[f["contig"]].append({
            "start": lo, "end": hi, "strand": f["strand"], "locus_tags": ",".join(tags),
            "gene": names[0] if names else "", "protein_id": protein_ids[0] if protein_ids else "",
            "uniprot": ",".join(uniprot), "pseudo": "1" if pseudo else "0",
            "product": (products[0] if products else "").replace("\t", " "),
        })
    for rows in contigs.values():
        rows.sort(key=lambda r: (r["start"], r["end"]))
    return contigs


def fetch_feature_tables(ids):
    params = urllib.parse.urlencode({"db": "nuccore", "id": ids, "rettype": "ft", "retmode": "text"})
    for attempt in range(RETRIES):
        try:
            time.sleep(1)
            return urllib.request.urlopen(EFETCH, params.encode(), timeout=600).read().decode()
        except Exception as e:
            log(f"  retry {attempt + 1}: {e}")
            time.sleep(5 * (attempt + 1))
    log(f"  gave up on {ids}")
    return ""


def fetched_contigs(path):
    if not os.path.exists(path):
        return set()
    return {line.split("\t", 1)[0] for line in gzip.open(path, "rt") if not line.startswith("contig\t")}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--futalosine", action="store_true",
                        help=f"target the DUF4229 genomes that carry MqnE; files go to {FUTALOSINE_DIR}")
    parser.add_argument("--targets", help=f"default {TARGETS_TSV}")
    parser.add_argument("--output", help=f"default {CDS_TSV}")
    parser.add_argument("--threads", type=int, default=THREADS)
    args = parser.parse_args()
    targets_tsv = args.targets or in_folder(TARGETS_TSV, args.futalosine)
    cds_tsv = args.output or in_folder(CDS_TSV, args.futalosine)

    rows = write_targets(targets_tsv, args.futalosine)
    wanted = sorted({c for r in rows for c in r["contigs"]})
    done = fetched_contigs(cds_tsv)
    todo = [c for c in wanted if c not in done]
    log(f"{len(wanted):,} contigs carry these DUF4229 proteins, {len(todo):,} still to fetch")

    new_file = not os.path.exists(cds_tsv)
    batches = [",".join(todo[i:i + BATCH_SIZE]) for i in range(0, len(todo), BATCH_SIZE)]
    with gzip.open(cds_tsv, "at") as out, ThreadPoolExecutor(args.threads) as pool:
        if new_file:
            out.write("\t".join(CDS_COLUMNS) + "\n")
        for n, text in enumerate(pool.map(fetch_feature_tables, batches), 1):
            for contig, cds in parse_feature_table(text).items():
                for idx, c in enumerate(cds, 1):
                    out.write("\t".join([contig, str(idx)] + [str(c[k]) for k in CDS_COLUMNS[2:]]) + "\n")
            if n % 10 == 0 or n == len(batches):
                log(f"  {n}/{len(batches)} batches")

    done = fetched_contigs(cds_tsv)
    missing = [c for c in wanted if c not in done]
    log(f"done; {len(missing)} contigs returned no CDS" + (f": {', '.join(missing[:20])}" if missing else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
