"""Where every gene sits on the contigs that carry MqnE.

    python3 _scripts/evoweaver/genomewide/fetch_contig_tables.py      # after fetch_families.py

Step 6. UniProt gives a protein's locus tags and contigs but not its position, and EvoWeaver's
gene-organization methods need an index and a strand. NCBI feature tables supply both, so every
contig carrying MqnE is fetched and its CDS features are parsed into (position, strand, locus tag).

Two cases the parser has to handle: a CDS joined across the origin of a circular contig, which would
otherwise span the whole genome, and a CDS whose locus tag sits on an enclosing gene line instead of
its own. Output is appended and already-fetched contigs are skipped, so an interrupted run resumes.

Writes mqne_contig_cds.tsv.gz.
"""

import gzip
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
import urllib.parse
import urllib.request

from fetch_families import read_genes

OUT = "_data/evoweaver/genomewide/"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 20

kept = set(open(OUT + "proteomes_kept.txt").read().split())


def parse_ft(text):
    # NCBI feature table -> {contig: [(start, end, strand, locus_tags, protein_id), ...]} for every CDS
    contigs, contig = {}, None
    feat, spans, tags, pid = None, [], [], ""
    gene_span, gene_tags = None, []

    def close():
        if feat == "CDS" and spans:
            lo, hi = min(a for a, b in spans), max(b for a, b in spans)
            if hi - lo > 100000:           # a join across the origin of a circular contig: order by its first part
                lo, hi = spans[0]
            # the CDS takes the locus tag of the gene line that contains it, even with other features in between
            inherited = gene_tags if gene_span and gene_span[0] <= lo and hi <= gene_span[1] else []
            contigs[contig].append((lo, hi, strand, tags or inherited, pid))

    for line in text.splitlines():
        if line.startswith(">Feature"):
            close()
            acc = line.split("|")[1] if "|" in line else line.split()[1]
            contig, feat, spans, gene_span, gene_tags = acc.split(".")[0], None, [], None, []
            contigs.setdefault(contig, [])
            continue
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] and parts[2]:          # new feature: start, end, key
            close()
            a, b = int(re.sub(r"\D", "", parts[0])), int(re.sub(r"\D", "", parts[1]))   # drops < > ^ marks
            feat, spans, tags, pid = parts[2], [(min(a, b), max(a, b))], [], ""
            strand = "1" if a > b else "0"
        elif len(parts) == 2 and parts[0] and feat:              # further parts of a joined feature
            a, b = int(re.sub(r"\D", "", parts[0])), int(re.sub(r"\D", "", parts[1]))
            spans.append((min(a, b), max(a, b)))
        elif len(parts) >= 5 and parts[3] in ("locus_tag", "old_locus_tag", "gene"):
            tags = tags + [parts[4]]
            if feat == "gene":
                gene_span = (min(x for x, _ in spans), max(y for _, y in spans))
                gene_tags = tags
        elif len(parts) >= 5 and parts[3] == "protein_id":
            pid = parts[4].split("|")[1] if "|" in parts[4] else parts[4]
    close()
    return contigs


def fetch_ft(ids):
    params = urllib.parse.urlencode({"db": "nuccore", "id": ids, "rettype": "ft", "retmode": "text"})
    for attempt in range(5):
        try:
            time.sleep(1)
            return urllib.request.urlopen(EFETCH, params.encode(), timeout=600).read().decode()
        except Exception as e:
            print(f"  retry {attempt + 1}: {e}", flush=True)
            time.sleep(5 * (attempt + 1))
    print(f"  gave up on {ids}", flush=True)
    return ""


if __name__ == "__main__":
    wanted = sorted({c for r in read_genes("MqnE") if r["proteome"] in kept for c in r["contigs"]})
    path = OUT + "mqne_contig_cds.tsv.gz"
    done = set()
    if os.path.exists(path):
        done = {line.split("\t")[0] for line in gzip.open(path, "rt")}
    todo = [c for c in wanted if c not in done]
    print(f"{len(wanted)} contigs carry MqnE, {len(todo)} still to fetch", flush=True)
    batches = [",".join(todo[i:i + BATCH]) for i in range(0, len(todo), BATCH)]
    with gzip.open(path, "at") as out, ThreadPoolExecutor(3) as pool:      # NCBI allows 3 requests/s without a key
        for n, text in enumerate(pool.map(fetch_ft, batches), 1):
            for contig, cds in parse_ft(text).items():
                cds.sort()
                for idx, (start, end, strand, tags, pid) in enumerate(cds, 1):
                    out.write(f"{contig}\t{idx}\t{start}\t{end}\t{strand}\t{','.join(tags)}\t{pid}\n")
            print(f"  {n}/{len(batches)} batches", flush=True)
