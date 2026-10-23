"""Download each family's proteins, with the locus tags and contigs that place them on a genome.

    python3 _scripts/evoweaver/genomewide/fetch_families.py [PF00001 ...]   # extra Pfam accessions optional

Step 3, and the module the later steps import. For every protein in a family it keeps the accession,
its locus tags, the proteome it belongs to and its EMBL contigs. The tags and contigs are what let
build_labels.py work out where a gene sits relative to MqnE, which EvoWeaver's gene-organization
methods need; presence/absence alone would only need the proteome column.

FAMILIES holds the five families the run always scores: MqnE and MqnC, the candidate DUF4229, and
two controls, MqnA/MqnD as a known pathway partner and Fer4 as an unrelated family. Results are
cached one TSV per family, so a rerun skips what is already downloaded.

Writes _data/evoweaver/genomewide/genes/<family>.tsv.
"""

import os
import re
import sys
import time
import urllib.parse
import urllib.request

OUT = "_data/evoweaver/genomewide/"
GENES = OUT + "genes/"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"
# every protein that belongs to a bacterial or archaeal reference proteome
SCOPE = "(keyword:KW-1185) AND ((taxonomy_id:2) OR (taxonomy_id:2157))"
FIELDS = "accession,gene_oln,gene_orf,xref_proteomes,xref_embl"

FAMILIES = {
    "MqnE": "xref:ncbifam-TIGR03700",
    "MqnC": "xref:ncbifam-TIGR03699",
    "DUF4229": "xref:pfam-PF14012",
    "MqnAD": "xref:pfam-PF02621",
    "Fer4": "xref:pfam-PF00037",
}


def fetch_genes(name, query, tries=5):
    # one row per protein: accession, locus tags, proteome, EMBL contigs; saved so reruns skip it
    path = GENES + name + ".tsv"
    if os.path.exists(path):
        return path
    params = {"query": f"({query}) AND {SCOPE}", "fields": FIELDS, "format": "tsv", "size": 500}
    url, pages = SEARCH + "?" + urllib.parse.urlencode(params), []
    while url:                                   # UniProt pages results; the Link header points to the next page
        for attempt in range(tries):
            try:
                r = urllib.request.urlopen(url, timeout=300)
                body = r.read().decode()
                break
            except Exception as e:
                print(f"  {name} retry {attempt + 1}: {e}", flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError(f"UniProt failed for {query}")
        pages.append(body if not pages else body.split("\n", 1)[1])
        nxt = re.search(r'<([^>]+)>; rel="next"', r.headers.get("Link") or "")
        url = nxt.group(1) if nxt else None
    text = "".join(pages)
    os.makedirs(GENES, exist_ok=True)
    with open(path + ".part", "w") as f:
        f.write(text)
    os.rename(path + ".part", path)
    return path


def read_genes(name):
    rows = []
    for line in open(GENES + name + ".tsv").read().splitlines()[1:]:
        if line.count("\t") != 4:      # UniProt's TSV has stray lines holding only a comma
            continue
        acc, oln, orf, prot, embl = line.split("\t")
        tags = (oln + " " + orf).split()
        ups = [p.strip().split(":")[0] for p in prot.split(";") if p.strip().startswith("UP")]
        contigs = [c for c in embl.split(";") if c]
        for up in ups:
            rows.append({"acc": acc, "proteome": up, "tags": tags, "contigs": contigs})
    return rows


def proteomes_with(name):
    return {r["proteome"] for r in read_genes(name)}


if __name__ == "__main__":
    families = dict(FAMILIES)
    for pf in sys.argv[1:]:
        families[pf] = f"xref:pfam-{pf}"
    for name, query in families.items():
        fetch_genes(name, query)
        print(f"{name}: {len(proteomes_with(name))} proteomes", flush=True)
