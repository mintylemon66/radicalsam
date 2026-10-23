"""Every bacterial and archaeal reference proteome in UniProt, with its BUSCO score.

    python3 _scripts/evoweaver/genomewide/fetch_proteomes.py      # from the project folder

Step 1 of the genome-wide EvoWeaver run. EvoWeaver compares two families across a fixed set of
genomes, so the set has to be defined once, up front, and reused by every later step. Reference
proteomes are UniProt's one-per-species selection, which keeps over-sequenced genera from flooding
the comparison. The release is read back from the response header so the run can be dated.

Writes _data/evoweaver/genomewide/proteomes.tsv.
"""

import os
import urllib.parse
import urllib.request

OUT = "_data/evoweaver/genomewide/"
URL = "https://rest.uniprot.org/proteomes/stream"

# every bacterial and archaeal reference proteome, with its taxonomy and BUSCO completeness
params = {
    "query": "(reference:true) AND ((taxonomy_id:2) OR (taxonomy_id:2157))",
    "fields": "upid,organism_id,lineage,protein_count,cpd,busco",
    "format": "tsv",
}
os.makedirs(OUT, exist_ok=True)
r = urllib.request.urlopen(URL + "?" + urllib.parse.urlencode(params), timeout=900)
text = r.read().decode()
open(OUT + "proteomes.tsv", "w").write(text)
print(f"{text.count(chr(10)) - 1} reference proteomes, UniProt release {r.headers.get('X-UniProt-Release')}")
