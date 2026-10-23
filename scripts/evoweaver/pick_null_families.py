"""The random Pfam families DUF4229 is scored against.

    python3 _scripts/evoweaver/genomewide/pick_null_families.py      # after fetch_families.py

Step 4. A high EvoWeaver score means nothing without knowing what score an arbitrary family gets, so
the run needs a null set. Families are drawn from the local Pfam-A library with a fixed seed, and
anything named after quinone biosynthesis is excluded so a true partner cannot land in the null by
accident.

A family much rarer or much commoner than DUF4229 would be easy to tell apart on prevalence alone,
so candidates are screened on UniProt protein count to roughly match it; build_labels.py then applies
the real 0.5-2x window on genome counts. The first 400 that pass are downloaded genome-wide.

Writes null_candidates.tsv and the per-family gene tables.
"""

import os
import random
import re
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from fetch_families import SCOPE, fetch_genes, proteomes_with

OUT = "_data/evoweaver/genomewide/"
PFAM = "_hmm/pfamA/Pfam-A.hmm"
SEED = 42
N_CANDIDATES = 400
DUF4229_PROTEOMES = 3380          # from fetch_families.py
# a cheap first screen on UniProt protein counts; a family with 0.5-2x DUF4229's genome count has at least
# 0.5x as many proteins, and fewer than 4x unless it averages over two copies per genome. build_labels.py then
# applies the real 0.5-2x window on genome counts.
LOW, HIGH = 0.5 * DUF4229_PROTEOMES, 4 * DUF4229_PROTEOMES
SKIP = re.compile(r"menaquinone|futalosine|quinone|PF14012|PF02621|PF00037", re.I)


def pfam_list():
    fams, acc, name = [], None, None
    for line in open(PFAM):
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("ACC"):
            acc = line.split()[1].split(".")[0]
        elif line.startswith("DESC"):
            fams.append((acc, name, line[6:].strip()))
    return fams


def uniprot_count(acc):
    q = f"(xref:pfam-{acc}) AND {SCOPE}"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": q, "size": 0})
    for _ in range(3):
        try:
            return int(urllib.request.urlopen(url, timeout=60).headers["X-Total-Results"])
        except Exception:
            pass
    return -1


if __name__ == "__main__":
    fams = [f for f in pfam_list() if not SKIP.search(" ".join(f))]
    random.Random(SEED).shuffle(fams)
    cand_file = OUT + "null_candidates.tsv"
    if not os.path.exists(cand_file):
        chosen = []
        with ThreadPoolExecutor(6) as pool:
            for start in range(0, len(fams), 60):
                block = fams[start:start + 60]
                for f, n in zip(block, pool.map(uniprot_count, [a for a, _, _ in block])):
                    if LOW <= n <= HIGH:
                        chosen.append((*f, n))
                print(f"screened {start + len(block)} families, {len(chosen)} in range", flush=True)
                if len(chosen) >= N_CANDIDATES:
                    break
        with open(cand_file, "w") as out:
            out.write("acc\tname\tdesc\tuniprot_entries\n")
            for row in chosen[:N_CANDIDATES]:
                out.write("\t".join(map(str, row)) + "\n")
    accs = [line.split("\t")[0] for line in open(cand_file).read().splitlines()[1:]]

    with ThreadPoolExecutor(4) as pool:
        for done, acc in enumerate(pool.map(lambda a: fetch_genes(a, f"xref:pfam-{a}") and a, accs), 1):
            print(f"{acc}: {len(proteomes_with(acc))} proteomes ({done}/{len(accs)})", flush=True)
