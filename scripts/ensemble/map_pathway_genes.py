"""Map each system's pathway genes to the Pfam family an HMM can detect.

Positives for the ensemble are the pathway genes of each system, not just the one cognate partner, since
one partner per system would leave too few positives to fit seven coefficients.

Gene sets come from KEGG modules where one exists, and otherwise from the system's named gene cluster.
Each gene symbol is resolved to a Pfam by asking UniProt for reviewed bacterial proteins with that symbol
and taking the family that most of them carry. The support column says how strong that agreement is, so a
symbol that resolves ambiguously can be dropped rather than guessed at.

The futalosine genes are mapped too, because the held-out system has to be scored, and never trained on.
"""

import collections
import json
import time
import urllib.parse
import urllib.request

OUT = "_data/ensemble/"
MIN_SUPPORT = 0.6      # at least this share of reviewed entries must agree on the family
MIN_ENTRIES = 3

# the anchor's own gene, which is the thing being scored against and so can never be a positive
ANCHOR_GENE = {"biotin": "bioB", "lipoyl": "lipA", "molybdopterin": "moaA", "thiamine": "thiH",
               "pqq": "pqqE", "mycofactocin": "mftC", "hydrogenase": "hydE", "queuosine": "queE",
               "pfl": "pflA", "futalosine": "mqnE"}

# PF04055 is the radical SAM domain every anchor carries, so it can never separate a partner from anything
# else. The precursor peptides are 25 to 40 residues, too short for a profile HMM to call reliably.
EXCLUDE_PFAM = {"PF04055"}
EXCLUDE_GENE = {"pqqA", "mftA"}

KEGG_MODULE = {
    "biotin": "M00123",
    "lipoyl": "M00881",
    "thiamine": "M00127",
    "futalosine": "M00930",      # held out, mapped only so it can be scored
}

# systems with no KEGG module: the named gene cluster, from the operon each pathway is known by
GENE_CLUSTER = {
    # M00880 lists descriptions rather than gene symbols, so the molybdopterin genes are named directly
    "molybdopterin": ["moaA", "moaB", "moaC", "moaD", "moaE", "mogA", "moeA", "mobA"],
    "pqq": ["pqqA", "pqqB", "pqqC", "pqqD", "pqqE", "pqqF"],
    "mycofactocin": ["mftA", "mftB", "mftC", "mftD", "mftE", "mftF"],
    "queuosine": ["queC", "queD", "queE", "queF"],
    "pfl": ["pflA", "pflB"],
}


def kegg_genes(module):
    """Gene symbols listed in a KEGG module's ORTHOLOGY block."""
    text = urllib.request.urlopen(f"https://rest.kegg.jp/get/{module}", timeout=60).read().decode()
    genes, inside = [], False
    for line in text.splitlines():
        if line.startswith("ORTHOLOGY"):
            inside = True
        elif inside and line[:1] not in (" ", "\t"):
            break
        if inside:
            parts = line.replace("ORTHOLOGY", "").split()
            if len(parts) >= 2 and parts[0].startswith("K"):
                for sym in parts[1].rstrip(";").split(","):
                    if sym and sym[0].islower():
                        genes.append(sym)
    return sorted(set(genes))


def pfam_for(symbol, tries=3):
    """The Pfam family most reviewed bacterial proteins of this gene symbol carry."""
    q = f"(gene:{symbol}) AND (reviewed:true) AND ((taxonomy_id:2) OR (taxonomy_id:2157))"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": q, "fields": "accession,xref_pfam", "format": "tsv", "size": 200})
    for _ in range(tries):
        try:
            rows = urllib.request.urlopen(url, timeout=120).read().decode().splitlines()[1:]
            break
        except Exception:
            time.sleep(5)
    else:
        return None
    counts = collections.Counter()
    for r in rows:
        parts = r.split("\t")
        if len(parts) > 1:
            for fam in {p for p in parts[1].split(";") if p.startswith("PF")}:
                counts[fam] += 1
    if not rows or not counts:
        return None
    fam, n = counts.most_common(1)[0]
    return {"pfam": fam, "entries": len(rows), "support": n / len(rows)}


if __name__ == "__main__":
    sets = {s: kegg_genes(m) for s, m in KEGG_MODULE.items()} | GENE_CLUSTER
    out = {}
    print(f"{'system':14s} {'gene':10s} {'pfam':10s} {'entries':>7s} {'support':>8s}  kept")
    for system, genes in sets.items():
        out[system] = {}
        for g in genes:
            hit = pfam_for(g)
            if hit is None:
                print(f"{system:14s} {g:10s} {'-':10s} {0:7d} {'-':>8s}  no reviewed entries")
                continue
            keep = (hit["entries"] >= MIN_ENTRIES and hit["support"] >= MIN_SUPPORT
                    and hit["pfam"] not in EXCLUDE_PFAM and g not in EXCLUDE_GENE
                    and g != ANCHOR_GENE.get(system))
            if keep:
                out[system][g] = hit["pfam"]
            print(f"{system:14s} {g:10s} {hit['pfam']:10s} {hit['entries']:7d} {hit['support']:8.0%}  {'yes' if keep else 'no'}")
    json.dump(out, open(OUT + "pathway_pfams.json", "w"), indent=1)
    print("\npositives per system:", {s: len(v) for s, v in out.items()})
    print("total:", sum(len(v) for s, v in out.items() if s != "futalosine"), "excluding the held-out system")
