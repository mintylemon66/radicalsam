#!/usr/bin/env python3
"""
Does the DUF4229 / MqnE / MqnP gene cluster sit inside a biosynthetic gene
cluster (BGC), i.e. one of the antibiotic-producing "factories" that make
Streptomyces famous?

This is Philipp's antiSMASH question, answered WITHOUT uploading anything:
antiSMASH maintains a public database of BGCs it has already computed for
~56,000 genomes. We look up our contigs there and test whether our cluster's
genomic coordinates fall inside any reported BGC.

Interpretation:
  overlap      -> DUF4229 may be tied to secondary metabolism; a more
                  complicated (and more careful) story.
  no overlap   -> the cluster is ordinary primary metabolism (a housekeeping
                  vitamin K2 operon), which is the cleaner reading and the one
                  the paper currently assumes.

Read-only HTTP queries. Nothing is submitted or published.
"""
import csv
import json
import os
import re
import time
import urllib.request
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
TSV = os.path.join(AUD, "snap", "neighbors_seq.tsv")
API = "https://antismash-db.secondarymetabolites.org/api/v2.0/search"
MQNP = re.compile(r"mqnp|menaquinone.*prenyltransferase", re.I)


def dom(p):
    s = set()
    if not os.path.exists(p):
        return s
    for line in open(p):
        if not line.startswith("#"):
            f = line.split()
            if len(f) >= 8:
                s.add(f[0])
    return s


def query(category, value, paginate=50):
    body = json.dumps({"query": {"search": "cluster", "return_type": "json",
                                 "terms": {"term_type": "expr",
                                           "category": category, "value": value}},
                       "paginate": paginate}).encode()
    req = urllib.request.Request(API, data=body,
                                 headers={"Content-Type": "application/json",
                                          "User-Agent": "radicalSAM-audit/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read().decode())
    except Exception:
        return {}


def main():
    duf = dom(os.path.join(AUD, "snap", "PF14012_vs_neighbors.domtbl"))
    ubia = dom(os.path.join(AUD, "cat_domtbl", "PF01040.domtbl"))

    rows = defaultdict(list)
    grp = {}
    with open(TSV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            grp[r["query_acc"]] = r["group"]
            rows[r["query_acc"]].append(r)

    def gi(r):
        try:
            return int(r["gene_index"])
        except (TypeError, ValueError):
            return None

    # locate our clusters and their genomic span
    clusters = []
    for acc, rr in rows.items():
        if grp.get(acc) != "mqn":
            continue
        inw = [r for r in rr if gi(r) is not None and abs(gi(r)) <= 10]
        d = [r for r in inw if r["n_protein_id"] in duf]
        p = [r for r in inw if r["n_protein_id"] in ubia and MQNP.search(r["n_product"] or "")]
        if not (d and p):
            continue
        coords = [int(r["n_start"]) for r in d + p] + [int(r["n_end"]) for r in d + p]
        coords += [int(rr[0]["q_start"]), int(rr[0]["q_end"])]
        clusters.append({"acc": acc, "contig": rr[0]["contig"],
                         "start": min(coords), "end": max(coords)})
    print(f"clusters with DUF4229 + MqnE + MqnP: {len(clusters)}")

    checked = overlap = notfound = 0
    results = []
    for c in clusters:
        base = c["contig"].split(".")[0]
        res = query("acc", base)
        found = res.get("clusters") or []
        # antiSMASH accessions are often NZ_-prefixed
        if not found:
            res = query("acc", "NZ_" + base)
            found = res.get("clusters") or []
        if not found:
            notfound += 1
            continue
        checked += 1
        hit = None
        for b in found:
            try:
                bs, be = int(b.get("start_pos")), int(b.get("end_pos"))
            except (TypeError, ValueError):
                continue
            if bs <= c["end"] and be >= c["start"]:
                hit = (b.get("term") or b.get("description") or "?", bs, be)
                break
        if hit:
            overlap += 1
            results.append((c["acc"], c["contig"], c["start"], c["end"], hit))
        time.sleep(0.3)

    print(f"\ncontigs found in antiSMASH-DB : {checked}")
    print(f"contigs not in antiSMASH-DB   : {notfound}")
    print(f"clusters OVERLAPPING a BGC    : {overlap}")
    if checked:
        print(f"overlap rate among checkable  : {100.0*overlap/checked:.1f}%")

    if results:
        print("\nOVERLAPS:")
        for acc, contig, s, e, (term, bs, be) in results[:15]:
            print(f"  {acc} {contig} cluster {s}-{e} inside BGC '{term}' {bs}-{be}")
    elif checked:
        print("\n  -> The DUF4229 / MqnE / MqnP cluster does NOT sit inside any")
        print("     annotated biosynthetic gene cluster. It behaves as ordinary")
        print("     primary metabolism (a housekeeping vitamin K2 operon), which")
        print("     is the cleaner reading and supports the paper's framing.")


if __name__ == "__main__":
    main()
