#!/usr/bin/env python3
"""Build the sequence package Prof. Liu asked for in the 2026-07-26 lab meeting.

His request, verbatim from the meeting: "send the sequence and the strain ID, so the
accession number of the bacterium and the sequence of the MqnE and MqnC and the DUF one.
And then we can put the vectors together."

So for each candidate strain this assembles, from one and the same organism:
  - organism name and strain designation
  - the nucleotide accession the genes sit on (the "accession number of the bacterium")
  - MqnE     protein accession + sequence   (gets the pull-down tag)
  - MqnC     protein accession + sequence
  - DUF4229  protein accession + sequence   (the small one, about 100 aa)
  - the operon arrangement: gene order, strands, intergenic gaps, any CDS overlap

Strains are ranked so the lab can take the top one and go. Ranking prefers a
lab-friendly genus, a tight same-strand MqnE-DUF4229 operon, presence of the AsnC/Lrp
regulator, and a complete set of three full-length sequences.

Run:  python3 scripts/build_construct_package.py
Writes into wetlab_constructs/ only. Reads everything else read-only.
"""
import collections
import csv
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "wetlab_constructs")

NEIGH = os.path.join(ROOT, "audit", "neighbors_seq.tsv")
NEIGH_FAA = os.path.join(ROOT, "audit", "neighbors_seq.faa")
DOMTBL = os.path.join(ROOT, "audit", "pfam_wide.domtbl")
TARGETS = os.path.join(ROOT, "data", "mqn_targets.fasta")
CLASSIF = os.path.join(ROOT, "audit", "mqn_seq_classification.tsv")
META = os.path.join(ROOT, "data", "metadata_org.tsv")
CATS = os.path.join(ROOT, "ssn", "genome_categories.tsv")

WINDOW_BP = 8000
WINDOW_GENES = 10

# Genera the lab already works with, in preference order.
PREFERRED_GENERA = ["Streptomyces", "Micromonospora", "Kitasatospora"]

# Genera carrying a classical Men pathway signature (mencheck, NOTEBOOK 43). A menaquinone
# readout in these is ambiguous because they can make MK by the classical route too.
AVOID_GENERA = {"Microbacterium", "Mycolicibacterium", "Mycobacteroides",
                "Actinoalloteichus", "Streptomonospora", "Acidithrix", "Nocardioides"}


def read_fasta(path):
    seqs, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            head = line[1:].strip()
            if head.startswith(("sp|", "tr|")):
                acc = head.split("|")[1]
            else:
                acc = head.split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def in_window(gap, gi):
    return gap <= WINDOW_BP and abs(gi) <= WINDOW_GENES


def main():
    os.makedirs(OUT, exist_ok=True)

    classif = {r["accession"]: r["seq_assignment"]
               for r in csv.DictReader(open(CLASSIF), delimiter="\t")}

    organism, phylum = {}, {}
    for r in csv.DictReader(open(META), delimiter="\t"):
        m = re.match(r"^(?:sp|tr)\|([^|]+)\|", r["Sequence_ID"])
        a = m.group(1) if m else r["Sequence_ID"]
        organism[a] = r.get("Organism_Name", "")
        phylum[a] = r.get("Phylum", "")

    # regulator flag per anchor, from the genome-category sort
    regulator = {}
    if os.path.exists(CATS):
        for r in csv.DictReader(open(CATS), delimiter="\t"):
            regulator[r["accession"]] = r.get("regulator", "?")

    # DUF4229 protein ids per anchor, within the canonical window
    duf_of = collections.defaultdict(set)
    for line in open(DOMTBL):
        if line.startswith("#"):
            continue
        p = line.split()
        if p[4].split(".")[0] != "PF14012":
            continue
        d = p[-1].split("|")
        if len(d) < 4 or d[1] != "mqn":
            continue
        try:
            gi = int(d[2].split("=")[1])
            gap = int(d[3].split("=")[1])
        except (IndexError, ValueError):
            continue
        if in_window(gap, gi):
            duf_of[d[0]].add(p[0])

    # neighbour rows, for coordinates and strands
    rows = collections.defaultdict(list)
    contig_of = {}
    qcoord = {}
    for r in csv.DictReader(open(NEIGH), delimiter="\t"):
        if r["group"] != "mqn":
            continue
        rows[r["query_acc"]].append(r)
        contig_of[r["query_acc"]] = r["contig"]
        qcoord[r["query_acc"]] = (int(r["q_start"]), int(r["q_end"]), r["q_strand"])

    # anchors grouped by organism
    by_org = collections.defaultdict(lambda: {"MqnE": [], "MqnC": []})
    for q in rows:
        c = classif.get(q)
        if c in ("MqnE", "MqnC"):
            by_org[organism.get(q, "")][c].append(q)

    target_seqs = read_fasta(TARGETS)
    neigh_seqs = read_fasta(NEIGH_FAA)

    candidates = []
    for org, anchors in by_org.items():
        if not org or not anchors["MqnE"] or not anchors["MqnC"]:
            continue
        genus = org.split()[0]
        if genus in AVOID_GENERA:
            continue
        # MqnE anchors that carry DUF4229 in their window
        for e in anchors["MqnE"]:
            dufs = duf_of.get(e)
            if not dufs:
                continue
            c = anchors["MqnC"][0]
            duf = sorted(dufs)[0]
            if not (target_seqs.get(e) and target_seqs.get(c) and neigh_seqs.get(duf)):
                continue
            # arrangement of the DUF4229 gene relative to MqnE
            row = next((r for r in rows[e] if r["n_protein_id"] == duf), None)
            if row is None:
                continue
            qs, qe, qstr = qcoord[e]
            ns, ne = int(row["n_start"]), int(row["n_end"])
            overlap = max(0, min(qe, ne) - max(qs, ns) + 1)
            same = row["q_strand"] == row["n_strand"]
            try:
                gap = int(row["bp_gap"])
            except ValueError:
                gap = -1
            gi = int(row["gene_index"]) if row["gene_index"] not in ("", None) else 0
            # transcription-normalised offset: on the minus strand the coordinate
            # order runs opposite to the direction of transcription (NOTEBOOK 50)
            tgi = -gi if qstr == "-" else gi

            pref = PREFERRED_GENERA.index(genus) if genus in PREFERRED_GENERA else len(PREFERRED_GENERA)
            score = (pref,
                     0 if same else 1,
                     0 if regulator.get(e) == "Y" else 1,
                     gap if gap >= 0 else 99999)
            candidates.append({
                "score": score, "organism": org, "genus": genus,
                "phylum": phylum.get(e, ""), "contig": contig_of[e],
                "mqnE_acc": e, "mqnE_len": len(target_seqs[e]),
                "mqnC_acc": c, "mqnC_len": len(target_seqs[c]),
                "duf_acc": duf, "duf_len": len(neigh_seqs[duf]),
                "duf_gene_offset_transcription": tgi,
                "bp_gap": gap, "same_strand": "Y" if same else "N",
                "cds_overlap_bp": overlap,
                "asnc_lrp_regulator": regulator.get(e, "?"),
                "mqnE_seq": target_seqs[e], "mqnC_seq": target_seqs[c],
                "duf_seq": neigh_seqs[duf],
            })

    candidates.sort(key=lambda d: d["score"])
    # one entry per organism
    seen, ranked = set(), []
    for d in candidates:
        if d["organism"] in seen:
            continue
        seen.add(d["organism"])
        ranked.append(d)

    cols = ["rank", "organism", "genus", "phylum", "contig_accession",
            "mqnE_acc", "mqnE_len", "mqnC_acc", "mqnC_len", "duf4229_acc", "duf4229_len",
            "duf4229_offset_by_transcription", "bp_gap", "same_strand", "cds_overlap_bp",
            "asnc_lrp_regulator"]
    man = os.path.join(OUT, "construct_manifest.tsv")
    with open(man, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for i, d in enumerate(ranked, 1):
            w.writerow([i, d["organism"], d["genus"], d["phylum"], d["contig"],
                        d["mqnE_acc"], d["mqnE_len"], d["mqnC_acc"], d["mqnC_len"],
                        d["duf_acc"], d["duf_len"], d["duf_gene_offset_transcription"],
                        d["bp_gap"], d["same_strand"], d["cds_overlap_bp"],
                        d["asnc_lrp_regulator"]])

    # per-strain FASTA for the top strains, plus one combined file
    top = ranked[:10]
    combined = os.path.join(OUT, "top10_all_sequences.faa")
    with open(combined, "w") as allfh:
        for i, d in enumerate(top, 1):
            safe = re.sub(r"[^A-Za-z0-9]+", "_", d["organism"]).strip("_")
            path = os.path.join(OUT, "rank%02d_%s.faa" % (i, safe))
            with open(path, "w") as fh:
                for role, acc, seq in (("MqnE", d["mqnE_acc"], d["mqnE_seq"]),
                                       ("MqnC", d["mqnC_acc"], d["mqnC_seq"]),
                                       ("DUF4229", d["duf_acc"], d["duf_seq"])):
                    head = ">%s|%s|%s|%s" % (role, acc, d["organism"].replace(" ", "_"), d["contig"])
                    for fh2 in (fh, allfh):
                        fh2.write(head + "\n")
                        for j in range(0, len(seq), 60):
                            fh2.write(seq[j:j + 60] + "\n")

    print("=" * 78)
    print("CONSTRUCT SEQUENCE PACKAGE")
    print("=" * 78)
    print("Strains carrying MqnE + MqnC + DUF4229 with all three sequences: %d" % len(ranked))
    print("Excluded genera with a classical Men signature: %s" % ", ".join(sorted(AVOID_GENERA)))
    print()
    print("Top 10 by lab-friendliness, operon tightness and regulator presence:")
    print()
    print("%-4s %-42s %-9s %-8s %-6s %-5s" % ("rank", "organism", "MqnE", "DUF4229", "gap", "same"))
    for i, d in enumerate(top, 1):
        print("%-4d %-42s %-9s %-8s %-6s %-5s" % (
            i, d["organism"][:42], d["mqnE_acc"], d["duf_acc"][:8], d["bp_gap"], d["same_strand"]))
    print()
    print("Wrote %s" % os.path.relpath(man, ROOT))
    print("Wrote %s and %d per-strain FASTA files" % (os.path.relpath(combined, ROOT), len(top)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
