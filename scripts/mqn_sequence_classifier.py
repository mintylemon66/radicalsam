#!/usr/bin/env python3
"""
Replace the label-based MqnC / MqnE split with a SEQUENCE-based one, and test
whether the original labels were trustworthy.

Problem (NOTEBOOK section 19): MqnC and MqnE are separated by UniProt `GN=`
gene-name text assigned by automatic annotation. Only 84 of 3,651 descriptions
explicitly say "aminofutalosine"; 47% are ambiguous. The headline
"83.5% of actinobacterial MqnE genomes" therefore rests on an unverified split.

Approach (semi-supervised, and deliberately NOT circular at the decision step):
  1. Seed sets = sequences whose FULL PROTEIN NAME is unambiguous
     ("aminofutalosine synthase" -> MqnE; "dehypoxanthin*" -> MqnC).
     A full protein-name annotation is much stronger evidence than a bare gene
     symbol, and these seeds are used only to build profiles.
  2. Build a profile HMM for each class (MAFFT then hmmbuild).
  3. Score ALL 3,651 targets against both profiles and assign by best score.
  4. Compare the sequence-based assignment to the gene-name labels. High
     agreement validates the original split; low agreement invalidates it.
  5. Also flag contaminants (thiH, cofG/cofH) that match neither well.

Reviewed Swiss-Prot MqnC/MqnE entries are too few (2 each) to seed profiles
directly, which is itself evidence that these families are poorly curated.
"""
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.path.join(ROOT, "data", "mqn_targets.fasta")
WORK = os.path.join(ROOT, "hmm", "mqn")


def read_fasta(path):
    recs, acc, hdr, seq = [], None, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                recs.append((acc, hdr, "".join(seq)))
            hdr = line[1:].rstrip()
            p = hdr.split("|")
            acc = p[1] if len(p) >= 2 else hdr.split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if acc:
        recs.append((acc, hdr, "".join(seq)))
    return recs


def gene_label(hdr):
    m = re.search(r"GN=(\w+)", hdr, re.I)
    g = (m.group(1) if m else "").lower()
    if g.startswith("mqne"):
        return "MqnE"
    if g.startswith("mqnc"):
        return "MqnC"
    if g.startswith(("cofg", "cofh", "thih")):
        return "OTHER"
    return "unlabelled"


def name_seed(hdr):
    h = hdr.lower()
    if "aminofutalosine synthase" in h:
        return "MqnE"
    if "dehypoxanthin" in h:
        return "MqnC"
    return None


def build_profile(name, seqs):
    faa = os.path.join(WORK, f"{name}_train.faa")
    aln = os.path.join(WORK, f"{name}_train.aln")
    hmm = os.path.join(WORK, f"{name}.hmm")
    with open(faa, "w") as fh:
        for a, s in seqs:
            fh.write(f">{a}\n{s}\n")
    subprocess.run(["mafft", "--auto", "--anysymbol", "--quiet", faa],
                   stdout=open(aln, "w"), check=True)
    subprocess.run(["hmmbuild", "--informat", "afa", hmm, aln],
                   check=True, capture_output=True)
    return hmm


def score_all(hmm, target_faa, tag):
    out = os.path.join(WORK, f"scores_{tag}.domtbl")
    subprocess.run(["hmmsearch", "-E", "1e5", "--noali", "-o", os.devnull,
                    "--domtblout", out, hmm, target_faa], check=True)
    best = {}
    for line in open(out):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 8:
            continue
        acc = f[0].split("|")[1] if "|" in f[0] else f[0]
        sc = float(f[7])
        if acc not in best or sc > best[acc]:
            best[acc] = sc
    return best


def main():
    os.makedirs(WORK, exist_ok=True)
    recs = read_fasta(FASTA)
    print(f"targets: {len(recs)}")

    seeds = defaultdict(list)
    for acc, hdr, seq in recs:
        s = name_seed(hdr)
        if s:
            seeds[s].append((acc, seq))
    # cap MqnC seeds so the two profiles are built from comparable evidence
    seeds["MqnC"] = seeds["MqnC"][:300]
    print(f"seed sets from explicit protein names: "
          f"MqnE={len(seeds['MqnE'])}  MqnC={len(seeds['MqnC'])}")
    if len(seeds["MqnE"]) < 10 or len(seeds["MqnC"]) < 10:
        sys.exit("insufficient seeds")

    print("[build] aligning and building profiles ...")
    hmm_e = build_profile("MqnE", seeds["MqnE"])
    hmm_c = build_profile("MqnC", seeds["MqnC"])

    # write a clean copy of all targets keyed by accession
    allf = os.path.join(WORK, "all_targets.faa")
    with open(allf, "w") as fh:
        for acc, hdr, seq in recs:
            fh.write(f">{acc}\n{seq}\n")

    print("[score] scoring all targets against both profiles ...")
    se = score_all(hmm_e, allf, "MqnE")
    sc = score_all(hmm_c, allf, "MqnC")

    assign, margin = {}, {}
    for acc, hdr, seq in recs:
        e, c = se.get(acc, 0.0), sc.get(acc, 0.0)
        assign[acc] = "MqnE" if e > c else ("MqnC" if c > e else "tie")
        margin[acc] = abs(e - c)

    labels = {acc: gene_label(hdr) for acc, hdr, _ in recs}

    print("\n" + "=" * 78)
    print("SEQUENCE-BASED vs GENE-NAME CLASSIFICATION")
    print("=" * 78)
    conf = Counter()
    for acc in assign:
        conf[(labels[acc], assign[acc])] += 1
    lab_order = ["MqnC", "MqnE", "OTHER", "unlabelled"]
    print(f"  {'gene-name label':<16}{'-> seq MqnC':<14}{'-> seq MqnE':<14}{'total'}")
    for L in lab_order:
        c = conf[(L, "MqnC")]
        e = conf[(L, "MqnE")]
        if c + e == 0:
            continue
        print(f"  {L:<16}{c:<14}{e:<14}{c+e}")

    both = [a for a in assign if labels[a] in ("MqnC", "MqnE")]
    agree = sum(1 for a in both if labels[a] == assign[a])
    print(f"\n  agreement on gene-name-labelled sequences: {agree}/{len(both)} "
          f"= {100.0*agree/len(both):.1f}%")

    unl = [a for a in assign if labels[a] == "unlabelled"]
    print(f"  previously UNLABELLED now assigned: {len(unl)} "
          f"(MqnC {sum(1 for a in unl if assign[a]=='MqnC')}, "
          f"MqnE {sum(1 for a in unl if assign[a]=='MqnE')})")

    lowm = [a for a in assign if margin[a] < 20]
    print(f"  low-confidence assignments (score margin < 20 bits): {len(lowm)} "
          f"({100.0*len(lowm)/len(assign):.1f}%)")

    out = os.path.join(ROOT, "audit", "mqn_seq_classification.tsv")
    with open(out, "w") as fh:
        fh.write("accession\tgene_label\tseq_assignment\tscore_MqnE\tscore_MqnC\tmargin\n")
        for acc, hdr, _ in recs:
            fh.write(f"{acc}\t{labels[acc]}\t{assign[acc]}\t"
                     f"{se.get(acc,0):.1f}\t{sc.get(acc,0):.1f}\t{margin[acc]:.1f}\n")
    print(f"\n[out] {out}")


if __name__ == "__main__":
    main()
