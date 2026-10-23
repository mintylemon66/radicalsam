#!/usr/bin/env python3
"""
Cross-family specificity of DUF4229 (PF14012), Methods Section 2.9.

Question: is DUF4229 a partner of MqnE specifically, or is it a generic
neighbour of radical SAM enzymes in general? The test compares five anchor
families that all belong to the radical SAM superfamily:

    MqnE   aminofutalosine synthase          (futalosine pathway)
    MqnC   cyclic dehypoxanthinyl futalosine synthase (futalosine pathway)
    bg     radical SAM background            (superfamily, non futalosine)
    BioB   biotin synthase                   (unrelated pathway)
    LipA   lipoyl synthase                   (unrelated pathway)

Because DUF4229 is nearly restricted to Actinomycetota, the comparison that
matters is the taxonomy-matched one: actinobacterial genomes only. A pooled
test of actinobacterial MqnE against the union of the three non-futalosine
actinobacterial anchor sets is reported as the headline specificity number.

This script rebuilds audit/duf4229_cross_family_specificity.tsv from raw data.
It reads nothing that is itself a summary of the analysis.

Inputs (all raw):
  audit/pfam_wide.domtbl              full Pfam-A scan of mqn + bg neighbourhoods
  audit/pfam_wide_generalize.domtbl   full Pfam-A scan of BioB + LipA neighbourhoods
  audit/neighbors_seq.tsv             reconstructed mqn + bg neighbourhoods
  audit/generalize_neighbors_seq.tsv  reconstructed BioB + LipA neighbourhoods
  audit/mqn_seq_classification.tsv    two-way sequence classifier, MqnE vs MqnC
  data/metadata_org.tsv               organism and phylum per UniProt accession

One note on the window. The project shorthand for the operonic window is
"bp_gap <= 8000 and abs(gene_index) <= 10 and gene_index != 0". The last clause
is a proxy for the rule Methods 2.9 actually states, which is that the anchor
CDS is excluded from its own neighbourhood. The two are not identical: the
neighbourhood reconstruction indexed genes relative to the FIRST CDS overlapping
the query coordinates, so in the 1,160 genomes where the upstream gene overlaps
the query start by a few bases, the anchor sits at gene_index 1 and gene_index 0
is a genuine neighbour. This script excludes the anchor by identity, which is
the stated rule, and reports how many genomes that distinction affects. For
PF14012 it is 2 genomes, and they are part of the canonical 342.

Usage:
  python3 scripts/cross_family_specificity.py            print only
  python3 scripts/cross_family_specificity.py --write    also rewrite the TSV

Exit code 0 if every canonical value reproduces, 1 otherwise. With --write the
TSV is rewritten only when all checks pass.
"""
import csv
import os
import re
import sys
from collections import defaultdict

try:
    from scipy.stats import fisher_exact
except ImportError:
    sys.exit("scipy is required: python3 -m pip install scipy")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")
DATA = os.path.join(ROOT, "data")

DOMTBL_MQN = os.path.join(AUD, "pfam_wide.domtbl")
DOMTBL_GEN = os.path.join(AUD, "pfam_wide_generalize.domtbl")
NEIGH_MQN = os.path.join(AUD, "neighbors_seq.tsv")
NEIGH_GEN = os.path.join(AUD, "generalize_neighbors_seq.tsv")
CLASSIFY = os.path.join(AUD, "mqn_seq_classification.tsv")
META = os.path.join(DATA, "metadata_org.tsv")
OUT_TSV = os.path.join(AUD, "duf4229_cross_family_specificity.tsv")

# ---- canonical conventions -------------------------------------------------
PFAM_ACC = "PF14012"          # DUF4229
PFAM_GA = 27.0                # Pfam gathering threshold for PF14012
MAX_GAP = 8000                # intergenic gap, bp
MAX_GENE_OFFSET = 10          # absolute gene offset from the anchor CDS
ACTINO_RE = re.compile(r"actinomycetota|actinobacteri", re.I)

# ---- canonical values this script must reproduce ---------------------------
CANON_TABLE = {
    # family: (all_denom, all_pos, actino_denom, actino_pos)
    "MqnE": (1733, 284, 329, 280),
    "MqnC": (1915, 58, 566, 55),
    "background_radicalSAM": (1211, 0, 212, 0),
    "BioB_biotin_synthase": (1345, 0, 259, 0),
    "LipA_lipoyl_synthase": (1439, 0, 300, 0),
}
CANON_DERIVED = {
    "futalosine_all_pos": 342,
    "futalosine_all_denom": 3648,
    "actino_MqnE_pct": 85.1,
    "actino_MqnC_pct": 9.7,
    "pooled_nonMqnE_pos": 0,
    "pooled_nonMqnE_denom": 771,
    # stated in Methods 2.9, over the futalosine plus background genomes
    "methods_tax_coverage_pct": 97.1,
}
CANON_P = {
    "BioB_biotin_synthase": 4.9e-117,
    "LipA_lipoyl_synthase": 5.2e-128,
    "background_radicalSAM": 4.7e-103,
    "pooled": 3.0e-211,
}

FAMILY_ORDER = ["MqnE", "MqnC", "background_radicalSAM",
                "BioB_biotin_synthase", "LipA_lipoyl_synthase"]
GROUP_OF_FAMILY = {
    "background_radicalSAM": "bg",
    "BioB_biotin_synthase": "BioB",
    "LipA_lipoyl_synthase": "LipA",
}


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------
def load_taxonomy():
    """UniProt accession -> phylum string, parsed from sp|ACC|NAME or tr|ACC|NAME."""
    tax = {}
    with open(META, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r.get("Sequence_ID", "") or ""
            parts = sid.split("|")
            acc = parts[1] if len(parts) >= 2 else sid
            tax[acc] = (r.get("Phylum") or "").strip()
    return tax


def load_classifier():
    """UniProt accession -> MqnE or MqnC, from the two-way sequence classifier."""
    call = {}
    with open(CLASSIFY, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            call[r["accession"]] = r["seq_assignment"]
    return call


def load_anchor_cds(paths):
    """
    query accession -> protein id of the anchor CDS itself.

    Methods 2.9 excludes the anchor CDS from its own neighbourhood. The usual
    proxy for that is gene_index != 0, but the neighbourhood reconstruction
    anchored the index on the FIRST CDS overlapping the query coordinates, and
    when the gene immediately upstream overlaps the query start by a few bases
    that first CDS is the upstream neighbour rather than the anchor. In those
    genomes the anchor sits at gene_index 1 and gene_index 0 is a real
    neighbour. Identifying the anchor by coordinate overlap instead of by index
    implements the Methods rule exactly and is index-shift proof.
    """
    anchor = {}
    anchor_gi = {}
    best = {}
    for path in paths:
        with open(path, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            ix = {c: i for i, c in enumerate(hdr)}
            for f in rd:
                try:
                    acc = f[ix["query_acc"]]
                    qs, qe = int(f[ix["q_start"]]), int(f[ix["q_end"]])
                    ns, ne = int(f[ix["n_start"]]), int(f[ix["n_end"]])
                except (ValueError, IndexError):
                    continue
                overlap = min(qe, ne) - max(qs, ns)
                gi = f[ix["gene_index"]]
                try:
                    gi = int(gi)
                except (TypeError, ValueError):
                    gi = 10 ** 9
                key = (overlap, -abs(gi))
                if acc not in best or key > best[acc]:
                    best[acc] = key
                    anchor[acc] = f[ix["n_protein_id"]]
                    anchor_gi[acc] = gi
    off_by_one = sum(1 for g in anchor_gi.values() if g != 0)
    return anchor, off_by_one


def parse_desc(desc):
    """
    Split the domtbl description field, which encodes
    query_accession|group|gi=<gene_index>|gap=<bp>.
    Returns (accession, group, gene_index, bp_gap) or None.
    """
    p = desc.split("|")
    if len(p) < 4:
        return None
    try:
        gi = int(p[2].split("=", 1)[1])
        gap = int(p[3].split("=", 1)[1])
    except (ValueError, IndexError):
        return None
    return p[0], p[1], gi, gap


def scan_domtbl(path, anchor):
    """
    One pass over a Pfam-A domtbl.

    Returns (universe, positives, skipped_anchor, index_zero) where
      universe[group]  = every anchor genome represented in the file
      positives[group] = genomes with a PF14012 hit inside the canonical window
      skipped_anchor   = number of PF14012 hits discarded as the anchor CDS
      index_zero[group]= genomes whose only qualifying hit sits at gene_index 0
                         on a CDS that is NOT the anchor, which the shorthand
                         rule gene_index != 0 would wrongly drop
    """
    universe = defaultdict(set)
    positives = defaultdict(set)
    nonzero = defaultdict(set)
    skipped_anchor = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            parsed = parse_desc(f[22])
            if parsed is None:
                continue
            acc, group, gi, gap = parsed
            universe[group].add(acc)

            fam_acc = f[4]
            fam = fam_acc.split(".")[0] if fam_acc != "-" else f[3]
            if fam != PFAM_ACC:
                continue
            # the scan was run with --cut_ga, this re-asserts the same threshold
            try:
                seq_score = float(f[7])
                dom_score = float(f[13])
            except ValueError:
                continue
            if seq_score < PFAM_GA or dom_score < PFAM_GA:
                continue
            # canonical window
            if gap > MAX_GAP or abs(gi) > MAX_GENE_OFFSET:
                continue
            if f[0] == anchor.get(acc):
                skipped_anchor += 1
                continue
            positives[group].add(acc)
            if gi != 0:
                nonzero[group].add(acc)
    index_zero = {g: positives[g] - nonzero[g] for g in positives}
    return universe, positives, skipped_anchor, index_zero


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def pct(k, n):
    return (100.0 * k / n) if n else 0.0


def frac(k, n):
    return f"{k}/{n} = {pct(k, n):.2f}%"


def fisher_greater(k1, n1, k2, n2):
    _, p = fisher_exact([[k1, n1 - k1], [k2, n2 - k2]], alternative="greater")
    return p


def odds_ratio(k1, n1, k2, n2):
    """Haldane-Anscombe corrected odds ratio."""
    a, b = k1 + 0.5, n1 - k1 + 0.5
    c, d = k2 + 0.5, n2 - k2 + 0.5
    return (a / b) / (c / d)


def check(label, got, want, tol=0):
    ok = abs(got - want) <= tol
    return ok, f"  [{'PASS' if ok else 'FAIL'}] {label:<52} got {got!r}  expected {want!r}"


# ---------------------------------------------------------------------------
def main():
    write = "--write" in sys.argv[1:]

    for p in (DOMTBL_MQN, DOMTBL_GEN, NEIGH_MQN, NEIGH_GEN, CLASSIFY, META):
        if not os.path.exists(p):
            sys.exit(f"missing required input: {p}")

    print("=" * 78)
    print("CROSS-FAMILY SPECIFICITY OF DUF4229 (PF14012), Methods 2.9")
    print("=" * 78)
    print(f"  window        : bp_gap <= {MAX_GAP} and abs(gene_index) <= {MAX_GENE_OFFSET},")
    print("                  anchor CDS excluded")
    print(f"  detection     : HMMER hmmsearch --cut_ga vs Pfam-A, {PFAM_ACC} GA = {PFAM_GA:.0f}")
    print("  MqnE vs MqnC  : audit/mqn_seq_classification.tsv, column seq_assignment")
    print("  actinobacteria: regex /actinomycetota|actinobacteri/i on data/metadata_org.tsv")

    tax = load_taxonomy()
    call = load_classifier()
    print(f"\n[load] taxonomy rows {len(tax)}, classifier calls {len(call)}")

    anchor, off_by_one = load_anchor_cds([NEIGH_MQN, NEIGH_GEN])
    print(f"[load] anchor CDS identified by coordinate overlap for "
          f"{len(anchor)} genomes")
    print(f"[load] genomes whose anchor is NOT at gene_index 0: "
          f"{frac(off_by_one, len(anchor))}")

    uni_a, pos_a, skip_a, iz_a = scan_domtbl(DOMTBL_MQN, anchor)
    uni_b, pos_b, skip_b, iz_b = scan_domtbl(DOMTBL_GEN, anchor)
    universe = dict(uni_a)
    universe.update(uni_b)
    positives = dict(pos_a)
    positives.update(pos_b)
    index_zero = dict(iz_a)
    index_zero.update(iz_b)
    print(f"[scan] {os.path.basename(DOMTBL_MQN)}: groups "
          + ", ".join(f"{g} n={len(universe[g])}" for g in sorted(uni_a)))
    print(f"[scan] {os.path.basename(DOMTBL_GEN)}: groups "
          + ", ".join(f"{g} n={len(universe[g])}" for g in sorted(uni_b)))
    print(f"[scan] {PFAM_ACC} hits discarded as the anchor CDS itself: {skip_a + skip_b}")
    iz_all = sorted(a for s in index_zero.values() for a in s)
    print(f"[scan] genomes kept only because the anchor was located by coordinate")
    print(f"       rather than assumed to be gene_index 0: {len(iz_all)}"
          + (f" ({', '.join(iz_all)})" if iz_all else ""))

    # ---- split the futalosine set into MqnE and MqnC -----------------------
    denom, pos = {}, {}
    mqn_all = universe.get("mqn", set())
    mqn_pos = positives.get("mqn", set())
    unclassified = sorted(a for a in mqn_all if call.get(a) not in ("MqnE", "MqnC"))
    denom["MqnE"] = {a for a in mqn_all if call.get(a) == "MqnE"}
    denom["MqnC"] = {a for a in mqn_all if call.get(a) == "MqnC"}
    pos["MqnE"] = {a for a in mqn_pos if call.get(a) == "MqnE"}
    pos["MqnC"] = {a for a in mqn_pos if call.get(a) == "MqnC"}
    for fam, grp in GROUP_OF_FAMILY.items():
        denom[fam] = universe.get(grp, set())
        pos[fam] = positives.get(grp, set())
    print(f"[split] futalosine anchors unclassified by the two-way classifier: "
          f"{len(unclassified)}")

    def is_actino(acc):
        return bool(ACTINO_RE.search(tax.get(acc, "")))

    def mapped(accs):
        return sum(1 for a in accs if tax.get(a))

    for fam in FAMILY_ORDER:
        d = denom[fam]
        print(f"[tax] {fam:<24} phylum assigned to {frac(mapped(d), len(d))}")
    core = universe.get("mqn", set()) | universe.get("bg", set())
    core_pct = pct(mapped(core), len(core))
    print(f"[tax] futalosine plus background (the Methods 2.9 figure): "
          f"{frac(mapped(core), len(core))}")

    rows = []
    for fam in FAMILY_ORDER:
        d, p = denom[fam], pos[fam]
        ad = {a for a in d if is_actino(a)}
        ap = {a for a in p if is_actino(a)}
        rows.append((fam, len(d), len(p), len(ad), len(ap)))

    # ---- per-family table ---------------------------------------------------
    print("\n" + "=" * 78)
    print("PER-FAMILY DUF4229 CO-OCCURRENCE")
    print("=" * 78)
    print(f"  {'family':<24}{'all genomes':<26}{'actinobacteria only':<26}")
    for fam, dn, pn, adn, apn in rows:
        print(f"  {fam:<24}{frac(pn, dn):<26}{frac(apn, adn):<26}")

    fut_pos = rows[0][2] + rows[1][2]
    fut_den = rows[0][1] + rows[1][1]
    print(f"\n  futalosine set pooled (MqnE + MqnC)      {frac(fut_pos, fut_den)}")

    # ---- Fisher, actinobacteria only ---------------------------------------
    by_fam = {r[0]: r for r in rows}
    e_pos, e_den = by_fam["MqnE"][4], by_fam["MqnE"][3]

    print("\n" + "=" * 78)
    print("ONE-SIDED FISHER, ACTINOBACTERIAL GENOMES ONLY, EACH FAMILY vs MqnE")
    print("=" * 78)
    print(f"  reference: MqnE {frac(e_pos, e_den)}")
    print(f"  {'family':<24}{'actino DUF4229+':<26}{'p':<14}{'odds ratio'}")
    pvals = {}
    for fam in FAMILY_ORDER[1:]:
        _, _, _, adn, apn = by_fam[fam]
        p = fisher_greater(e_pos, e_den, apn, adn)
        pvals[fam] = p
        orr = odds_ratio(e_pos, e_den, apn, adn)
        print(f"  {fam:<24}{frac(apn, adn):<26}{p:<14.1e}{orr:>10.1f}x")

    pool_pos = sum(by_fam[f][4] for f in GROUP_OF_FAMILY)
    pool_den = sum(by_fam[f][3] for f in GROUP_OF_FAMILY)
    p_pool = fisher_greater(e_pos, e_den, pool_pos, pool_den)
    or_pool = odds_ratio(e_pos, e_den, pool_pos, pool_den)
    print("\n  POOLED non-MqnE actinobacterial radical SAM (background + BioB + LipA)")
    print(f"    MqnE          {frac(e_pos, e_den)}")
    print(f"    pooled non-MqnE {frac(pool_pos, pool_den)}")
    print(f"    Fisher one-sided p = {p_pool:.1e}, odds ratio = {or_pool:.1f}x")

    # The pooled denominator sums the three anchor sets, which is what the
    # canonical 771 is. A handful of UniProt accessions serve as anchors in more
    # than one set, so the number of DISTINCT actinobacterial genomes is smaller.
    # Reported here so the double counting is visible rather than hidden.
    pool_union = set()
    for fam in GROUP_OF_FAMILY:
        pool_union |= {a for a in denom[fam] if is_actino(a)}
    dup = pool_den - len(pool_union)
    print(f"    note: {pool_den} is the sum of the three anchor sets; "
          f"{len(pool_union)} distinct genomes, {dup} counted in two sets")

    # ---- PASS / FAIL --------------------------------------------------------
    print("\n" + "=" * 78)
    print("PASS / FAIL AGAINST CANONICAL VALUES")
    print("=" * 78)
    results = []
    for fam, dn, pn, adn, apn in rows:
        cd, cp, cad, cap = CANON_TABLE[fam]
        results.append(check(f"{fam} all-genome denominator", dn, cd))
        results.append(check(f"{fam} all-genome DUF4229 positive", pn, cp))
        results.append(check(f"{fam} actinobacterial denominator", adn, cad))
        results.append(check(f"{fam} actinobacterial DUF4229 positive", apn, cap))

    results.append(check("futalosine DUF4229 positive (MqnE + MqnC)",
                         fut_pos, CANON_DERIVED["futalosine_all_pos"]))
    results.append(check("futalosine denominator (MqnE + MqnC)",
                         fut_den, CANON_DERIVED["futalosine_all_denom"]))
    results.append(check("actinobacterial MqnE rate percent",
                         round(pct(e_pos, e_den), 1),
                         CANON_DERIVED["actino_MqnE_pct"], tol=0.05))
    results.append(check("actinobacterial MqnC rate percent",
                         round(pct(by_fam["MqnC"][4], by_fam["MqnC"][3]), 1),
                         CANON_DERIVED["actino_MqnC_pct"], tol=0.05))
    results.append(check("pooled non-MqnE actinobacterial positive",
                         pool_pos, CANON_DERIVED["pooled_nonMqnE_pos"]))
    results.append(check("pooled non-MqnE actinobacterial denominator",
                         pool_den, CANON_DERIVED["pooled_nonMqnE_denom"]))
    results.append(check("phylum coverage percent, futalosine plus background",
                         round(core_pct, 1), CANON_DERIVED["methods_tax_coverage_pct"],
                         tol=0.05))

    for fam in ("BioB_biotin_synthase", "LipA_lipoyl_synthase",
                "background_radicalSAM"):
        got = float(f"{pvals[fam]:.1e}")
        results.append(check(f"Fisher p, MqnE vs {fam}", got, CANON_P[fam],
                             tol=abs(CANON_P[fam]) * 0.02))
    results.append(check("Fisher p, pooled", float(f"{p_pool:.1e}"),
                         CANON_P["pooled"], tol=abs(CANON_P["pooled"]) * 0.02))

    for ok, line in results:
        print(line)
    n_fail = sum(1 for ok, _ in results if not ok)
    print("-" * 78)
    if n_fail:
        print(f"  OVERALL: FAIL, {n_fail} of {len(results)} checks did not match")
    else:
        print(f"  OVERALL: PASS, all {len(results)} checks reproduce the canonical values")

    # ---- output -------------------------------------------------------------
    header = [
        "# DUF4229 (PF14012) co-occurrence across separate radical SAM gene families",
        "# Detection: HMMER PF14012 at Pfam gathering threshold (GA=27) on reconstructed neighbourhoods.",
        "# Sources: audit/pfam_wide.domtbl (MqnE/MqnC/bg), audit/pfam_wide_generalize.domtbl (BioB/LipA).",
        "# window = operonic +/-8 kb (bp_gap <= 8000); actino = Actinomycetota per data/metadata_org.tsv",
    ]
    body = ["\t".join(["family", "all_genomes_denom", "all_genomes_DUF4229pos",
                       "all_rate_pct", "actino_denom", "actino_DUF4229pos",
                       "actino_rate_pct"])]
    for fam, dn, pn, adn, apn in rows:
        body.append("\t".join([fam, str(dn), str(pn), f"{pct(pn, dn):.2f}",
                               str(adn), str(apn), f"{pct(apn, adn):.2f}"]))
    footer = [
        "# Fisher (actino-only, one-sided): "
        f"MqnE vs BioB p={pvals['BioB_biotin_synthase']:.1e}; "
        f"vs LipA p={pvals['LipA_lipoyl_synthase']:.1e}; "
        f"vs bg p={pvals['background_radicalSAM']:.1e}",
        f"# Pooled: MqnE {e_pos}/{e_den} vs {pool_pos}/{pool_den} actinobacterial "
        f"non-MqnE radical SAM (bg+BioB+LipA), p={p_pool:.1e}",
    ]
    text = "\n".join(header + body + footer) + "\n"

    print("\n" + "=" * 78)
    print(f"TABLE AS IT WOULD BE WRITTEN TO {os.path.relpath(OUT_TSV, ROOT)}")
    print("=" * 78)
    print(text, end="")

    if write:
        if n_fail:
            sys.exit("\nrefusing to write: canonical checks failed")
        with open(OUT_TSV, "w") as fh:
            fh.write(text)
        print(f"\n[write] rewrote {OUT_TSV}")
    else:
        print("\n[dry run] nothing written, pass --write to rewrite the TSV")

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
