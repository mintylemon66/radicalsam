#!/usr/bin/env python3
"""
Generalized text-recall sweep: which neighbour families does product-name
annotation BURY? (Methods 2.10, "Finding families hidden by annotation".)

scripts/hmm_recount.py answers this for DUF4229 alone. Until now there was no
saved script for the generalized sweep, so audit/buried_partners_clade_specific.tsv
could not be regenerated from source. This is that script.

WHAT IT DOES
  1. Sweeps EVERY Pfam family present in audit/pfam_wide.domtbl and counts, in
     the canonical operonic window, how many futalosine (target) and how many
     background anchor genomes carry at least one member.
  2. Computes the TEXT-RECALL metric, defined exactly as hmm_recount.py defines
     it for DUF4229: of all proteins the HMM calls members of the family (at the
     Pfam gathering threshold, the ground truth), what fraction would a
     case-insensitive product-name search for the family name have found?
     buried_pct = 100 - text_recall_pct.
     The metric is validated by DUF4229 re-deriving at exactly 56.9%.
  3. Emits the five flagged dark families with the columns of the existing TSV.
  4. Recomputes the clade-matched controls (Actinomycetota and Bacillota) with
     one-sided Fisher's exact tests.
  5. Prints PASS/FAIL against every canonical number of the study.

Read-only unless --write is given.

Usage:
    python3 scripts/text_recall_sweep.py                 # report only
    python3 scripts/text_recall_sweep.py --write         # also rewrite the TSVs
    python3 scripts/text_recall_sweep.py --exclude-index-zero   # variant, see below
"""
import argparse
import csv
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "audit")

DOM_MQN = os.path.join(AUD, "pfam_wide.domtbl")             # mqn + bg neighbourhoods
DOM_GEN = os.path.join(AUD, "pfam_wide_generalize.domtbl")  # BioB + LipA neighbourhoods
TSV_MQN = os.path.join(AUD, "neighbors_seq.tsv")
TSV_GEN = os.path.join(AUD, "generalize_neighbors_seq.tsv")
CLASSIFY = os.path.join(AUD, "mqn_seq_classification.tsv")  # TWO-WAY classifier
META = os.path.join(ROOT, "data", "metadata_org.tsv")

OUT_FLAGGED = os.path.join(AUD, "buried_partners_clade_specific.tsv")
OUT_SWEEP = os.path.join(AUD, "text_recall_sweep_all_families.tsv")

# ---------------------------------------------------------------------------
# Canonical window
# ---------------------------------------------------------------------------
MAX_GAP = 8000
MAX_GENE_INDEX = 10

# On the "gene_index != 0" clause of the window convention.
# -----------------------------------------------------------------------------
# The convention is written as: gap <= 8000 bp AND abs(gene_index) <= 10 AND
# gene_index != 0, the last clause meaning "exclude the anchor CDS itself".
# In this project's neighbour tables the anchor CDS is never a row at all:
# 0 of the 301,922 rows of audit/neighbors_seq.tsv have neighbour coordinates
# equal to the anchor coordinates. gene_index == 0 is instead carried by 6,798
# rows that are genuine flanking genes which OVERLAP the anchor CDS by a few
# bases (which is also why their bp_gap is 0). Discarding them therefore throws
# away real neighbours rather than the anchor.
# The published counts keep them. With gene_index == 0 rows kept, DUF4229 is
# 342/3,648 target-wide and 280/329 in actinobacterial MqnE, the canonical
# values. With a literal gene_index != 0 filter it is 340/3,648 and 278/329,
# and the same-strand denominator drops from 565/566 to 561/562, which is
# exactly the "difference of four rows whose origin is not yet identified"
# recorded in Methods 2.10. That origin is identified here.
# Default therefore keeps gene_index == 0. Pass --exclude-index-zero to see the
# literal-filter variant; it fails the 342 and 280 checks by construction.

ACTINO_RE = re.compile(r"actinomycetota|actinobacteri", re.I)
FIRMICUTE = "Bacillota"

# The five families flagged as annotation-buried partners, in the order the
# existing TSV lists them.
FLAGGED = [("DUF4229", "PF14012"),
           ("YtxC", "PF08812"),
           ("DUF1450", "PF07293"),
           ("YuiB", "PF14068"),
           ("DUF1462", "PF07315")]

# A family's name is only usable as a text query when the name IS the identifier
# an annotator would write. Families named for a function score near zero recall
# trivially, because GenBank records the function, not the Pfam identifier. This
# heuristic ("DUF" numbers and Y-gene symbols) is used only to rank the printed
# sweep, never to select the flagged five, which are the published set.
DARK_NAME_RE = re.compile(r"^(DUF\d+|Y[a-z]{2,4}[A-Z]?\d?)$")

MIN_GENOMES_FOR_SWEEP = 15   # same floor as scripts/pfam_wide_scan.py

# Dominance threshold for the dominant_phylum column. A family whose top phylum
# holds less than this share of its mapped member proteins is reported "mixed".
PHYLUM_DOMINANCE = 95.0


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def in_window(gene_index, bp_gap, exclude_index_zero):
    if bp_gap is None or gene_index is None:
        return False
    if bp_gap > MAX_GAP or abs(gene_index) > MAX_GENE_INDEX:
        return False
    if exclude_index_zero and gene_index == 0:
        return False
    return True


def parse_domtbl(path, keep_accessions=None):
    """Parse an hmmsearch --domtblout produced against the neighbour FASTA.

    The FASTA description encodes query_accession|group|gi=<gene_index>|gap=<bp>,
    so the domtbl alone carries the anchor genome, the group and the window
    coordinates. Returns:
      fam_name[pfam_acc]                 -> family name (the HMM query name)
      fam_proteins[pfam_acc]             -> set of member protein ids, NO window
                                            filter (this is the text-recall
                                            denominator)
      fam_rows[pfam_acc]                 -> list of (protein, anchor, group,
                                            gene_index, bp_gap), one per domain
                                            hit line, NO window filter
    """
    fam_name, fam_proteins, fam_rows = {}, defaultdict(set), defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            acc = f[4].split(".")[0] if f[4] != "-" else f[3]
            if keep_accessions is not None and acc not in keep_accessions:
                continue
            desc = f[22]
            parts = desc.split("|")
            if len(parts) < 4:
                continue
            try:
                anchor, group = parts[0], parts[1]
                gi = int(parts[2].split("=", 1)[1])
                gap = int(parts[3].split("=", 1)[1])
            except (IndexError, ValueError):
                continue
            fam_name[acc] = f[3]
            fam_proteins[acc].add(f[0])
            fam_rows[acc].append((f[0], anchor, group, gi, gap))
    return fam_name, fam_proteins, fam_rows


def load_neighbor_table(path, want_products=True):
    """Returns (group_of_anchor, products_of_protein, rows_by_anchor_protein)."""
    group_of = {}
    products = defaultdict(set)
    pair_rows = defaultdict(list)   # (anchor, protein) -> [(q_strand, n_strand, bp_gap)]
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            anchor = r["query_acc"]
            group_of[anchor] = r["group"]
            pid = r["n_protein_id"]
            if want_products:
                products[pid].add(r["n_product"] or "")
            try:
                gap = int(r["bp_gap"])
            except (TypeError, ValueError):
                gap = None
            pair_rows[(anchor, pid)].append((r["q_strand"], r["n_strand"], gap))
    return group_of, products, pair_rows


def load_taxonomy():
    """UniProt accession -> phylum. Sequence_ID looks like sp|ACC|NAME."""
    tax = {}
    with open(META, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sid = r.get("Sequence_ID", "")
            parts = sid.split("|")
            acc = parts[1] if len(parts) >= 2 else sid
            tax[acc] = r.get("Phylum", "") or ""
    return tax


def load_classification():
    """UniProt accession -> MqnE or MqnC, from the TWO-WAY sequence classifier."""
    out = {}
    with open(CLASSIFY, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["accession"]] = r["seq_assignment"]
    return out


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
def text_recall(family_name, member_proteins, products):
    """The metric, defined exactly as scripts/hmm_recount.py defines it.

    Ground truth is the HMM call at the Pfam gathering threshold. The query is a
    case-insensitive literal search for the family name in the GenBank product
    name (column n_product of audit/neighbors_seq.tsv).

      numerator   = members found by BOTH the text search and the HMM
      denominator = ALL HMM-detected members, with no window restriction and
                    across both groups

    hmm_recount.py computes len(prot_label & prot_hmm) / len(prot_hmm) for
    DUF4229 and gets 255/448 = 56.9%. This reproduces that number exactly.
    """
    needle = family_name.lower()
    both = 0
    for pid in member_proteins:
        for prod in products.get(pid, ()):
            if needle in prod.lower():
                both += 1
                break
    denom = len(member_proteins)
    pct = 100.0 * both / denom if denom else 0.0
    return both, denom, pct


def dominant_phylum(member_proteins, protein_anchors, tax):
    """Dominant phylum over ALL HMM-detected members, no window filter.

    Each member protein is assigned the phylum of the anchor genome or genomes
    whose neighbourhood it was pulled from; a protein seen in two phyla counts
    once for each. This reproduces the published Actinomycetota 442/448 for
    DUF4229 and Bacillota 247/259 for YtxC.
    """
    counts = Counter()
    for pid in member_proteins:
        for ph in {tax.get(a, "") or "?" for a in protein_anchors.get(pid, ())}:
            counts[ph] += 1
    total = sum(counts.values())
    if not total:
        return "unmapped", 0, 0
    named = [(ph, n) for ph, n in counts.items() if ph != "?"]
    if not named:
        return "unmapped", 0, total
    ph, n = max(named, key=lambda kv: kv[1])
    return ph, n, total


def operon_stats(window_rows, pair_rows):
    """Same-strand fraction and median intergenic gap.

    COUNTING CONVENTION (stated explicitly because it sets the denominator).
    The unit is one (anchor genome, neighbour protein) pair that clears the
    window, expanded over every row the neighbour table holds for that pair.
    Two consequences worth naming:
      - the domtbl is de-duplicated first, so a protein matched by two domains
        of the same family contributes once, not twice;
      - audit/neighbors_seq.tsv lists 85,419 of its 215,827 (anchor, protein)
        keys more than once, because neighbourhoods were reconstructed from
        overlapping efetch records, so such a pair contributes as many rows as
        the table holds for it.
    This is the convention that reproduces the published DUF4229 565/566 and
    YtxC 303/303 exactly. Counting unique pairs instead gives the same
    percentages on smaller denominators (DUF4229 341/342, YtxC 252/252).
    The median gap is taken over the same expanded row set.
    """
    pairs = sorted({(anchor, pid) for pid, anchor, _g, _gi, _gap in window_rows})
    same = total = 0
    gaps = []
    for anchor, pid in pairs:
        for q_strand, n_strand, gap in pair_rows.get((anchor, pid), ()):
            total += 1
            if q_strand and q_strand == n_strand:
                same += 1
            if gap is not None:
                gaps.append(gap)
    med = statistics.median(gaps) if gaps else None
    return same, total, med


def fmt_frac(num, den):
    if not den:
        return f"{num}/{den} = n/a"
    return f"{num}/{den} = {100.0 * num / den:.2f}%"


def fisher(k1, n1, k2, n2):
    """One-sided Fisher (greater) plus a Haldane-Anscombe corrected odds ratio."""
    from scipy.stats import fisher_exact
    _, p = fisher_exact([[k1, n1 - k1], [k2, n2 - k2]], alternative="greater")
    a, b = k1 + 0.5, n1 - k1 + 0.5
    c, d = k2 + 0.5, n2 - k2 + 0.5
    return p, (a / b) / (c / d)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true",
                    help="rewrite audit/buried_partners_clade_specific.tsv and the "
                         "full-sweep TSV; without it nothing is written")
    ap.add_argument("--exclude-index-zero", action="store_true",
                    help="apply gene_index != 0 literally; see the note in the "
                         "source. Fails the 342 and 280 checks by construction.")
    ap.add_argument("--top", type=int, default=20,
                    help="how many most-buried dark families to print (default 20)")
    ap.add_argument("--out-flagged", default=OUT_FLAGGED)
    ap.add_argument("--out-sweep", default=OUT_SWEEP)
    args = ap.parse_args()

    xz = args.exclude_index_zero
    for path in (DOM_MQN, DOM_GEN, TSV_MQN, TSV_GEN, CLASSIFY, META):
        if not os.path.exists(path):
            sys.exit(f"missing required input: {path}")

    print("=" * 78)
    print("TEXT-RECALL SWEEP: families buried by product-name annotation")
    print("=" * 78)
    print(f"  window          : intergenic gap <= {MAX_GAP} bp and "
          f"abs(gene_index) <= {MAX_GENE_INDEX}")
    print(f"  gene_index == 0 : {'DROPPED (literal filter)' if xz else 'KEPT (see source note)'}")
    print(f"  detection       : hmmsearch --cut_ga vs Pfam-A, {os.path.basename(DOM_MQN)}")
    print(f"  classifier      : two-way {os.path.basename(CLASSIFY)} (seq_assignment)")

    # ---- inputs -----------------------------------------------------------
    fam_name, fam_prot, fam_rows = parse_domtbl(DOM_MQN)
    print(f"\n[load] {len(fam_name)} Pfam families with at least one hit in "
          f"{os.path.basename(DOM_MQN)}")

    group_of, products, pair_rows = load_neighbor_table(TSV_MQN)
    n_mqn = sum(1 for g in group_of.values() if g == "mqn")
    n_bg = sum(1 for g in group_of.values() if g == "bg")
    print(f"[load] anchor genomes: futalosine (mqn) {n_mqn}, background (bg) {n_bg}")

    tax = load_taxonomy()
    cls = load_classification()
    print(f"[load] taxonomy for {len(tax)} accessions, classifier for {len(cls)}")

    # protein -> anchor genomes it was pulled from (used for the phylum column)
    protein_anchors = defaultdict(set)
    for acc, rows in fam_rows.items():
        for pid, anchor, _g, _gi, _gap in rows:
            protein_anchors[pid].add(anchor)

    # ---- per-family sweep -------------------------------------------------
    sweep = []
    win_rows_of = {}
    win_genomes_of = {}
    for acc, rows in fam_rows.items():
        wrows = [r for r in rows if in_window(r[3], r[4], xz)]
        win_rows_of[acc] = wrows
        tgt = {a for _p, a, g, _gi, _gap in wrows if g == "mqn"}
        bgs = {a for _p, a, g, _gi, _gap in wrows if g == "bg"}
        win_genomes_of[acc] = (tgt, bgs)
        name = fam_name[acc]
        both, denom, rec = text_recall(name, fam_prot[acc], products)
        ph, ph_n, ph_tot = dominant_phylum(fam_prot[acc], protein_anchors, tax)
        same, tot, med = operon_stats(wrows, pair_rows)
        sweep.append({
            "family": name, "pfam": acc,
            "target_n": len(tgt), "bg_n": len(bgs),
            "recall_num": both, "recall_den": denom, "recall_pct": rec,
            "buried_pct": 100.0 - rec,
            "phylum": ph, "phylum_n": ph_n, "phylum_tot": ph_tot,
            "same": same, "strand_tot": tot, "median_gap": med,
        })
    sweep.sort(key=lambda d: (-d["buried_pct"], -d["target_n"]))

    # ---- printed sweep ----------------------------------------------------
    print("\n" + "=" * 78)
    print("MOST-BURIED DARK FAMILIES NEAR THE FUTALOSINE ENZYMES")
    print("=" * 78)
    print("  Restricted to DUF/Y families (the metric is only interpretable where")
    print("  the product name IS the family identifier), present in at least")
    print(f"  {MIN_GENOMES_FOR_SWEEP} target genomes, and absent from the background.")
    print(f"  {'family':<12}{'pfam':<11}{'target':<14}{'background':<13}"
          f"{'text recall':<18}{'buried'}")
    shown = 0
    for d in sweep:
        if shown >= args.top:
            break
        if d["target_n"] < MIN_GENOMES_FOR_SWEEP or d["bg_n"] != 0:
            continue
        if not DARK_NAME_RE.match(d["family"]):
            continue
        tgt_s = f"{d['target_n']}/{n_mqn}"
        bg_s = f"{d['bg_n']}/{n_bg}"
        rec_s = f"{d['recall_num']}/{d['recall_den']} = {d['recall_pct']:.1f}%"
        print(f"  {d['family']:<12}{d['pfam']:<11}{tgt_s:<14}{bg_s:<13}"
              f"{rec_s:<18}{d['buried_pct']:.1f}%")
        shown += 1

    # ---- the flagged five -------------------------------------------------
    print("\n" + "=" * 78)
    print("FLAGGED DARK FAMILIES (the published set)")
    print("=" * 78)
    by_acc = {d["pfam"]: d for d in sweep}
    flagged_rows = []
    for name, acc in FLAGGED:
        d = by_acc.get(acc)
        if d is None:
            print(f"  {name} ({acc}): NO HITS in {os.path.basename(DOM_MQN)}")
            continue
        ph_str = (f"{d['phylum']}({d['phylum_n']}/{d['phylum_tot']})"
                  if d["phylum_tot"] and
                  100.0 * d["phylum_n"] / d["phylum_tot"] >= PHYLUM_DOMINANCE
                  else "mixed")
        strand_str = (f"{d['same']}/{d['strand_tot']}="
                      f"{100.0 * d['same'] / d['strand_tot']:.1f}%"
                      if d["strand_tot"] else "n/a")
        med_str = f"{d['median_gap']:g}" if d["median_gap"] is not None else "n/a"
        flagged_rows.append([
            d["family"], d["pfam"], f"{d['target_n']}/{n_mqn}", f"{d['bg_n']}/{n_bg}",
            f"{d['recall_pct']:.1f}", f"{d['buried_pct']:.1f}",
            ph_str, strand_str, med_str,
        ])
        print(f"\n  {d['family']} ({d['pfam']})")
        print(f"    target genomes in window   : {fmt_frac(d['target_n'], n_mqn)}")
        print(f"    background genomes in window: {fmt_frac(d['bg_n'], n_bg)}")
        print(f"    text recall                : {d['recall_num']}/{d['recall_den']}"
              f" = {d['recall_pct']:.1f}%")
        print(f"    buried by annotation       : "
              f"{d['recall_den'] - d['recall_num']}/{d['recall_den']}"
              f" = {d['buried_pct']:.1f}%")
        print(f"    dominant phylum            : {ph_str}"
              f"  (top phylum {d['phylum']} {d['phylum_n']}/{d['phylum_tot']}, "
              f"'mixed' below {PHYLUM_DOMINANCE:g}%)")
        print(f"    same-strand with anchor    : {strand_str}")
        print(f"    median intergenic gap      : {med_str} bp")

    # ---- clade-matched controls ------------------------------------------
    ok_checks = []

    def check(label, observed, expected):
        ok_checks.append((label, observed, expected, observed == expected))

    def check_p(label, observed, expected, rel=0.05):
        good = expected > 0 and abs(observed - expected) / expected <= rel
        ok_checks.append((label, f"{observed:.3g}", f"{expected:.3g}", good))

    mqn_accs = [a for a, g in group_of.items() if g == "mqn"]
    bg_accs = [a for a, g in group_of.items() if g == "bg"]

    def is_actino(a):
        return bool(ACTINO_RE.search(tax.get(a, "") or ""))

    def is_firmicute(a):
        return (tax.get(a, "") or "") == FIRMICUTE

    duf_pos = set().union(*win_genomes_of["PF14012"])
    d1450_pos = set().union(*win_genomes_of["PF07293"])
    ytxc_pos = set().union(*win_genomes_of["PF08812"])

    act_mqnE = [a for a in mqn_accs if is_actino(a) and cls.get(a) == "MqnE"]
    act_mqnC = [a for a in mqn_accs if is_actino(a) and cls.get(a) == "MqnC"]
    act_bg = [a for a in bg_accs if is_actino(a)]
    fir_mqn = [a for a in mqn_accs if is_firmicute(a)]
    fir_mqnE = [a for a in fir_mqn if cls.get(a) == "MqnE"]
    fir_mqnC = [a for a in fir_mqn if cls.get(a) == "MqnC"]
    fir_bg = [a for a in bg_accs if is_firmicute(a)]

    # BioB and LipA: separate scan, PF14012 only is needed from it
    _gn, _gp, _gr = parse_domtbl(DOM_GEN, keep_accessions={"PF14012"})
    gen_win = [r for r in _gr.get("PF14012", []) if in_window(r[3], r[4], xz)]
    gen_pos = {a for _p, a, _g, _gi, _gap in gen_win}
    gen_group = {}
    with open(TSV_GEN, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gen_group[r["query_acc"]] = r["group"]
    biob = [a for a, g in gen_group.items() if g == "BioB" and is_actino(a)]
    lipa = [a for a, g in gen_group.items() if g == "LipA" and is_actino(a)]

    def npos(accs, positives):
        return sum(1 for a in accs if a in positives)

    print("\n" + "=" * 78)
    print("CLADE-MATCHED CONTROL 1: Actinomycetota, DUF4229 (PF14012)")
    print("=" * 78)
    # The BioB and LipA anchors live in the other scan, so their DUF4229 calls
    # come from gen_pos rather than duf_pos.
    rows_act = [("MqnE (futalosine)", act_mqnE, duf_pos),
                ("MqnC (futalosine)", act_mqnC, duf_pos),
                ("background radical SAM", act_bg, duf_pos),
                ("BioB biotin synthase", biob, gen_pos),
                ("LipA lipoyl synthase", lipa, gen_pos)]
    for label, accs, positives in rows_act:
        print(f"  {label:<26}{fmt_frac(npos(accs, positives), len(accs))}")
    pooled = act_bg + biob + lipa
    pooled_pos = npos(act_bg, duf_pos) + npos(biob, gen_pos) + npos(lipa, gen_pos)
    print(f"  {'pooled non-MqnE':<26}{fmt_frac(pooled_pos, len(pooled))}"
          f"   (background + BioB + LipA)")
    p_pool, or_pool = fisher(npos(act_mqnE, duf_pos), len(act_mqnE), pooled_pos, len(pooled))
    p_bg, or_bg = fisher(npos(act_mqnE, duf_pos), len(act_mqnE), npos(act_bg, duf_pos), len(act_bg))
    print(f"  Fisher MqnE vs background : p = {p_bg:.2g}, "
          f"odds ratio (Haldane-Anscombe) = {or_bg:.1f}x")
    print(f"  Fisher MqnE vs pooled     : p = {p_pool:.2g}, "
          f"odds ratio (Haldane-Anscombe) = {or_pool:.1f}x")

    print("\n" + "=" * 78)
    print("CLADE-MATCHED CONTROL 2: Bacillota, enzyme-restricted")
    print("=" * 78)
    print(f"  denominators: background {len(fir_bg)}, futalosine {len(fir_mqn)}, "
          f"MqnE {len(fir_mqnE)}, MqnC {len(fir_mqnC)}")
    p_1450, or_1450 = fisher(npos(fir_mqnE, d1450_pos), len(fir_mqnE),
                             npos(fir_bg, d1450_pos), len(fir_bg))
    p_ytxc, or_ytxc = fisher(npos(fir_mqnC, ytxc_pos), len(fir_mqnC),
                             npos(fir_bg, ytxc_pos), len(fir_bg))
    print(f"  DUF1450 in MqnE  {fmt_frac(npos(fir_mqnE, d1450_pos), len(fir_mqnE))}"
          f"   vs background {fmt_frac(npos(fir_bg, d1450_pos), len(fir_bg))}"
          f"   p = {p_1450:.2g}, OR = {or_1450:.1f}x")
    print(f"  YtxC    in MqnC  {fmt_frac(npos(fir_mqnC, ytxc_pos), len(fir_mqnC))}"
          f"   vs background {fmt_frac(npos(fir_bg, ytxc_pos), len(fir_bg))}"
          f"   p = {p_ytxc:.2g}, OR = {or_ytxc:.1f}x")
    print("  cross-checks (the partners are enzyme-specific, not clade-specific only):")
    print(f"    DUF1450 in MqnC  {fmt_frac(npos(fir_mqnC, d1450_pos), len(fir_mqnC))}")
    print(f"    YtxC    in MqnE  {fmt_frac(npos(fir_mqnE, ytxc_pos), len(fir_mqnE))}")
    print(f"    DUF4229 in all Bacillota futalosine  "
          f"{fmt_frac(npos(fir_mqn, duf_pos), len(fir_mqn))}")

    # ---- PASS / FAIL ------------------------------------------------------
    duf = by_acc["PF14012"]
    check("DUF4229 near futalosine", f"{duf['target_n']}/{n_mqn}", f"342/3648")
    check("DUF4229 near background", f"{duf['bg_n']}/{n_bg}", "0/1211")
    check("DUF4229 text recall", f"{duf['recall_pct']:.1f}%", "56.9%")
    check("actinobacterial MqnE with DUF4229",
          f"{npos(act_mqnE, duf_pos)}/{len(act_mqnE)}", "280/329")
    check("actinobacterial MqnC with DUF4229",
          f"{npos(act_mqnC, duf_pos)}/{len(act_mqnC)}", "55/566")
    check("actinobacterial background with DUF4229",
          f"{npos(act_bg, duf_pos)}/{len(act_bg)}", "0/212")
    check("actinobacterial BioB with DUF4229",
          f"{npos(biob, gen_pos)}/{len(biob)}", "0/259")
    check("actinobacterial LipA with DUF4229",
          f"{npos(lipa, gen_pos)}/{len(lipa)}", "0/300")
    check("pooled non-MqnE with DUF4229", f"{pooled_pos}/{len(pooled)}", "0/771")
    check_p("Fisher, actino MqnE vs pooled non-MqnE", p_pool, 3.0e-211)
    check("Bacillota background genomes", str(len(fir_bg)), "234")
    check("Bacillota futalosine genomes", str(len(fir_mqn)), "890")
    check("Bacillota MqnE genomes", str(len(fir_mqnE)), "468")
    check("Bacillota MqnC genomes", str(len(fir_mqnC)), "422")
    check("DUF1450 in firmicute MqnE",
          f"{npos(fir_mqnE, d1450_pos)}/{len(fir_mqnE)}", "298/468")
    check("DUF1450 in firmicute background",
          f"{npos(fir_bg, d1450_pos)}/{len(fir_bg)}", "0/234")
    check_p("Fisher, DUF1450 MqnE vs firmicute background", p_1450, 2.8e-75)
    check("YtxC in firmicute MqnC",
          f"{npos(fir_mqnC, ytxc_pos)}/{len(fir_mqnC)}", "240/422")
    check("YtxC in firmicute background",
          f"{npos(fir_bg, ytxc_pos)}/{len(fir_bg)}", "0/234")
    check_p("Fisher, YtxC MqnC vs firmicute background", p_ytxc, 1.9e-62)

    print("\n" + "=" * 78)
    print("PASS / FAIL AGAINST THE CANONICAL NUMBERS")
    print("=" * 78)
    print(f"  {'check':<46}{'computed':<16}{'canonical':<16}{'verdict'}")
    n_fail = 0
    for label, obs, exp, good in ok_checks:
        if not good:
            n_fail += 1
        print(f"  {label:<46}{obs:<16}{exp:<16}{'PASS' if good else 'FAIL'}")
    print(f"\n  {len(ok_checks) - n_fail}/{len(ok_checks)} = "
          f"{100.0 * (len(ok_checks) - n_fail) / len(ok_checks):.1f}% of canonical "
          f"numbers reproduced")
    if n_fail:
        print(f"  {n_fail}/{len(ok_checks)} FAILED. Do not use this output.")
        if xz:
            print("  (--exclude-index-zero was given; see the window note in the source.)")
    else:
        print("  ALL PASS.")

    # ---- output -----------------------------------------------------------
    if not args.write:
        print("\n[dry run] nothing written. Re-run with --write to update:")
        print(f"          {args.out_flagged}")
        print(f"          {args.out_sweep}")
        return 1 if n_fail else 0

    header = ("# Sequence-vs-label re-audit across ALL Pfam families "
              "(Philipp's \"any protein buried by bad annotation?\")\n"
              f"# Regenerated by scripts/text_recall_sweep.py. Detection: HMMER --cut_ga on the "
              f"{len(products)} distinct neighbour proteins\n"
              f"# of {os.path.basename(TSV_MQN)}, scan {os.path.basename(DOM_MQN)}. "
              f"Metric validated: DUF4229 text-recall = "
              f"{duf['recall_num']}/{duf['recall_den']} = {duf['recall_pct']:.1f}%.\n"
              "# text_recall = fraction of sequence-detected members a product-name search "
              "for the family name recovers.\n"
              f"# Window: gap <= {MAX_GAP} bp and abs(gene_index) <= {MAX_GENE_INDEX}; "
              f"gene_index == 0 rows {'dropped' if xz else 'kept'} (overlapping neighbours, "
              "not the anchor CDS).\n"
              "# same_strand_operon counts (anchor, neighbour protein) pairs expanded over "
              "every neighbour-table row for the pair.\n")
    with open(args.out_flagged, "w", newline="") as fh:
        fh.write(header)
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "pfam", "target_all", "bg_all", "text_recall_pct",
                    "buried_pct", "dominant_phylum", "same_strand_operon",
                    "median_gap_bp"])
        for row in flagged_rows:
            w.writerow(row)
        fh.write("#\n# CLADE-MATCHED CONTROLS (taxonomy held constant)\n")
        fh.write(f"# Actinomycetota: DUF4229 MqnE "
                 f"{fmt_frac(npos(act_mqnE, duf_pos), len(act_mqnE))} vs background "
                 f"{fmt_frac(npos(act_bg, duf_pos), len(act_bg))} (p={p_bg:.1e}); "
                 f"MqnC {fmt_frac(npos(act_mqnC, duf_pos), len(act_mqnC))}\n")
        fh.write(f"# Actinomycetota pooled non-MqnE (background + BioB + LipA): "
                 f"{fmt_frac(pooled_pos, len(pooled))}, p={p_pool:.1e}\n")
        fh.write(f"# Bacillota, ENZYME-RESTRICTED (two-way mqn_seq_classification): "
                 f"DUF1450 -> MqnE {fmt_frac(npos(fir_mqnE, d1450_pos), len(fir_mqnE))} "
                 f"vs {fmt_frac(npos(fir_bg, d1450_pos), len(fir_bg))} (p={p_1450:.1e}); "
                 f"YtxC -> MqnC {fmt_frac(npos(fir_mqnC, ytxc_pos), len(fir_mqnC))} "
                 f"vs {fmt_frac(npos(fir_bg, ytxc_pos), len(fir_bg))} (p={p_ytxc:.1e})\n")
        fh.write(f"# DUF4229 is {fmt_frac(npos(fir_mqn, duf_pos), len(fir_mqn))} in "
                 f"Bacillota futalosine genomes; YtxC and DUF1450 are absent from "
                 f"Actinomycetota. Mutually exclusive by clade.\n")
    print(f"\n[out] {args.out_flagged}")

    with open(args.out_sweep, "w", newline="") as fh:
        fh.write("# Full text-recall sweep over every Pfam family in "
                 f"{os.path.basename(DOM_MQN)}. scripts/text_recall_sweep.py\n")
        fh.write("# text_recall_num/den are protein counts with NO window filter; "
                 "target/bg are anchor-genome counts WITHIN the window.\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "pfam", "target_pos", "target_n", "bg_pos", "bg_n",
                    "text_recall_num", "text_recall_den", "text_recall_pct",
                    "buried_pct", "dominant_phylum", "dominant_phylum_n",
                    "dominant_phylum_den", "same_strand", "same_strand_den",
                    "median_gap_bp"])
        for d in sweep:
            w.writerow([d["family"], d["pfam"], d["target_n"], n_mqn, d["bg_n"], n_bg,
                        d["recall_num"], d["recall_den"], f"{d['recall_pct']:.1f}",
                        f"{d['buried_pct']:.1f}", d["phylum"], d["phylum_n"],
                        d["phylum_tot"], d["same"], d["strand_tot"],
                        "" if d["median_gap"] is None else f"{d['median_gap']:g}"])
    print(f"[out] {args.out_sweep}  ({len(sweep)} families)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
