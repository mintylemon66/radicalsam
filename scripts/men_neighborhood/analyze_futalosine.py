#!/usr/bin/env python3
"""Count mqnP, mqnE and the classical men genes near DUF4229 in genomes WITH MqnE and test the enrichment.

    python3 _scripts/men_neighborhood/analyze_futalosine.py      # after call_men_genes.py --futalosine

The companion of analyze_windows.py for the 931 kept reference proteomes that
carry both DUF4229 and MqnE (NCBIfam TIGR03700), whose files the earlier
scripts write to FUTALOSINE_DIR when given --futalosine. The question is
whether DUF4229 sits near mqnP (NCBIfam NF041584), the UbiA-family membrane
prenyltransferase of the futalosine branch, as it sits near menA, the
prenyltransferase of the classical branch, in the genomes without MqnE.

Every count and test is the one in analyze_windows.py, whose functions this
script imports: DUF4229 windows (10 genes either side) that hold the gene
against control windows that hold it, the one-sided Fisher exact test, the
exact matched-set test, and the Fisher test again with one DUF4229 locus per
genus. mqnE is counted the same way, as a check on the paper's operon figures.
mqnP and mqnE belong to the same branch, so the loci are also split by whether
mqnE is in the window; an mqnP signal that only comes along with mqnE shows up
there. Genome-level presence is read from the main run's genome_men.tsv, which
covers every kept proteome.
"""

import argparse
import math
import random
import sys
from collections import Counter, defaultdict

from analyze_windows import (fisher_log10_p, fmt_p, genes_in, load_windows, matched_log10_p, orientation, pct,
                             proteomes_in, read_tsv)

OUT_DIR = "_data/men_neighborhood/"
FUTALOSINE_DIR = OUT_DIR + "futalosine/"
TARGETS_TSV = FUTALOSINE_DIR + "targets.tsv"
WINDOWS_TSV = FUTALOSINE_DIR + "windows.tsv"
CALLS_TSV = FUTALOSINE_DIR + "men_calls.tsv"
GENOME_TSV = OUT_DIR + "genome_men.tsv"           # written by fetch_genome_men.py for all kept proteomes
SUMMARY_TSV = FUTALOSINE_DIR + "summary_by_gene.tsv"
SUBSETS_TSV = FUTALOSINE_DIR + "summary_by_subset.tsv"
ORDERS_TSV = FUTALOSINE_DIR + "summary_by_order.tsv"
HITS_TSV = FUTALOSINE_DIR + "men_hits_in_windows.tsv"
LOCI_TSV = FUTALOSINE_DIR + "duf4229_loci.tsv"
RESULTS_MD = FUTALOSINE_DIR + "RESULTS_TABLES.md"

GENOMEWIDE_DIR = "_data/evoweaver/genomewide/"
KEPT_TXT = GENOMEWIDE_DIR + "proteomes_kept.txt"

PRENYL = "MqnP_or_MenA"        # either menaquinone prenyltransferase; the closer copy counts
OTHER_UBIA = "UbiA_other"      # a Pfam UbiA (PF01040) protein called neither MqnP nor MenA
LABELS = ["MqnP", "MenA", PRENYL, "UbiA", OTHER_UBIA, "MqnE", "MqnC",
          "MenE", "MenB", "MenD", "MenC", "MenF", "MenH", "MenG", "any_classical", "MenJ"]
SUBSET_LABELS = ["MqnP", "MenA", PRENYL, "UbiA", OTHER_UBIA, "MqnE"]
ARRANGEMENT_GENES = ["MqnE", "MqnP", "MqnC", "MenA"]
ANCHORS = [("UP000001973", "Streptomyces coelicolor A3(2)")]
SEED = 4229
CONTROLS_PER_LOCUS = 5                      # as in build_windows.py


def log(msg=""):
    print(msg, flush=True)


def labels_in(window):
    """genes_in() of analyze_windows.py, plus PRENYL and OTHER_UBIA."""
    found = genes_in(window)
    near = [found[g] for g in ("MqnP", "MenA") if g in found]
    if near:
        found[PRENYL] = min(near, key=lambda m: (abs(m["offset"]), m["offset"]))
    for m in sorted(window["members"], key=lambda m: (abs(m["offset"]), m["offset"])):
        if m["offset"] != 0 and m["ubia"] and m["call"] not in ("MqnP", "MenA"):
            found[OTHER_UBIA] = m
            break
    return found


def arrangement(window):
    """The DUF4229 locus as a row of genes, e.g. 'mqnE> . DUF4229> . mqnP<'.

    Only ARRANGEMENT_GENES are named and '.' stands for any other gene; the row
    runs from the farthest named gene upstream to the farthest downstream, and
    each arrow gives the strand relative to DUF4229.
    """
    named = {m["offset"]: m for m in window["members"] if m["call"] in ARRANGEMENT_GENES}
    if not named:
        return "DUF4229> (no mqn or menA gene within 10)"
    lo, hi = min(min(named), 0), max(max(named), 0)
    cells = []
    for offset in range(lo, hi + 1):
        if offset == 0:
            cells.append("DUF4229>")
        elif offset in named:
            m = named[offset]
            gene = m["call"][0].lower() + m["call"][1:]
            cells.append(gene + (">" if m["rel_strand"] == "same" else "<"))
        else:
            cells.append(".")
    return " ".join(cells)


def describe(hits):
    """Offset, strand and adjacency of the nearest copy at each locus that has one."""
    offsets = Counter(m["offset"] for m in hits)
    same = sum(1 for m in hits if m["rel_strand"] == "same")
    adjacent = Counter(orientation(m["offset"], m["rel_strand"]) for m in hits if abs(m["offset"]) == 1)
    return {"top_offsets": ", ".join(f"{o:+d}: {n}" for o, n in offsets.most_common(4)),
            "same_strand": same, "n_hits": len(hits), "adjacent": dict(adjacent),
            "within_2": sum(1 for m in hits if abs(m["offset"]) <= 2),
            "downstream": sum(1 for m in hits if m["offset"] > 0)}


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    targets = read_tsv(TARGETS_TSV)
    lineage = {t["proteome"]: t for t in targets}
    calls = {r["protein_id"]: r for r in read_tsv(CALLS_TSV)}
    windows = load_windows(WINDOWS_TSV, calls, {})
    genome = {r["proteome"]: r for r in read_tsv(GENOME_TSV)}

    duf = {w["locus"]: w for w in windows.values() if w["kind"] == "duf"}
    controls = defaultdict(list)
    for w in windows.values():
        if w["kind"] == "control":
            controls[w["locus"]].append(w)
    found = {wid: labels_in(w) for wid, w in windows.items()}
    duf_found = {locus: found[w["id"]] for locus, w in duf.items()}
    n_duf = len(duf)
    n_ctrl = sum(len(v) for v in controls.values())
    target_genomes = {t["proteome"] for t in targets}
    log(f"{len(targets):,} DUF4229 proteins in {len(target_genomes):,} genomes with MqnE; "
        f"{n_duf:,} placed on their contigs ({len({w['proteome'] for w in duf.values()}):,} genomes); "
        f"{n_ctrl:,} control windows")
    full = 2 * max(abs(m["offset"]) for w in duf.values() for m in w["members"]) + 1
    for kind in ("duf", "control"):
        group = [w for w in windows.values() if w["kind"] == kind]
        short = sum(1 for w in group if len(w["members"]) < full)
        log(f"  {kind} windows cut short by a contig end: {pct(short, len(group))}")
    log(f"  loci with fewer than {CONTROLS_PER_LOCUS} control windows: "
        f"{sum(1 for l in duf if len(controls[l]) < CONTROLS_PER_LOCUS):,}")

    # one locus per genus, for the dereplicated test (same draw as analyze_windows.py)
    rng = random.Random(args.seed)
    by_genus = defaultdict(list)
    for locus, w in sorted(duf.items()):
        by_genus[lineage[w["proteome"]]["genus"]].append(locus)
    genus_loci = [rng.choice(loci) for genus, loci in sorted(by_genus.items())]

    def count(label, loci):
        a = sum(1 for locus in loci if label in duf_found[locus])
        ctrl = [w for locus in loci for w in controls[locus]]
        c = sum(1 for w in ctrl if label in found[w["id"]])
        return a, len(loci), c, len(ctrl)

    def matched_sets(label, loci):
        sets = []
        for locus in loci:
            d = int(label in duf_found[locus])
            t = d + sum(1 for w in controls[locus] if label in found[w["id"]])
            sets.append((d, t, len(controls[locus])))
        return sets

    def test(label, loci, genus):
        a, n1, c, n2 = count(label, loci)
        rate1, rate2 = (a / n1 if n1 else 0), (c / n2 if n2 else 0)
        m_log_p, observed, expected = matched_log10_p(matched_sets(label, loci))
        ga, gn1, gc, gn2 = count(label, genus)
        row = {"gene": label, "duf_windows": a, "duf_total": n1, "control_windows": c, "control_total": n2,
               "duf_rate": rate1, "control_rate": rate2, "enrichment": rate1 / rate2 if rate2 else math.inf,
               "fisher_log10_p": fisher_log10_p(a, n1, c, n2), "matched_log10_p": m_log_p,
               "matched_expected": expected, "genus_duf": ga, "genus_total": gn1, "genus_control": gc,
               "genus_control_total": gn2, "genus_fisher_log10_p": fisher_log10_p(ga, gn1, gc, gn2)}
        row.update(describe([duf_found[locus][label] for locus in loci if label in duf_found[locus]]))
        return row

    # ---- per gene
    all_loci = sorted(duf)
    summary = [test(label, all_loci, genus_loci) for label in LABELS]
    with open(SUMMARY_TSV, "w") as out:
        out.write("gene\tduf_windows_with_gene\tduf_windows\tcontrol_windows_with_gene\tcontrol_windows\t"
                  "duf_rate\tcontrol_rate\tenrichment\tfisher_log10_p\tmatched_log10_p\tmatched_expected\t"
                  "genus_duf_with_gene\tgenera\tgenus_control_with_gene\tgenus_control_windows\tgenus_fisher_log10_p\t"
                  "top_offsets\tsame_strand\tadjacent_orientation\twithin_2_genes\tdownstream\n")
        for s in summary:
            out.write(f"{s['gene']}\t{s['duf_windows']}\t{s['duf_total']}\t{s['control_windows']}\t{s['control_total']}\t"
                      f"{s['duf_rate']:.4f}\t{s['control_rate']:.5f}\t{s['enrichment']:.1f}\t{s['fisher_log10_p']:.1f}\t"
                      f"{s['matched_log10_p']:.1f}\t{s['matched_expected']:.1f}\t{s['genus_duf']}\t{s['genus_total']}\t"
                      f"{s['genus_control']}\t{s['genus_control_total']}\t{s['genus_fisher_log10_p']:.1f}\t"
                      f"{s['top_offsets']}\t{s['same_strand']}\t"
                      f"{';'.join(f'{k}={v}' for k, v in sorted(s['adjacent'].items()))}\t{s['within_2']}\t"
                      f"{s['downstream']}\n")

    log()
    log(f"{n_duf:,} DUF4229 loci, {n_ctrl:,} control windows, {len(genus_loci)} genera (one locus each); "
        f"log10 p in summary_by_gene.tsv")
    log("gene          DUF4229 windows        control windows       enrichment  Fisher p   matched p  "
        "genera (DUF4229 vs control)        top offsets")
    for s in summary:
        log(f"{s['gene']:<13} {pct(s['duf_windows'], s['duf_total']):<22} {pct(s['control_windows'], s['control_total']):<21} "
            f"{s['enrichment']:>8.1f}x  {fmt_p(s['fisher_log10_p']):>9}  {fmt_p(s['matched_log10_p']):>9}  "
            f"{pct(s['genus_duf'], s['genus_total'])} vs {pct(s['genus_control'], s['genus_control_total'])}   "
            f"{s['top_offsets']}")

    # ---- position and strand of the two prenyltransferases and of mqnE
    log()
    log("position of the nearest copy (offsets in DUF4229's direction of transcription)")
    by_label = {s["gene"]: s for s in summary}
    for label in ("MqnP", "MenA", PRENYL, OTHER_UBIA, "MqnE"):
        s = by_label[label]
        n = s["n_hits"]
        if not n:
            log(f"  {label}: none within 10 genes")
            continue
        log(f"  {label}: {n:,} loci; same strand {pct(s['same_strand'], n)}; downstream (offset > 0) "
            f"{pct(s['downstream'], n)}; within 2 genes {pct(s['within_2'], n)}; adjacent "
            f"{sum(s['adjacent'].values())} ({', '.join(f'{k} {v}' for k, v in sorted(s['adjacent'].items()))}); "
            f"top offsets {s['top_offsets']}")
    mqne_hits = [duf_found[l]["MqnE"] for l in all_loci if "MqnE" in duf_found[l]]
    if mqne_hits:
        at_minus2 = sum(1 for m in mqne_hits if m["offset"] == -2)
        upstream = sum(1 for m in mqne_hits if m["offset"] < 0)
        log(f"  check against the paper (DUF4229 two genes downstream of mqnE, same strand 341 of 342): "
            f"mqnE at -2 in {pct(at_minus2, len(mqne_hits))}, upstream of DUF4229 (DUF4229 downstream) in "
            f"{pct(upstream, len(mqne_hits))}, same strand {pct(sum(m['rel_strand'] == 'same' for m in mqne_hits), len(mqne_hits))}")

    # ---- does mqnP come with mqnE?
    log()
    log("mqnP and mqnE together (DUF4229 loci; control windows in brackets)")
    table = Counter((("MqnP" in duf_found[l]), ("MqnE" in duf_found[l])) for l in all_loci)
    ctrl_table = Counter((("MqnP" in found[w["id"]]), ("MqnE" in found[w["id"]]))
                         for l in all_loci for w in controls[l])
    for p in (True, False):
        for e in (True, False):
            log(f"  mqnP {'yes' if p else 'no ':<3}  mqnE {'yes' if e else 'no ':<3}  {table[(p, e)]:>5,}  "
                f"[{ctrl_table[(p, e)]:,}]")
    both = [l for l in all_loci if "MqnP" in duf_found[l] and "MqnE" in duf_found[l]]
    if both:
        same_side = sum(1 for l in both if (duf_found[l]["MqnP"]["offset"] > 0) == (duf_found[l]["MqnE"]["offset"] > 0))
        closer = sum(1 for l in both if abs(duf_found[l]["MqnP"]["offset"]) < abs(duf_found[l]["MqnE"]["offset"]))
        gap = Counter(duf_found[l]["MqnP"]["offset"] - duf_found[l]["MqnE"]["offset"] for l in both)
        log(f"  loci with both: {len(both):,}; mqnP on the same side of DUF4229 as mqnE {pct(same_side, len(both))}; "
            f"mqnP closer to DUF4229 than mqnE {pct(closer, len(both))}; "
            f"mqnP offset minus mqnE offset, commonest: {', '.join(f'{d:+d}: {n}' for d, n in gap.most_common(5))}")

    # ---- the same tests in subsets of loci
    has_mqnp = {up for up in target_genomes if int(genome[up]["n_MqnP"]) > 0}
    has_mena = {up for up in target_genomes if int(genome[up]["n_MenA"]) > 0}
    subsets = [
        ("all loci", all_loci),
        ("genomes with MqnP (NF041584, UniProtKB)", [l for l in all_loci if duf[l]["proteome"] in has_mqnp]),
        ("genomes without MqnP", [l for l in all_loci if duf[l]["proteome"] not in has_mqnp]),
        ("mqnE within 10 genes of DUF4229", [l for l in all_loci if "MqnE" in duf_found[l]]),
        ("no mqnE within 10 genes", [l for l in all_loci if "MqnE" not in duf_found[l]]),
    ]
    log()
    log("subsets (DUF4229 windows against their own control windows)")
    with open(SUBSETS_TSV, "w") as out:
        out.write("subset\tloci\tgenomes\tgene\tduf_windows_with_gene\tduf_windows\tcontrol_windows_with_gene\t"
                  "control_windows\tfisher_log10_p\tmatched_log10_p\ttop_offsets\tsame_strand\n")
        for name, loci in subsets:
            ups = {duf[l]["proteome"] for l in loci}
            chosen = set(loci)
            sub_genus = [l for l in genus_loci if l in chosen]
            log(f"  {name}: {len(ups):,} genomes, {len(loci):,} loci")
            for label in SUBSET_LABELS:
                s = test(label, loci, sub_genus)
                out.write(f"{name}\t{len(loci)}\t{len(ups)}\t{label}\t{s['duf_windows']}\t{s['duf_total']}\t"
                          f"{s['control_windows']}\t{s['control_total']}\t{s['fisher_log10_p']:.1f}\t"
                          f"{s['matched_log10_p']:.1f}\t{s['top_offsets']}\t{s['same_strand']}\n")
                log(f"    {label:<13} {pct(s['duf_windows'], s['duf_total']):<22} control "
                    f"{pct(s['control_windows'], s['control_total']):<20} Fisher p {fmt_p(s['fisher_log10_p']):>9}  "
                    f"matched p {fmt_p(s['matched_log10_p']):>9}")

    # ---- genome level
    kept = set(open(KEPT_TXT).read().split())
    with_mqnc = proteomes_in("MqnC") & kept
    with_pf02621 = proteomes_in("MqnAD") & kept
    log()
    log("genome level (UniProtKB InterPro matches, genome_men.tsv)")
    log(f"  MqnP NF041584 anywhere in the genome: {pct(len(has_mqnp), len(target_genomes))}")
    log(f"  MenA TIGR00751 anywhere in the genome: {pct(len(has_mena), len(target_genomes))}")
    log(f"  MenD and MenB: {pct(sum(1 for up in target_genomes if genome[up]['classical'] == '1'), len(target_genomes))}")
    log(f"  MqnC TIGR03699: {pct(len(target_genomes & with_mqnc), len(target_genomes))}; "
        f"PF02621 (MqnA/MqnD): {pct(len(target_genomes & with_pf02621), len(target_genomes))}")
    log(f"  phyla: {dict(Counter(lineage[up]['phylum'] for up in target_genomes).most_common())}")

    # ---- per order
    by_order = defaultdict(list)
    for locus, w in duf.items():
        by_order[lineage[w["proteome"]]["order"] or "(none)"].append(locus)
    log()
    log("order                  loci  genera  mqnP in window   mqnE in window   menA in window   mqnP or menA     "
        "control mqnP or menA   genomes with MqnP")
    with open(ORDERS_TSV, "w") as out:
        out.write("order\tloci\tgenera\tMqnP_in_window\tMqnE_in_window\tMenA_in_window\tMqnP_or_MenA_in_window\t"
                  "control_MqnP_or_MenA\tcontrol_windows\tgenomes_with_MqnP\tgenomes\tcommonest_arrangement\t"
                  "arrangement_loci\n")
        for order, loci in sorted(by_order.items(), key=lambda x: -len(x[1])):
            genera = len({lineage[duf[l]["proteome"]]["genus"] for l in loci})
            n = {g: sum(1 for l in loci if g in duf_found[l]) for g in ("MqnP", "MqnE", "MenA", PRENYL)}
            ctrl = [w for l in loci for w in controls[l]]
            cp = sum(1 for w in ctrl if PRENYL in found[w["id"]])
            ups = {duf[l]["proteome"] for l in loci}
            top, top_n = Counter(arrangement(duf[l]) for l in loci).most_common(1)[0]
            out.write(f"{order}\t{len(loci)}\t{genera}\t{n['MqnP']}\t{n['MqnE']}\t{n['MenA']}\t{n[PRENYL]}\t{cp}\t"
                      f"{len(ctrl)}\t{len(ups & has_mqnp)}\t{len(ups)}\t{top}\t{top_n}\n")
            log(f"{order:<22} {len(loci):>5}  {genera:>6}  {pct(n['MqnP'], len(loci)):<15}  {pct(n['MqnE'], len(loci)):<15}  "
                f"{pct(n['MenA'], len(loci)):<15}  {pct(n[PRENYL], len(loci)):<15}  {pct(cp, len(ctrl)):<21}  "
                f"{pct(len(ups & has_mqnp), len(ups))}")

    log()
    log("commonest arrangements (named: mqnE, mqnP, mqnC, menA; '.' any other gene; arrows relative to DUF4229)")
    for row, k in Counter(arrangement(duf[l]) for l in all_loci).most_common(10):
        log(f"  {k:>4}  {row}")

    # ---- every call in every window, and one row per locus
    with open(HITS_TSV, "w") as out:
        out.write("window_id\tkind\tproteome\tgenus\torder\tgene\toffset\trel_strand\torientation\tprotein_id\t"
                  "locus_tags\tgene_name\tproduct\n")
        for wid, w in sorted(windows.items()):
            for m in w["members"]:
                if m["offset"] == 0 or not (m["call"] or m["ubia"]):
                    continue
                t = lineage[w["proteome"]]
                for g in ([m["call"]] if m["call"] else []) + (["UbiA"] if m["ubia"] else []):
                    out.write(f"{wid}\t{w['kind']}\t{w['proteome']}\t{t['genus']}\t{t['order']}\t{g}\t{m['offset']}\t"
                              f"{m['rel_strand']}\t{orientation(m['offset'], m['rel_strand'])}\t{m['protein_id']}\t"
                              f"{m['locus_tags']}\t{m['gene_name']}\t{m['product']}\n")
    with open(LOCI_TSV, "w") as out:
        out.write("proteome\tduf_acc\tgenus\torder\tspecies\tcontig\tgenome_has_MqnP\t" +
                  "\t".join(f"{g}_offset" for g in LABELS) + "\tarrangement\tcontrol_windows\n")
        for locus in all_loci:
            w = duf[locus]
            t = lineage[w["proteome"]]
            f = duf_found[locus]
            cells = [f"{f[g]['offset']:+d}{'s' if f[g]['rel_strand'] == 'same' else 'o'}" if g in f else ""
                     for g in LABELS]
            out.write(f"{w['proteome']}\t{w['duf_acc']}\t{t['genus']}\t{t['order']}\t{t['species']}\t{w['contig']}\t"
                      f"{int(w['proteome'] in has_mqnp)}\t" + "\t".join(cells) +
                      f"\t{arrangement(w)}\t{len(controls[locus])}\n")

    # ---- product names behind the HMM calls, and an anchor locus
    log()
    log("GenBank product names of the nearest MqnP, MenA and MqnE calls in DUF4229 windows")
    for label in ("MqnP", "MenA", "MqnE"):
        names = Counter(duf_found[l][label]["product"] for l in all_loci if label in duf_found[l])
        log(f"  {label}: " + "; ".join(f"{p or '(none)'} {n}" for p, n in names.most_common(4)))
    for label, other in (("MqnP", "MenA"), ("MenA", "MqnP")):
        nearest = [duf_found[l][label]["protein_id"] for l in all_loci if label in duf_found[l]]
        also = sum(1 for pid in nearest if f"{other}:" in calls[pid]["all_men_hits"])
        log(f"  nearest {label} calls that also pass a {other} model: {pct(also, len(nearest))}")
    log()
    for up, name in ANCHORS:
        for locus in [l for l in all_loci if duf[l]["proteome"] == up]:
            w = duf[locus]
            log(f"anchor: {name} {up} {w['duf_acc']} on {w['contig']} (DUF4229 strand {w['center_strand']})")
            for m in sorted(w["members"], key=lambda m: m["offset"]):
                if abs(m["offset"]) <= 8:
                    log(f"  {m['offset']:+3d} {m['rel_strand']:<8} {m['locus_tags'] or m['gene_name']:<22} "
                        f"{m['call'] or '-':<5} {'UbiA' if m['ubia'] else '':<5} {m['product']}")
        if up not in {duf[l]["proteome"] for l in all_loci}:
            log(f"anchor {name} {up} not among the placed loci")

    write_results_md(summary, n_duf, n_ctrl, len(genus_loci), len(target_genomes), len(has_mqnp))
    return 0


def write_results_md(summary, n_duf, n_ctrl, n_genera, n_genomes, n_mqnp):
    with open(RESULTS_MD, "w") as out:
        out.write("# Menaquinone genes within 10 genes of DUF4229, genomes with MqnE\n\n")
        out.write(f"{n_genomes:,} genomes ({n_mqnp:,} with MqnP anywhere in the genome), {n_duf:,} DUF4229 loci, "
                  f"{n_ctrl:,} control windows, {n_genera} genera (one locus each).\n\n")
        out.write("| gene | DUF4229 windows | control windows | enrichment | Fisher p | matched p | genera, DUF4229 | "
                  "genera, control | top offsets | same strand |\n")
        out.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for s in summary:
            out.write(f"| {s['gene']} | {pct(s['duf_windows'], s['duf_total'])} | {pct(s['control_windows'], s['control_total'])} | "
                      f"{s['enrichment']:.1f} | {fmt_p(s['fisher_log10_p'])} | {fmt_p(s['matched_log10_p'])} | "
                      f"{pct(s['genus_duf'], s['genus_total'])} | {pct(s['genus_control'], s['genus_control_total'])} | "
                      f"{s['top_offsets']} | {s['same_strand']}/{s['n_hits']} |\n")


if __name__ == "__main__":
    sys.exit(main())
