#!/usr/bin/env python3
"""Count classical men genes near DUF4229 in genomes without MqnE and test the enrichment.

    python3 _scripts/men_neighborhood/analyze_windows.py         # after call_men_genes.py and fetch_genome_men.py

For each men gene the script counts the DUF4229 windows (WINDOW genes either
side) that hold at least one copy, the control windows that hold one, and
tests the difference two ways:

  pooled     one-sided Fisher exact test on the 2 x 2 table of windows
             (DUF4229 or control, with or without the gene);
  matched    an exact permutation test inside each locus's matched set
             (the DUF4229 window plus its control windows, all from the same
             contig): if DUF4229 sat at a random length-matched position, the
             DUF4229 window would hold the gene with probability t / (k + 1),
             where t of the k + 1 windows in the set hold it. The observed
             count is compared with the exact (Poisson binomial) distribution.

Related genomes are not independent, so the pooled test is repeated with one
DUF4229 locus per genus. The script also reports the offset and strand of each
gene relative to DUF4229, which genomes carry MenD and MenB anywhere in the
genome (genome_men.tsv), a per-order breakdown, the Corynebacterium
glutamicum ATCC 13032 anchor locus, and how the HMM calls compare with the
product names in the GenBank annotation. p-values are reported as log10 so
that very small values do not underflow.
"""

import argparse
import math
import random
import re
import sys
from collections import Counter, defaultdict

OUT_DIR = "_data/men_neighborhood/"
TARGETS_TSV = OUT_DIR + "targets.tsv"
WINDOWS_TSV = OUT_DIR + "windows.tsv"
CALLS_TSV = OUT_DIR + "men_calls.tsv"
GENOME_TSV = OUT_DIR + "genome_men.tsv"
SUMMARY_TSV = OUT_DIR + "summary_by_gene.tsv"
ORDERS_TSV = OUT_DIR + "summary_by_order.tsv"
CLASSES_TSV = OUT_DIR + "summary_by_genome_class.tsv"
HITS_TSV = OUT_DIR + "men_hits_in_windows.tsv"
LOCI_TSV = OUT_DIR + "duf4229_loci.tsv"
RESULTS_MD = OUT_DIR + "RESULTS_TABLES.md"

GENOMEWIDE_DIR = "_data/evoweaver/genomewide/"
KEPT_TXT = GENOMEWIDE_DIR + "proteomes_kept.txt"
PROTEOMES_TSV = GENOMEWIDE_DIR + "proteomes.tsv"
GENES_DIR = GENOMEWIDE_DIR + "genes/"

CLASSICAL = ["MenA", "MenE", "MenB", "MenD", "MenC", "MenF", "MenH", "MenG"]
OTHERS = ["MenJ", "MqnP", "MqnE", "MqnC", "UbiA"]   # UbiA = Pfam PF01040 at GA, any UbiA-family prenyltransferase
CLASS_LABELS = ["MenA", "MenE", "MenB", "MenD", "any_classical", "MqnE", "MqnC", "MqnP", "UbiA"]
ANCHOR = ("UP000000582", "Q8NT63")          # C. glutamicum ATCC 13032, Cgl0447 = cg0530
SEED = 4229
CONTROLS_PER_LOCUS = 5                      # as in build_windows.py

# product-name patterns, used only to compare with the HMM calls
TEXT_PATTERNS = {
    "MenA": r"naphthoate (octa|poly)?prenyltransferase",
    "MenE": r"succinylbenzo(ate|ic acid)[- ]+coa ligase",
    "MenB": r"naphtho(yl|ate)[- ]coa synthase|naphthoate synthase",
    "MenD": r"2-succinyl-5-enolpyruvyl|sephchc",
    "MenC": r"o-succinylbenzoate synthase",
    "MenF": r"isochorismate synthase",
    "MenH": r"2-succinyl-6-hydroxy|shchc",
    "MenG": r"demethylmenaquinone methyltransferase|menaquinone biosynthesis (c-)?methyltransferase|ubiquinone/menaquinone",
}


def log(msg=""):
    print(msg, flush=True)


def read_tsv(path):
    lines = open(path).read().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def proteomes_in(name):
    """Proteomes in a family table of _data/evoweaver/genomewide/genes/."""
    found = set()
    for line in open(GENES_DIR + name + ".tsv").read().splitlines()[1:]:
        if line.count("\t") == 4:
            found.update(p.strip().split(":")[0] for p in line.split("\t")[3].split(";") if p.strip().startswith("UP"))
    return found


# ---------------------------------------------------------------- statistics

def log_add(a, b):
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    hi, lo = max(a, b), min(a, b)
    return hi + math.log1p(math.exp(lo - hi))


def log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_log10_p(a, n1, c, n2):
    """One-sided Fisher exact test that group 1 (a of n1) has the higher rate than group 2 (c of n2)."""
    total, hits = n1 + n2, a + c
    log_p = -math.inf
    for x in range(a, min(n1, hits) + 1):
        log_p = log_add(log_p, log_choose(n1, x) + log_choose(n2, hits - x) - log_choose(total, hits))
    return log_p / math.log(10)


def matched_log10_p(sets):
    """Exact upper tail of sum of Bernoulli(t / (k + 1)) over matched sets, as log10.

    sets is a list of (duf_has_gene, t, k). Sets with t = 0 add nothing and sets
    with t = k + 1 add a fixed 1, so only 0 < t < k + 1 enter the convolution.
    """
    observed = sum(d for d, t, k in sets)
    fixed = sum(1 for d, t, k in sets if t == k + 1)
    probs = [t / (k + 1) for d, t, k in sets if 0 < t < k + 1]
    dist = [0.0]                         # log P(sum = i) over the informative sets
    for p in probs:
        lp, lq = math.log(p), math.log(1 - p)
        new = [-math.inf] * (len(dist) + 1)
        for i, v in enumerate(dist):
            if v == -math.inf:
                continue
            new[i] = log_add(new[i], v + lq)
            new[i + 1] = log_add(new[i + 1], v + lp)
        dist = new
    need = observed - fixed
    log_p = -math.inf
    for i in range(max(need, 0), len(dist)):
        log_p = log_add(log_p, dist[i])
    expected = fixed + sum(probs)
    return log_p / math.log(10), observed, expected


def fmt_p(log10_p):
    if log10_p == -math.inf:
        return "0"
    if log10_p > -3:
        return f"{10 ** log10_p:.3g}"
    exponent = math.floor(log10_p)
    return f"{10 ** (log10_p - exponent):.1f}e{exponent}"


def pct(a, n):
    return f"{a:,}/{n:,} ({100 * a / n:.1f}%)" if n else f"{a}/{n}"


# ---------------------------------------------------------------- data

def load_windows(path, calls, text_patterns):
    """window_id -> dict with kind, locus, proteome, members [(offset, rel_strand, protein_id, gene, text_genes)]."""
    windows = {}
    for row in read_tsv(path):
        w = windows.setdefault(row["window_id"], {
            "id": row["window_id"], "kind": row["kind"], "locus": row["proteome"] + "|" + row["duf_acc"],
            "proteome": row["proteome"],
            "duf_acc": row["duf_acc"], "contig": row["contig"], "center_strand": row["center_strand"],
            "members": []})
        call = calls.get(row["protein_id"], {})
        named = (row["product"] + " ").lower()
        text_genes = [g for g, pattern in text_patterns.items() if re.search(pattern, named)]
        w["members"].append({
            "offset": int(row["offset"]), "rel_strand": row["rel_strand"], "protein_id": row["protein_id"],
            "locus_tags": row["locus_tags"], "product": row["product"], "gene_name": row["gene"],
            "call": call.get("gene", ""), "ubia": call.get("ubia_pfam_score", "") != "",
            "curated_menB": call.get("curated_menB_pass") == "1", "curated_menD": call.get("curated_menD_pass") == "1",
            "text": text_genes})
    return windows


def genes_in(window):
    """label -> closest member (offset != 0) carrying it.

    Labels are the HMM gene calls, "UbiA" for the Pfam context call, "text:<gene>"
    for product-name matches, and "any_classical" for any of the eight CLASSICAL genes.
    """
    found = {}
    for m in sorted(window["members"], key=lambda m: (abs(m["offset"]), m["offset"])):
        if m["offset"] == 0:
            continue
        labels = [m["call"]] if m["call"] else []
        if m["call"] in CLASSICAL:
            labels.append("any_classical")
        if m["ubia"]:
            labels.append("UbiA")
        labels += ["text:" + g for g in m["text"]]
        for label in labels:
            found.setdefault(label, m)
    return found


def orientation(offset, rel_strand):
    if rel_strand == "same":
        return "co-directional"
    return "convergent" if offset > 0 else "divergent"


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()

    targets = read_tsv(TARGETS_TSV)
    lineage = {t["proteome"]: t for t in targets}
    calls = {r["protein_id"]: r for r in read_tsv(CALLS_TSV)}
    windows = load_windows(WINDOWS_TSV, calls, TEXT_PATTERNS)
    genome = {r["proteome"]: r for r in read_tsv(GENOME_TSV)}

    duf = {w["locus"]: w for w in windows.values() if w["kind"] == "duf"}
    controls = defaultdict(list)
    for w in windows.values():
        if w["kind"] == "control":
            controls[w["locus"]].append(w)
    found = {wid: genes_in(w) for wid, w in windows.items()}
    duf_found = {locus: found[w["id"]] for locus, w in duf.items()}
    n_duf = len(duf)
    n_ctrl = sum(len(v) for v in controls.values())
    target_genomes = {t["proteome"] for t in targets}
    placed_genomes = {w["proteome"] for w in duf.values()}
    log(f"{len(targets):,} DUF4229 proteins in {len(target_genomes):,} genomes without MqnE; "
        f"{n_duf:,} placed on their contigs ({len(placed_genomes):,} genomes); {n_ctrl:,} control windows")
    full = 2 * max(abs(m["offset"]) for w in duf.values() for m in w["members"]) + 1
    for kind in ("duf", "control"):
        group = [w for w in windows.values() if w["kind"] == kind]
        short = sum(1 for w in group if len(w["members"]) < full)
        log(f"  {kind} windows cut short by a contig end: {pct(short, len(group))}")
    log(f"  loci with fewer than {CONTROLS_PER_LOCUS} control windows: "
        f"{sum(1 for l in duf if len(controls[l]) < CONTROLS_PER_LOCUS):,}")

    # one locus per genus, for the dereplicated test
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
            ctrl = controls[locus]
            t = d + sum(1 for w in ctrl if label in found[w["id"]])
            sets.append((d, t, len(ctrl)))
        return sets

    # ---- per gene
    summary = []
    all_loci = sorted(duf)
    for label in CLASSICAL + ["any_classical"] + OTHERS + ["text:" + g for g in CLASSICAL]:
        a, n1, c, n2 = count(label, all_loci)
        log_p = fisher_log10_p(a, n1, c, n2)
        m_log_p, observed, expected = matched_log10_p(matched_sets(label, all_loci))
        ga, gn1, gc, gn2 = count(label, genus_loci)
        g_log_p = fisher_log10_p(ga, gn1, gc, gn2)
        rate1, rate2 = a / n1, (c / n2 if n2 else 0)
        enrichment = rate1 / rate2 if rate2 else math.inf
        hits = [duf_found[locus][label] for locus in all_loci if label in duf_found[locus]]
        offsets = Counter(m["offset"] for m in hits)
        same = sum(1 for m in hits if m["rel_strand"] == "same")
        adjacent = Counter(orientation(m["offset"], m["rel_strand"]) for m in hits if abs(m["offset"]) == 1)
        top = ", ".join(f"{o:+d}: {n}" for o, n in offsets.most_common(4))
        summary.append({
            "gene": label, "duf_windows": a, "duf_total": n1, "control_windows": c, "control_total": n2,
            "duf_rate": rate1, "control_rate": rate2, "enrichment": enrichment,
            "fisher_log10_p": log_p, "matched_log10_p": m_log_p, "matched_expected": expected,
            "genus_duf": ga, "genus_total": gn1, "genus_control": gc, "genus_control_total": gn2,
            "genus_fisher_log10_p": g_log_p, "top_offsets": top, "same_strand": same, "n_hits": len(hits),
            "adjacent": dict(adjacent), "within_2": sum(1 for m in hits if abs(m["offset"]) <= 2)})

    with open(SUMMARY_TSV, "w") as out:
        out.write("gene\tduf_windows_with_gene\tduf_windows\tcontrol_windows_with_gene\tcontrol_windows\t"
                  "duf_rate\tcontrol_rate\tenrichment\tfisher_log10_p\tmatched_log10_p\tmatched_expected\t"
                  "genus_duf_with_gene\tgenera\tgenus_control_with_gene\tgenus_control_windows\tgenus_fisher_log10_p\t"
                  "top_offsets\tsame_strand\tadjacent_orientation\twithin_2_genes\n")
        for s in summary:
            out.write(f"{s['gene']}\t{s['duf_windows']}\t{s['duf_total']}\t{s['control_windows']}\t{s['control_total']}\t"
                      f"{s['duf_rate']:.4f}\t{s['control_rate']:.5f}\t{s['enrichment']:.1f}\t{s['fisher_log10_p']:.1f}\t"
                      f"{s['matched_log10_p']:.1f}\t{s['matched_expected']:.1f}\t{s['genus_duf']}\t{s['genus_total']}\t"
                      f"{s['genus_control']}\t{s['genus_control_total']}\t{s['genus_fisher_log10_p']:.1f}\t"
                      f"{s['top_offsets']}\t{s['same_strand']}\t"
                      f"{';'.join(f'{k}={v}' for k, v in sorted(s['adjacent'].items()))}\t{s['within_2']}\n")

    log()
    log("gene      DUF4229 windows        control windows       enrichment  Fisher p   matched p  "
        "genera (DUF4229 vs control)        top offsets")
    for s in summary:
        log(f"{s['gene']:<9} {pct(s['duf_windows'], s['duf_total']):<22} {pct(s['control_windows'], s['control_total']):<21} "
            f"{s['enrichment']:>8.1f}x  {fmt_p(s['fisher_log10_p']):>9}  {fmt_p(s['matched_log10_p']):>9}  "
            f"{pct(s['genus_duf'], s['genus_total'])} vs {pct(s['genus_control'], s['genus_control_total'])}   {s['top_offsets']}")

    # ---- every men hit in a window, for the offset and strand figures
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

    # ---- genome level: MenD and MenB anywhere in the genome
    kept = set(open(KEPT_TXT).read().split())
    with_mqne = proteomes_in("MqnE") & kept
    with_duf = proteomes_in("DUF4229") & kept
    with_mqnc = proteomes_in("MqnC") & kept
    with_pf02621 = proteomes_in("MqnAD") & kept
    futalosine_duf = with_duf & with_mqne

    def classical_share(group):
        group = [up for up in group if up in genome]
        n = sum(1 for up in group if genome[up]["classical"] == "1")
        return n, len(group)

    log()
    log("genome level (UniProtKB InterPro matches to NCBIfam MenD TIGR00173 and MenB TIGR01929/NF006186/NF005637)")
    status = Counter()
    for up in target_genomes:
        g = genome[up]
        status[(int(g["n_MenD"]) > 0, int(g["n_MenB"]) > 0)] += 1
    log(f"  DUF4229 without MqnE: MenD and MenB {pct(status[(True, True)], len(target_genomes))}; "
        f"MenD only {status[(True, False)]}; MenB only {status[(False, True)]}; neither {status[(False, False)]}")
    n, d = classical_share(futalosine_duf)
    log(f"  DUF4229 with MqnE (futalosine): MenD and MenB {pct(n, d)}")
    log(f"  DUF4229 without MqnE that carry MqnC (TIGR03699): {len(target_genomes & with_mqnc)}; "
        f"PF02621 (MqnA/MqnD): {len(target_genomes & with_pf02621)}")
    for gene in ("MenA", "MenE", "MqnP"):
        has = sum(1 for up in target_genomes if int(genome[up][f"n_{gene}"]) > 0)
        log(f"  DUF4229 without MqnE with {gene} anywhere in the genome: {pct(has, len(target_genomes))}")

    # how often each kind of Actinomycetota genome carries DUF4229 at all
    actino = set()
    for line in open(PROTEOMES_TSV).read().splitlines()[1:]:
        parts = line.split("\t")
        if parts[0] in kept and "Actinomycetota" in parts[2]:
            actino.add(parts[0])
    kinds = {"MqnE (futalosine)": {up for up in actino if up in with_mqne},
             "no MqnE, MenD and MenB (classical)": {up for up in actino - with_mqne if genome[up]["classical"] == "1"},
             "no MqnE, MenD or MenB only": {up for up in actino - with_mqne if genome[up]["classical"] == "0"
                                            and (int(genome[up]["n_MenD"]) or int(genome[up]["n_MenB"]))},
             "no MqnE, neither": {up for up in actino - with_mqne
                                  if not int(genome[up]["n_MenD"]) and not int(genome[up]["n_MenB"])}}
    log(f"  kept Actinomycetota genomes carrying DUF4229, by pathway ({len(actino):,} genomes):")
    for name, ups in kinds.items():
        log(f"    {name}: {pct(len(ups & with_duf), len(ups))}")

    # ---- the same counts split by what the genome carries
    def genome_class(up):
        d, b = int(genome[up]["n_MenD"]) > 0, int(genome[up]["n_MenB"]) > 0
        if d and b:
            return "classical (MenD and MenB)"
        if d or b:
            return "MenD or MenB only"
        if up in with_mqnc or up in with_pf02621:
            return "neither, futalosine genes present"
        return "neither"

    by_class = defaultdict(list)
    for locus in all_loci:
        by_class[genome_class(duf[locus]["proteome"])].append(locus)
    log()
    log("by genome class (DUF4229 windows against their own control windows)")
    with open(CLASSES_TSV, "w") as out:
        out.write("genome_class\tgenomes\tgene\tduf_windows_with_gene\tduf_windows\tcontrol_windows_with_gene\t"
                  "control_windows\tfisher_log10_p\tmatched_log10_p\n")
        for name, loci in sorted(by_class.items(), key=lambda x: -len(x[1])):
            ups = {duf[l]["proteome"] for l in loci}
            log(f"  {name}: {len(ups):,} genomes, {len(loci):,} loci")
            for label in CLASS_LABELS:
                a, n1, c, n2 = count(label, loci)
                log_p = fisher_log10_p(a, n1, c, n2)
                m_log_p = matched_log10_p(matched_sets(label, loci))[0]
                out.write(f"{name}\t{len(ups)}\t{label}\t{a}\t{n1}\t{c}\t{n2}\t{log_p:.1f}\t{m_log_p:.1f}\n")
                log(f"    {label:<14} {pct(a, n1):<22} control {pct(c, n2):<20} Fisher p {fmt_p(log_p):>9}  "
                    f"matched p {fmt_p(m_log_p):>9}")

    # ---- per order
    by_order = defaultdict(list)
    for locus, w in duf.items():
        by_order[lineage[w["proteome"]]["order"] or "(none)"].append(locus)
    with open(ORDERS_TSV, "w") as out:
        out.write("order\tloci\tgenera\tMenA_in_window\tMenE_in_window\tany_classical_in_window\t"
                  "control_any_classical\tcontrol_windows\tgenomes_with_MenD_and_MenB\tgenomes\n")
        log()
        log("order                  loci  genera  MenA in window   MenE in window   any men in window   control any men   classical genomes")
        for order, loci in sorted(by_order.items(), key=lambda x: -len(x[1])):
            genera = len({lineage[duf[l]['proteome']]['genus'] for l in loci})
            a = sum(1 for l in loci if "MenA" in duf_found[l])
            e = sum(1 for l in loci if "MenE" in duf_found[l])
            anyc = sum(1 for l in loci if "any_classical" in duf_found[l])
            ctrl = [w for l in loci for w in controls[l]]
            cany = sum(1 for w in ctrl if "any_classical" in found[w["id"]])
            ups = {duf[l]["proteome"] for l in loci}
            cl = sum(1 for up in ups if genome[up]["classical"] == "1")
            out.write(f"{order}\t{len(loci)}\t{genera}\t{a}\t{e}\t{anyc}\t{cany}\t{len(ctrl)}\t{cl}\t{len(ups)}\n")
            log(f"{order:<22} {len(loci):>5}  {genera:>6}  {pct(a, len(loci)):<15}  {pct(e, len(loci)):<15}  "
                f"{pct(anyc, len(loci)):<18}  {pct(cany, len(ctrl)):<16}  {pct(cl, len(ups))}")

    # ---- per locus table
    with open(LOCI_TSV, "w") as out:
        out.write("proteome\tduf_acc\tgenus\torder\tspecies\tcontig\tclassical_genome\t" +
                  "\t".join(f"{g}_offset" for g in CLASSICAL + OTHERS) + "\tcontrol_windows\n")
        for locus in all_loci:
            w = duf[locus]
            t = lineage[w["proteome"]]
            f = duf_found[locus]
            cells = [f"{f[g]['offset']:+d}{'s' if f[g]['rel_strand'] == 'same' else 'o'}" if g in f else ""
                     for g in CLASSICAL + OTHERS]
            out.write(f"{w['proteome']}\t{w['duf_acc']}\t{t['genus']}\t{t['order']}\t{t['species']}\t{w['contig']}\t"
                      f"{genome[w['proteome']]['classical']}\t" + "\t".join(cells) + f"\t{len(controls[locus])}\n")

    # ---- anchor locus
    log()
    locus = "|".join(ANCHOR)
    if locus in duf:
        w = duf[locus]
        log(f"anchor: {lineage[w['proteome']]['species']} {w['proteome']} {w['duf_acc']} on {w['contig']} "
            f"(DUF4229 strand {w['center_strand']})")
        for m in sorted(w["members"], key=lambda m: m["offset"]):
            if abs(m["offset"]) <= 6:
                log(f"  {m['offset']:+3d} {m['rel_strand']:<8} {m['locus_tags']:<14} {m['call'] or '-':<5} "
                    f"{'UbiA' if m['ubia'] else '':<5} {m['product']}")
    else:
        log(f"anchor {locus} not placed")

    # ---- HMM calls against product names, and NCBIfam against the curated July profiles
    log()
    log("HMM call against product-name text, over every protein in every window")
    agreement = []
    seen = set()
    for w in windows.values():
        for m in w["members"]:
            if m["protein_id"] in seen or not m["protein_id"]:
                continue
            seen.add(m["protein_id"])
            agreement.append(m)
    for g in CLASSICAL:
        both = sum(1 for m in agreement if m["call"] == g and g in m["text"])
        hmm_only = sum(1 for m in agreement if m["call"] == g and g not in m["text"])
        text_only = sum(1 for m in agreement if m["call"] != g and g in m["text"])
        log(f"  {g}: both {both:,}, HMM only {hmm_only:,}, text only {text_only:,}")
    for g, key in (("MenB", "curated_menB"), ("MenD", "curated_menD")):
        both = sum(1 for m in agreement if m["call"] == g and m[key])
        nc_only = sum(1 for m in agreement if m["call"] == g and not m[key])
        cu_only = sum(1 for m in agreement if m["call"] != g and m[key])
        log(f"  {g} NCBIfam against the curated July profile: both {both:,}, NCBIfam only {nc_only:,}, curated only {cu_only:,}")

    write_results_md(summary, status, target_genomes, classical_share(futalosine_duf), n_duf, n_ctrl, len(genus_loci))
    return 0


def write_results_md(summary, status, target_genomes, futalosine, n_duf, n_ctrl, n_genera):
    with open(RESULTS_MD, "w") as out:
        out.write("# Men genes within 10 genes of DUF4229, genomes without MqnE\n\n")
        out.write(f"{n_duf:,} DUF4229 loci, {n_ctrl:,} control windows, {n_genera} genera (one locus each).\n\n")
        out.write("| gene | DUF4229 windows | control windows | enrichment | Fisher p | matched p | genera, DUF4229 | genera, control | top offsets | same strand |\n")
        out.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for s in summary:
            out.write(f"| {s['gene']} | {pct(s['duf_windows'], s['duf_total'])} | {pct(s['control_windows'], s['control_total'])} | "
                      f"{s['enrichment']:.1f} | {fmt_p(s['fisher_log10_p'])} | {fmt_p(s['matched_log10_p'])} | "
                      f"{pct(s['genus_duf'], s['genus_total'])} | {pct(s['genus_control'], s['genus_control_total'])} | "
                      f"{s['top_offsets']} | {s['same_strand']}/{s['n_hits']} |\n")
        out.write(f"\nGenomes without MqnE carrying MenD and MenB: {pct(status[(True, True)], len(target_genomes))}. "
                  f"DUF4229 genomes with MqnE carrying MenD and MenB: {pct(*futalosine)}.\n")


if __name__ == "__main__":
    sys.exit(main())
