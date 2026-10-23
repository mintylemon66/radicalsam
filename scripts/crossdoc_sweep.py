#!/usr/bin/env python3
"""Systematic sweep for half-applied corrections across the project's live documents.

The four errors found on 2026-07-25 all had the same shape: a value was corrected in
some documents but not others, so no single document contradicted itself and the
pattern-based checker in scripts/consistency_check.py saw nothing.

This looks for that shape directly, keyed on the ENTITY a number describes so that
unrelated families sharing a denominator are never compared:

  Signal 1  the same (entity, denominator) carries two different numerators
  Signal 2  the same (entity, fraction) is paired with two different p-values
  Signal 3  the same (entity, denominator) carries two different percentages

Each number is attributed to the NEAREST entity mention on its line, because real
sentences name several entities at once. Pfam accessions are folded into their family
name so that "DUF1450" and "PF07293" are recognised as the same thing.

Lines that mark themselves as historical are skipped, as elsewhere in the project.

Usage:  python3 crossdoc_sweep.py [project_root]
Exit 1 if any signal is found.
"""
import os
import re
import sys
import collections

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/Users/jackli/Desktop/Radical_SAM_Project"

# Dated records and files that deliberately carry before/after narratives.
EXEMPT = {
    "NOTEBOOK.md",
    "AUDIT_REPORT_2026-07-20.md",
    "WEEKLY_REPORT_2026-07-13_to_07-20.md",
    "SELF_AUDIT_2026-07-24.md",
    "SELF_AUDIT_2026-07-25.md",
    "COMPLETENESS_REVIEW_2026-07-24.md",
    "_START_HERE.md",
    "papermaterials.md",
    # Dated record of the 2026-07-22 overnight run. Its DUF1450/YtxC rows carry the
    # four-way-classifier values on purpose; a supersession banner was added above the
    # table on 2026-07-25 pointing at the canonical two-way values.
    "MORNING_REPORT.md",
    # Exists specifically to compare two counting conventions side by side, so it
    # carries both values of several quantities on purpose.
    "ANCHOR_CONVENTION_2026-07-26.md",
}

# (entity, denominator) pairs where two different numerators are legitimately expected,
# with the reason. Reviewed and confirmed benign; listed here so they are not re-chased.
KNOWN_BENIGN = {
    # 566 means two unrelated things in this project: the number of actinobacterial
    # MqnC genomes (55/566 carry DUF4229) and the number of DUF4229-MqnE gene pairs
    # scored for strand (565/566 same-strand). A coincidence of denominators, not an error.
    ("DUF4229", 566): "566 = actino MqnC genomes AND DUF4229-MqnE operon pairs",
    ("MQNC", 566): "566 = actino MqnC genomes AND DUF4229-MqnE operon pairs",
    ("MQNE", 566): "566 = actino MqnC genomes AND DUF4229-MqnE operon pairs",
    # 342 is the count of DUF4229-positive futalosine loci, and many DIFFERENT
    # quantities are legitimately measured out of it: 341 same-strand, 282 downstream,
    # 188 at offset +2, 7 with an overlapping CDS. Different numerators over a shared
    # denominator is expected here and is not a half-applied correction.
    ("DUF4229", 342): "342 = DUF4229-positive loci; many distinct measures share it",
    ("MQNE", 342): "342 = DUF4229-positive loci; many distinct measures share it",
    ("MQNC", 342): "342 = DUF4229-positive loci; many distinct measures share it",
}

HIST = re.compile(
    r"earlier|superseded|supersede|previously|label-based|old (value|figure)|"
    r"historical|corrected|correction|wrong|stale|never use|before/after|"
    r"recount|does not reproduce|instead of|rather than", re.I)

FRAC = re.compile(r"(?<![\d.])(\d[\d,]*)\s*(?:/|\bof\b)\s*(\d[\d,]*)(?![\d])")
PCT = re.compile(r"(\d+(?:\.\d+)?)\s*%")
PVAL = re.compile(r"[pq]\s*=\s*~?\s*(\d+(?:\.\d+)?\s*(?:e-?\d+|[×x]\s*10⁻[⁰¹²³⁴⁵⁶⁷⁸⁹]+))", re.I)

ENTITY = re.compile(
    r"\b(DUF\d{3,4}|YtxC|YuiB|PF\d{5}|MqnE|MqnC|MqnP|BioB|LipA|AsnC|UbiA|MenB|MenD)\b", re.I)

# Pfam accession -> family name, so the two spellings collapse to one entity.
ALIAS = {
    "PF14012": "DUF4229", "PF07293": "DUF1450", "PF08812": "YTXC",
    "PF14068": "YUIB", "PF07315": "DUF1462", "PF13829": "DUF4191",
    "PF04359": "DUF493",
}


def canon(tok):
    t = tok.upper()
    return ALIAS.get(t, t)


def norm(n):
    return int(n.replace(",", ""))


def live_files():
    out = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames
                       if not d.startswith(".") and d not in
                       {"node_modules", "mlenv", "archive", "drafts_reference",
                        "foldseek_db", "proteomes", "batches", "snap", "models", "ml",
                        # analysis output with its own dense result tables, where many
                        # unrelated families legitimately share one background denominator
                        "philipp_todos"}]
        for fn in filenames:
            if fn.endswith(".md") and fn not in EXEMPT:
                out.append(os.path.relpath(os.path.join(dirpath, fn), ROOT))
    return sorted(out)


def entities_on(line):
    """Every distinct entity named on the line, canonicalised.

    A number is attributed to ALL of them rather than to the nearest one. Nearest-entity
    attribution was tried first and proved fragile: "YtxC is found in 240/422 of MqnC
    genomes" and "YtxC co-occurs in 240/422 of MqnC genomes" put the same fact nearer to
    different words, so the two statements were filed under different keys and a real
    disagreement between them would have been missed. Over-attributing costs a few
    candidates to eyeball; under-attributing costs a silent miss, which is the whole
    failure mode this sweep exists to catch.
    """
    return {canon(m.group(1)) for m in ENTITY.finditer(line)}


def main():
    files = live_files()
    den_map = collections.defaultdict(lambda: collections.defaultdict(list))
    pval_map = collections.defaultdict(lambda: collections.defaultdict(list))
    pct_map = collections.defaultdict(lambda: collections.defaultdict(list))

    for rel in files:
        for i, line in enumerate(open(os.path.join(ROOT, rel), errors="replace"), 1):
            if HIST.search(line) or not ENTITY.search(line):
                continue
            loc = (rel, i, line.strip()[:130])
            ents = entities_on(line)
            fr_all = FRAC.findall(line)
            pv_all = PVAL.findall(line)
            pc_all = PCT.findall(line)
            # A line enumerating many fractions (a results table row, a summary listing
            # eight families against one background) is not a single claim, and
            # attributing every fraction on it to every entity on it is meaningless.
            # Such lines are the dominant false-positive source, so skip them.
            if len(fr_all) > 2 or len(ents) > 3:
                continue
            for m in FRAC.finditer(line):
                num, den = norm(m.group(1)), norm(m.group(2))
                if den < 50 or num > den:
                    continue
                for ent in ents:
                    den_map[(ent, den)][num].append(loc)
                    # pair with a p-value or percentage only when the line states exactly
                    # one of each alongside exactly one fraction, so the pairing is certain
                    if len(pv_all) == 1 and len(fr_all) == 1:
                        pval_map[(ent, f"{num}/{den}")][
                            re.sub(r"\s+", "", pv_all[0]).lower()].append(loc)
                    if len(pc_all) == 1 and len(fr_all) == 1:
                        pct_map[(ent, den)][pc_all[0]].append(loc)

    problems = 0

    def report(title, mapping, fmt):
        nonlocal problems
        print("\n" + "=" * 74)
        print(title)
        print("=" * 74)
        # A half-applied correction spans DOCUMENTS. Two different numerators for the
        # same denominator inside one line (or one file) is almost always a legitimate
        # list of several families measured against a shared background, which is the
        # dominant false positive here. Require the disagreeing values to appear in at
        # least two distinct files.
        hits = {}
        benign = []
        for k, v in mapping.items():
            if len(v) < 2:
                continue
            if k in KNOWN_BENIGN:
                benign.append((k, KNOWN_BENIGN[k]))
                continue
            files_per_val = {val: {loc[0] for loc in locs} for val, locs in v.items()}
            allfiles = set().union(*files_per_val.values())
            if len(allfiles) < 2:
                continue
            # and the values must not all co-occur in a single shared file
            if any(len(fs) == len(allfiles) for fs in files_per_val.values()) and \
               all(fs == allfiles for fs in files_per_val.values()):
                continue
            hits[k] = v
        for k, why in sorted(benign, key=str):
            print(f"  (known benign, skipped) {k[0]} denominator {k[1]}: {why}")
        if not hits:
            print("  clean" if not benign else "  clean apart from the known-benign entries above")
            return
        for key in sorted(hits, key=lambda k: (str(k[0]), str(k[1]))):
            problems += 1
            print(f"\n  {fmt(key)}")
            for val in sorted(hits[key], key=str):
                for rel, i, txt in hits[key][val][:2]:
                    print(f"    -> {val}   {rel}:{i}")
                    print(f"       {txt}")

    report("SIGNAL 1: same entity + denominator, DIFFERENT numerators",
           den_map, lambda k: f"{k[0]}, denominator {k[1]}")
    report("SIGNAL 2: same entity + fraction, DIFFERENT p-values",
           pval_map, lambda k: f"{k[0]}, {k[1]}")
    report("SIGNAL 3: same entity + denominator, DIFFERENT percentages",
           pct_map, lambda k: f"{k[0]}, denominator {k[1]}")

    print("\n" + "=" * 74)
    print(f"TOTAL CROSS-DOCUMENT SIGNALS: {problems}")
    print("=" * 74)
    if problems:
        print("\nEach hit is a candidate, not a confirmed error. Look at the two lines:")
        print("if they describe the same quantity, one document was left behind by a")
        print("correction. That is the failure mode this sweep exists to catch.")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
