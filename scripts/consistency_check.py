#!/usr/bin/env python3
"""
Cross-document consistency check.

Scans every project markdown file for figures that were SUPERSEDED by the
2026-07-20 audit, and for internal contradictions (e.g. a corrected count
sitting next to a stale p-value). Run this before any submission.

Exit code 1 if any stale figure is found.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Files that are historical records and are ALLOWED to contain old numbers.
EXEMPT = {"NOTEBOOK.md", "AUDIT_REPORT_2026-07-20.md",
          "WEEKLY_REPORT_2026-07-13_to_07-20.md",
          # dated record, superseded by BIWEEKLY_REPORT_2026-07-20_to_08-04.md
          "WEEKLY_REPORT_2026-07-20_to_24.md"}

# A stale figure is ACCEPTABLE on a line that explicitly marks it as historical
# or as pending recomputation. These markers are how we keep honest provenance
# without the checker crying wolf.
ALLOW_MARKERS = re.compile(
    r"earlier|superseded|supersede|previously|label-based|not yet recomputed|"
    r"PURGE|stray/wrong|old (value|figure)|was 27 genera|historical", re.I)

# (regex, why it is stale, what to use instead)
STALE = [
    (r"5\.5% of (vitamin K|futalosine)", "old label-based co-occurrence rate", "9.4% (342/3,648)"),
    (r"200\s*/\s*3,?647", "old numerator/denominator", "342/3,648"),
    (r"3\.8\s*[×x]\s*10⁻²⁶", "old Fisher p", "3.9 × 10⁻⁴⁵"),
    (r"3\.8e-26", "old Fisher p", "3.9e-45"),
    (r"~?140\s*[×x]\b", "old odds ratio", "~251×"),
    (r"0\s*/\s*1,?206", "old background denominator", "0/1,211"),
    (r"1,?206[- ](enzyme|genome|background)", "old background size", "1,211"),
    (r"9\.5% of MqnE", "old MqnE rate", "85.1% (280/329, actinobacteria)"),
    (r"1\.8% of MqnC", "old MqnC rate", "9.7% (55/566, actinobacteria)"),
    # The audit-era label-based MqnE/MqnC split, superseded by the sequence-based
    # classifier (NOTEBOOK 22). These are one generation newer than the 9.5%/1.8%
    # figures above and were still live in README and the Philipp email until 2026-07-25.
    (r"278\s*/\s*333", "label-based MqnE split", "280/329 = 85.1%"),
    (r"83\.5%", "label-based MqnE rate", "85.1% (280/329)"),
    (r"6\.3e-114|6\.3\s*[×x]\s*10⁻¹¹⁴", "label-based MqnE p-value", "4.1e-120"),
    (r"same[- ]strand (in )?12/12", "old operon sample", "565/566 (99.8%)"),
    (r"12/12 DUF4229", "old operon sample", "565/566 pairs"),
    (r"27 (distinct|independent) genera", "old genus count", "45 genera"),
    (r"0 of 774 background genera", "superseded by taxonomy-matched control", "0/212 actinobacterial background"),
    (r"(?<![\d,])3,628(?![\d])|(?<![\d,])3,678(?![\d])", "stray/incorrect intermediate count", "3,648 or 3,651"),
]

# Numbers that MUST co-occur consistently within a single file.
# (corrected_value, stale_value, corrected_counterpart, why)
# A line holding the stale value AND its corrected counterpart is a
# before/after comparison and is not a contradiction.
PAIRED = [
    (r"342\s*/\s*3,?648", r"3\.8e-26|3\.8\s*[×x]\s*10⁻²⁶", r"3\.9e-45|3\.9\s*[×x]\s*10⁻⁴⁵",
     "corrected count next to OLD p-value"),
    (r"342\s*/\s*3,?648", r"140\s*[×x]", r"251\s*[×x]",
     "corrected count next to OLD odds ratio"),
    (r"45 distinct genera", r"1 genus of 27|27 independent", r"45 genera|45 distinct genera",
     "corrected genus count next to OLD genus claim"),
]



def check_naked_rates(files, root):
    """Flag HEADLINE percentages that appear without their numerator/denominator.

    Rationale: the one substantive analytical error in this project was quoting
    a subgroup rate (52.8%, 93.3%) that came from a positive-enriched stratified
    sample without reweighting it. A bare percentage hides its denominator and
    therefore hides its sampling design. A percentage reported as "280/329 =
    85.1%" cannot hide that. This flags emphasised (bold) percentages that carry
    no counts on the same line, which are the ones a reader will quote.
    """
    # the percentage must sit INSIDE a single bold span, with no other
    # bold delimiters between, so "**A** ... 50% ... **B**" does not match
    pct_bold = re.compile(r"\*\*(?:(?!\*\*).)*?\d+(?:\.\d+)?\s*%(?:(?!\*\*).)*?\*\*")
    # counts may be written "342/3,648" OR "342 of 3,648"
    frac = re.compile(r"\d[\d,]*\s*(?:/|\bof\b)\s*\d[\d,]*")
    # percentages used as PARAMETERS (thresholds, cutoffs) are not rates
    param = re.compile(r"threshold|cutoff|identity|co-occurrence threshold|window", re.I)
    hits = []
    for f in files:
        if f in EXEMPT:
            continue
        for i, line in enumerate(open(os.path.join(root, f), errors="replace"), 1):
            if (pct_bold.search(line) and not frac.search(line)
                    and not param.search(line) and not ALLOW_MARKERS.search(line)):
                hits.append((f, i, line.strip()[:96]))
    print("\n" + "=" * 72)
    print("NAKED HEADLINE RATE CHECK (bold % with no numerator/denominator)")
    print("=" * 72)
    if not hits:
        print("  clean: every emphasised percentage carries its counts")
    else:
        for f, i, t in hits:
            print(f"  {f}:{i}  {t}")
        print(f"\n  {len(hits)} emphasised percentage(s) without counts.")
        print("  Add the fraction (e.g. '280/329 = 85.1%') so the denominator,")
        print("  and therefore the sampling design, is visible.")
    return len(hits)


def check_rate_arithmetic(files, root):
    """Verify that a stated numerator/denominator actually equals its percentage.

    Rationale: the existing naked-rate check enforces that an emphasised
    percentage CARRIES counts, but not that the counts AGREE with it. That gap
    let a half-completed update through on 2026-07-24: when the firmicute lead
    numbers were restandardized onto the two-way classifier (NOTEBOOK 46), the
    percentages and p-values in Results 3.6 were updated to the new classifier
    but the numerators were left on the old four-way one, producing
    "294 of 468 = 63.7%" (294/468 is 62.8%). Every other live document had been
    updated, so nothing contradicted it and no stale-figure pattern matched.
    Arithmetic is the only check that catches a half-applied correction.
    """
    # "342/3,648 = 9.38%" and "280 of 329 = 85.1%"
    fwd = re.compile(r"(\d[\d,]*)\s*(?:/|\bof\b)\s*(\d[\d,]*)\s*=\s*(\d+(?:\.\d+)?)\s*%")
    # "9.38% (342/3,648)"
    rev = re.compile(r"(\d+(?:\.\d+)?)\s*%\s*\(\s*(\d[\d,]*)\s*(?:/|\bof\b)\s*(\d[\d,]*)\s*\)")
    # "250 of 342 positive genomes (73.1%)" and "73.1% of positives (250 of 342)".
    # These shapes slipped a half-applied correction through on 2026-08-04: the
    # fractions were updated by search-and-replace but the inline percentages in
    # parentheses several words away were not, and neither pattern above matched.
    par = re.compile(r"(\d[\d,]*)\s*(?:/|\bof\b)\s*(\d[\d,]*)\b[^()%]{0,30}?\(\s*(\d+(?:\.\d+)?)\s*%\s*\)")
    # No second "%" may appear in the gap, or a line carrying two statistics pairs
    # the first percentage with the second fraction and reports a phantom mismatch.
    parrev = re.compile(r"(\d+(?:\.\d+)?)\s*%[^()%]{0,30}?\(\s*(\d[\d,]*)\s*(?:/|\bof\b)\s*(\d[\d,]*)\s*\)")
    hits = []
    for f in files:
        if f in EXEMPT:
            continue
        for i, line in enumerate(open(os.path.join(root, f), errors="replace"), 1):
            # Quoting a known-wrong fraction in order to document its correction is
            # legitimate, and is marked the same way as elsewhere in this checker.
            if ALLOW_MARKERS.search(line):
                continue
            found = ([(m.group(1), m.group(2), m.group(3)) for m in fwd.finditer(line)] +
                     [(m.group(2), m.group(3), m.group(1)) for m in rev.finditer(line)] +
                     [(m.group(1), m.group(2), m.group(3)) for m in par.finditer(line)] +
                     [(m.group(2), m.group(3), m.group(1)) for m in parrev.finditer(line)])
            for num_s, den_s, pct_s in found:
                num = int(num_s.replace(",", ""))
                den = int(den_s.replace(",", ""))
                if den == 0:
                    continue
                stated = float(pct_s)
                actual = 100.0 * num / den
                # tolerate rounding to the precision the percentage is written at
                dec = len(pct_s.split(".")[1]) if "." in pct_s else 0
                tol = 0.5 * (10 ** -dec) + 1e-6
                if abs(actual - stated) > tol:
                    hits.append((f, i, f"{num_s}/{den_s} = {pct_s}%", actual))
    print("\n" + "=" * 72)
    print("RATE ARITHMETIC CHECK (does numerator/denominator equal the stated %?)")
    print("=" * 72)
    if not hits:
        print("  clean: every stated fraction matches its percentage")
    else:
        for f, i, txt, actual in hits:
            print(f"  {f}:{i}  '{txt}'  but {txt.split('=')[0].strip()} = {actual:.2f}%")
        print(f"\n  {len(hits)} fraction(s) disagree with their own percentage.")
        print("  This is the signature of a half-applied correction: check whether")
        print("  the numerator or the percentage is the one that was updated.")
    return len(hits)


def main():
    problems = 0
    files = sorted(f for f in os.listdir(ROOT) if f.endswith(".md"))
    print("=" * 72)
    print("STALE FIGURE SCAN")
    print("=" * 72)
    for f in files:
        if f in EXEMPT:
            continue
        path = os.path.join(ROOT, f)
        text = open(path, errors="replace").read()
        hits = []
        for pat, why, use in STALE:
            for m in re.finditer(pat, text):
                line = text[:m.start()].count("\n") + 1
                line_text = text.split("\n")[line - 1]
                if ALLOW_MARKERS.search(line_text):
                    continue          # deliberate historical reference
                # A line carrying BOTH the stale value and its replacement is a
                # before/after comparison, which is legitimate.
                use_key = re.sub(r"[^0-9a-zA-Z.]", "", use.split("(")[0])
                if use_key and use_key in re.sub(r"[^0-9a-zA-Z.]", "", line_text):
                    continue
                hits.append((line, m.group(0), why, use))
        if hits:
            problems += len(hits)
            print(f"\n  {f}")
            for line, txt, why, use in sorted(hits):
                print(f"    L{line:<5} '{txt}'  [{why}] -> use {use}")
    if not problems:
        print("\n  no stale figures found")

    print("\n" + "=" * 72)
    print("INTERNAL CONTRADICTION SCAN")
    print("=" * 72)
    contra = 0
    for f in files:
        if f in EXEMPT:
            continue
        text = open(os.path.join(ROOT, f), errors="replace").read()
        for a, b, newpat, why in PAIRED:
            # only a contradiction if the stale value appears on a line that is
            # NOT marked as a historical reference
            stale_lines = [ln for ln in text.split("\n")
                           if re.search(b, ln) and not ALLOW_MARKERS.search(ln)
                           and not re.search(a, ln)
                           and not re.search(newpat, ln)]  # both values = comparison
            if re.search(a, text) and stale_lines:
                print(f"  {f}: {why}")
                for ln in stale_lines[:2]:
                    print(f"      {ln.strip()[:110]}")
                contra += 1
    if not contra:
        print("  no internal contradictions found")

    naked = check_naked_rates(files, ROOT)
    arith = check_rate_arithmetic(files, ROOT)

    print("\n" + "=" * 72)
    print("EM DASH CHECK (project style rule)")
    print("=" * 72)
    bad = 0
    for f in files:
        n = open(os.path.join(ROOT, f), errors="replace").read().count("—")
        if n:
            print(f"  {f}: {n} em dashes")
            bad += n
    if not bad:
        print("  clean")

    total = problems + contra + bad + naked + arith
    print("\n" + "=" * 72)
    print(f"TOTAL ISSUES: {total}")
    print("=" * 72)
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
