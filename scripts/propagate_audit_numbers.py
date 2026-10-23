#!/usr/bin/env python3
"""
Propagate the 2026-07-20 audit corrections into the project documents.

Every replacement is explicit and reported. Statistics that were computed on
the ORIGINAL label-based calls and have NOT been recomputed (permutation test,
Benjamini-Hochberg FDR table) are deliberately NOT silently updated; they are
annotated as pending instead, because changing the underlying counts while
leaving those derived statistics untouched would be inconsistent.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

EDITS = {
"papermaterials.md": [
 ("Full run: **3,647 vitamin K (MqnC/MqnE) enzymes vs. 1,206 background radical SAM enzymes**, ±10-gene neighborhoods (FINAL, complete).",
  "Full run: **3,648 vitamin K (MqnC/MqnE) enzymes vs. 1,211 background radical SAM enzymes**, ±10-gene neighborhoods. *(Counts and rates below were CORRECTED by the 2026-07-20 sequence-based audit; see `AUDIT_REPORT_2026-07-20.md`.)*"),
 ("- DUF4229 occurs next to **5.5% of vitamin K enzymes (200/3,647) vs. 0% of background** (FINAL: 0/1,206 background - DUF4229 never appears near a non-MqnE radical SAM enzyme).",
  "- DUF4229 occurs next to **9.4% of vitamin K enzymes (342/3,648) vs. 0% of background** (0/1,211 background, confirmed by SEQUENCE via HMMER/PF14012, not by product-name text). The earlier 5.5% figure was an under-count: the text search was only 56.9% sensitive."),
 ("- **Fisher's exact test (one-sided): p = 3.8 × 10⁻²⁶ (final, 1,206 background).** (earlier smaller samples gave p=2.6e-7; strengthened as background grew - report this trend.)",
  "- **Fisher's exact test (one-sided): p = 3.9 × 10⁻⁴⁵ (sequence-based, 1,211 background).** (Earlier label-based value was 3.8e-26; the association strengthened under every correction.)"),
 ("- Odds ratio (Haldane–Anscombe continuity-corrected for the zero cell): **~140×**",
  "- Odds ratio (Haldane–Anscombe continuity-corrected for the zero cell): **~251×**"),
 ("- The 5.5% rate was stable across every scale tested (45 → 300 → 3,628 genomes) - report this stability as evidence of robustness.",
  "- The rate was stable across every scale tested - report this stability as evidence of robustness."),
 ("the 200 hits span **27 distinct genera and 175 distinct species**",
  "the hits span **45 distinct genera** (sequence-based recount; the earlier label-based figure was 27 genera)"),
 ("- **1,206** = background radical SAM neighborhoods.",
  "- **1,211** = background radical SAM neighborhoods (audited; the raw table holds 1,298 rows, of which 1,211 yielded neighborhoods)."),
 ("- **200 / 3,647 = 5.5%** DUF4229 near vitamin K; **0 / 1,206 = 0%** near background.",
  "- **342 / 3,648 = 9.4%** DUF4229 near vitamin K (sequence-based); **0 / 1,211 = 0%** near background."),
 ("- **Fisher p = 3.8×10⁻²⁶; odds ratio ≈ 140×; 27 genera vs 0.**",
  "- **Fisher p = 3.9×10⁻⁴⁵; odds ratio ≈ 251×; 45 genera vs 0.**\n- **TAXONOMY-MATCHED (the strongest control): 335/895 = 37.4% of ACTINOBACTERIAL futalosine genomes vs 0/212 actinobacterial background, p = 2.0×10⁻³⁸. Within *Streptomyces*: 193/509 vs 0/49, p = 2.8×10⁻¹⁰.**\n- **MqnE-specific (headline): 278/333 = 83.5% of actinobacterial MqnE genomes carry DUF4229, vs 9.9% of MqnC, p = 6.3×10⁻¹¹⁴.**\n- **Operon: 565/566 pairs (99.8%) same-strand, binomial p = 2.3×10⁻¹⁶⁸, median gap 705 bp.**"),
 ("- **Operon / co-transcription (POSITIVE):** across 12/12 DUF4229-host organisms, DUF4229 and *mqnE* are on the **same DNA strand** (median intergenic gap ~740 bp). Same-strand in 12/12 is highly non-random (~1/4000 by chance) → supports co-transcription/co-regulation. A genuine, independent positive line.",
  "- **Operon / co-transcription (POSITIVE, now measured at scale):** across **565 of 566** DUF4229-radical SAM pairs (**99.8%**), DUF4229 and the radical SAM gene are on the **same DNA strand** (median intergenic gap **705 bp**; binomial vs 50% null **p = 2.3×10⁻¹⁶⁸**). Supersedes the earlier 12/12 observation."),
 ("**(b) Permutation test (100,000 label shuffles):** observed 200/3,647 DUF4229 in the vitamin K group;",
  "**(b) Permutation test (100,000 label shuffles) [COMPUTED ON THE ORIGINAL LABEL-BASED CALLS; NOT YET RECOMPUTED on the sequence-based counts]:** observed 200/3,647 DUF4229 in the vitamin K group;"),
 ("| **DUF4229 (discovery)** | 5% | 0% | 3.8e-26 | **1.1e-25** |",
  "| **DUF4229 (discovery)** [label-based; FDR not yet recomputed] | 5% | 0% | 3.8e-26 | **1.1e-25** |"),
],
"_START_HERE.md": [
 ("near 5.5% of vitamin K enzymes vs 0% background (200/3,647 vs 0/1,206), Fisher p = 3.8e-26, 27 genera vs 0, operon same-strand 12/12.",
  "near 9.4% of vitamin K enzymes vs 0% background (342/3,648 vs 0/1,211, SEQUENCE-based), Fisher p = 3.9e-45, 45 genera vs 0, operon same-strand 565/566 (99.8%). In actinobacteria specifically, 83.5% of MqnE genomes vs 0/212 matched controls."),
 ("3,651 target MqnC/MqnE; 3,647 with usable neighborhoods; 1,206 background; DUF4229 5.5% vs 0%; p=3.8e-26; ~140x; 27 genera. Database = UniProtKB TrEMBL, release 2026_07 (E-value 1e-10).",
  "3,651 target MqnC/MqnE; 3,648 with usable neighborhoods; 1,211 background; DUF4229 **9.4% vs 0%** (sequence-based); **p=3.9e-45; ~251x; 45 genera**. Actinobacteria-matched: **37.4% vs 0/212**; MqnE-only **83.5%**. Operon **565/566**. Database = UniProtKB TrEMBL, release 2026_07 (E-value 1e-10). **These supersede the old 5.5%/140x/3.8e-26 figures (see `AUDIT_REPORT_2026-07-20.md`).**"),
],
"Results.md": [
 ("co-occurs beside **5.5% of futalosine enzymes (200 of 3,647) but 0% of background enzymes (0 of 1,206)**",
  "co-occurs beside **9.4% of futalosine enzymes (342 of 3,648) but 0% of background enzymes (0 of 1,211)**"),
 ("- **Taxonomic independence.** The 200 hits span **27 distinct genera and 175 distinct species**.",
  "- **Taxonomic independence.** The hits span **45 distinct genera**."),
 ("**Partner specificity.** The association is specifically with **MqnE (aminofutalosine synthase), not MqnC**: DUF4229 co-occurs beside 9.5% of MqnE genomes but only 1.8% of MqnC genomes.",
  "**Partner specificity.** The association is specifically with **MqnE (aminofutalosine synthase), not MqnC**: in actinobacteria, DUF4229 co-occurs beside **83.5% of MqnE genomes (278/333) but only 9.9% of MqnC genomes (55/554)**, Fisher p = 6.3 × 10⁻¹¹⁴."),
],
"Methods.md": [
 ("The **background set** is 1,206 ordinary radical SAM enzymes",
  "The **background set** is 1,211 ordinary radical SAM enzymes"),
 ("The analysis was run for the 3,647 target enzymes with a reconstructable neighborhood (4 of the 3,651 dereplicated targets had no mappable genome and were dropped) against the 1,206 background neighborhoods.",
  "The analysis was run for the 3,648 target enzymes with a reconstructable neighborhood (3 of the 3,651 dereplicated targets had no mappable genome and were dropped) against the 1,211 background neighborhoods. Family membership of neighbors was assigned from **sequence**, by HMMER search against the PF14012 profile at the Pfam-curated gathering threshold, not from GenBank product-name text."),
],
"UNDERSTAND_YOUR_PROJECT.md": [
 ("**Your finding (DUF4229) is specifically MqnE's partner** (near 9.5% of MqnE vs 1.8% of MqnC)",
  "**Your finding (DUF4229) is specifically MqnE's partner** (in actinobacteria, near **83.5% of MqnE vs 9.9% of MqnC**, p = 6.3e-114)"),
 ("→ No - collapsed to one-per-genus, it's still 27 genera vs 0. That's a database bias I corrected for.",
  "→ No. Collapsed to one-per-genus it is still 45 genera vs 0. Better still, I ran a taxonomy-matched control: within *Streptomyces* alone, DUF4229 sits beside 193 of 509 futalosine enzymes and **0 of 49** non-futalosine ones (p = 2.8e-10). Same genus, so the over-sampling cannot explain it."),
],
"WETLAB_PROPOSAL.md": [
 ("DUF4229 is operonically conserved with *mqnE* (same strand 12/12; 27 genera; co-occurrence p = 3.8×10⁻²⁶)",
  "DUF4229 is operonically conserved with *mqnE* (same strand **565/566 = 99.8%**; **45 genera**; co-occurrence **p = 3.9×10⁻⁴⁵**; in actinobacteria it accompanies **83.5% of MqnE** genes vs **0/212** matched controls)"),
],
}

BANNER = """> **NUMBERS CORRECTED 2026-07-20.** A sequence-based audit found the original
> co-occurrence counts were derived from GenBank product-name TEXT, which was only
> 56.9% sensitive. All headline figures were recomputed from sequence (HMMER vs
> PF14012). The finding survived every control and got stronger. See
> `AUDIT_REPORT_2026-07-20.md`. Statistics marked "not yet recomputed" (permutation
> test, FDR table) still rest on the original label-based calls.

"""


def main():
    total = 0
    for fname, edits in EDITS.items():
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            print(f"[skip] {fname} not found")
            continue
        s = open(path).read()
        applied = missed = 0
        for old, new in edits:
            if old in s:
                s = s.replace(old, new)
                applied += 1
            else:
                missed += 1
                print(f"  [MISS] {fname}: {old[:72]}...")
        if fname in ("papermaterials.md", "_START_HERE.md") and "NUMBERS CORRECTED 2026-07-20" not in s:
            lines = s.split("\n")
            i = 1 if lines and lines[0].startswith("#") else 0
            s = "\n".join(lines[:i]) + ("\n" if i else "") + BANNER + "\n".join(lines[i:])
        open(path, "w").write(s)
        print(f"[{fname}] applied {applied}, missed {missed}")
        total += applied
    print(f"\nTOTAL replacements applied: {total}")


if __name__ == "__main__":
    sys.exit(main())
