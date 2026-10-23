# DUF4229 sits beside classical menaquinone genes in Actinomycetota without MqnE

Run 18 October 2026. Every number below is in `logs/analyze_windows.log` and the `summary_*.tsv` files.
`latex/main_sts.tex` was not touched.

## Answer

Yes. In the 2,205 kept genomes that carry DUF4229 and no MqnE (2,203 Actinomycetota), DUF4229 sits next to
genes of the classical (Men) pathway, above all *menA*, the 1,4-dihydroxy-2-naphthoate octaprenyltransferase.

- **menA within 10 genes of DUF4229: 1,789 of 2,228 loci (80.3%)**, against 17 of 10,221 control windows (0.17%):
  483-fold, one-sided Fisher exact p = 10^-1704.5, matched-set exact p = 10^-1278.4. Both are far below the
  usual reporting floor, so write p < 10^-300.
- **Adjacent:** menA is the very next gene at 858 loci (38.5%; controls 2 of 10,221) and within two genes at
  1,235 (55.4%; controls 3 of 10,221).
- **Position and strand:** menA lies downstream of DUF4229 at 1,770 of 1,789 loci and on the opposite strand
  at 1,753 (98.0%). Where they are adjacent the two genes face each other (convergent, 853 of 858).
- **Across lineages:** with one locus per genus, menA is within 10 genes in 206 of 265 genera (77.7%), against
  1 of 1,252 control windows (p = 10^-198).
- **Any of the eight classical genes** (MenA, B, C, D, E, F, G, H) within 10 genes: 2,030 of 2,228 loci (91.1%)
  against 329 of 10,221 control windows (3.2%).
- **Most of these genomes run the classical pathway.** 2,035 of 2,205 (92.3%) encode both MenD and MenB, and 2,149
  (97.5%) encode MenA. Of the 931 DUF4229 genomes that do have MqnE, 6 (0.6%) encode MenD and MenB.
- **DUF4229 is common in both branches.** Among kept Actinomycetota genomes it is in 931 of 950 with MqnE
  (98.0%) and in 2,034 of 2,278 without MqnE that encode MenD and MenB (89.3%). That is why its presence does
  not track MqnE across genomes (the p = 0.21 null).
- **Anchor case reproduced.** *C. glutamicum* ATCC 13032 (proteome UP000000582, GenBank BA000036, same
  genome as BX927147): DUF4229 Cgl0447 (= cg0530) on the + strand, *menA* Cgl0448 (= cg0531) at +1 on the
  opposite strand, Cgl0449 (= cg0532) at +2, *menE* Cgl0450 (= cg0533) at +3 on the opposite strand. On the
  BX927147 proteins themselves, cg0531 scores 297 bits on TIGR00751 and cg0533 scores 564 bits on NF005877.

## Sentences she could use

Results, after the lineage explanation of the p = 0.21 null:

> Most of the 2,205 genomes without MqnE make menaquinone by the classical route instead: 2,035 of them
> encode MenD and MenB. DUF4229 lies within ten genes of \textit{menA} at 1,789 of their 2,228 DUF4229 loci
> (80.3\%), against 17 of 10,221 windows around length-matched genes from the same contigs (0.2\%; Fisher's
> exact test, $p < 10^{-300}$), most often as the next gene downstream, on the opposite strand and facing it.

Discussion:

> DUF4229 therefore sits near menaquinone genes in both branches of the pathway: near \textit{mqnE} where
> menaquinone is made by the futalosine route, and near \textit{menA}, as in \textit{C. glutamicum}, where it
> is made by the classical route.

("Near" on purpose: with MqnE, DUF4229 is most often two genes away, as the 17 October edit of
"immediately downstream" records.)

One wording point in the current text (line 196): counted from cg0530, *menE* (cg0533) is the third gene
(cg0531, cg0532, cg0533). "Two genes away" is right only when counted from *menA*.

## Counts per gene

DUF4229 windows (n = 2,228) against control windows (n = 10,221). "Genera" uses one DUF4229 locus per genus
(265 genera) and that locus's own control windows (1,252).

| gene | DUF4229 windows | control windows | fold | Fisher p | genera, DUF4229 vs control | usual position (strand relative to DUF4229) |
|---|---|---|---|---|---|---|
| MenA | 1,789 (80.3%) | 17 (0.17%) | 483 | 10^-1704.5 | 206/265 vs 1/1,252 | +1 (858), +2 (364), +3 (112); opposite 1,753/1,789 |
| MenB | 1,000 (44.9%) | 80 (0.8%) | 57 | 10^-727.2 | 151/265 vs 8/1,252 | +3 (370), +5 (184); opposite 974/1,000 |
| MenD | 828 (37.2%) | 102 (1.0%) | 37 | 10^-550.9 | 117/265 vs 9/1,252 | -2 (269), +6 to +9 |
| MenE | 740 (33.2%) | 92 (0.9%) | 37 | 10^-484.8 | 79/265 vs 4/1,252 | +2 to +4; opposite 717/740 |
| MenC | 701 (31.5%) | 79 (0.8%) | 41 | 10^-463.9 | 111/265 vs 10/1,252 | +4, +5; same 689/701 |
| MenG | 585 (26.3%) | 163 (1.6%) | 16.5 | 10^-309.5 | 93/265 vs 9/1,252 | +2 (153), +1 (64); same 351/585 |
| MenF | 392 (17.6%) | 58 (0.6%) | 31 | 10^-236.4 | 75/265 vs 4/1,252 | -4 (130) |
| MenH | 0 | 0 | | | | NCBIfam TIGR03695 matches no actinobacterial protein here |
| any of the 8 | 2,030 (91.1%) | 329 (3.2%) | 28 | 10^-1704.2 | 242/265 vs 24/1,252 | |

The matched-set test gives the same picture for every gene (MenA 10^-1278, MenB 10^-631, MenD 10^-493,
MenE 10^-458; all in `summary_by_gene.tsv`). Restricting to the 2,035 genomes that encode MenD and MenB:
menA within 10 genes at 1,680 of 2,056 loci (81.7%) against 17 of 9,491 control windows.

**Typical loci by order** (`summary_by_order.tsv`; arrows are strands relative to DUF4229):

| order | loci | genera | menA within 10 | any men gene within 10 | commonest arrangement |
|---|---|---|---|---|---|
| Micrococcales | 736 | 135 | 722 (98.1%) | 730 (99.2%) | menD< . DUF4229> menA< . menB< menC> |
| Mycobacteriales | 624 | 20 | 538 (86.2%) | 539 (86.4%) | DUF4229> . menA< (menA at +2 to +6) |
| Pseudonocardiales | 275 | 30 | 262 (95.3%) | 264 (96.0%) | DUF4229> . menA< . menE< menB< menD> |
| Streptosporangiales | 238 | 33 | 0 | 226 (95.0%) | DUF4229> . menG> (menG at +1 or +2, same strand) |
| Propionibacteriales | 183 | 8 | 161 (88.0%) | 161 (88.0%) | DUF4229> . menA< . . menC> menD> |

In Streptosporangiales the neighbor is *menG* (demethylmenaquinone methyltransferase), at +1 or +2 in 217 of
238 loci and followed by a geranylgeranyl reductase in 158 of those 217. Counting either *menA* or *menG*
within two genes: 1,452 of 2,228 loci (65.2%) against 38 of 10,221 controls (0.37%). These are the last two
steps of the classical pathway: MenA attaches the polyprenyl tail, and MenG methylates the resulting
membrane-bound demethylmenaquinone.

## Genomes "without MqnE" that are really futalosine genomes

55 of the 2,205 genomes (57 loci) encode neither MenD nor MenB but do carry the other futalosine genes
(MqnC in 49, MqnA/MqnD domain PF02621 in 55, MqnP in 54); 28 are *Streptomyces*, the rest mostly
Micromonosporales (*Paractinoplanes*, *Actinoplanes*, *Catenuloplanes*). Their MqnE is simply not matched in
UniProtKB. In them no classical gene sits near DUF4229 (2 of 57 windows, same as controls), and the
futalosine prenyltransferase *mqnP* is within 10 genes in 21 of 57 (36.8%, controls 2 of 236). So the
"2,203 Actinomycetota without MqnE" in the text includes about 55 futalosine genomes; the classical result
above does not depend on them.

## How each call was made

- **Gene order.** NCBI feature tables (efetch rettype=ft) for all 2,226 EMBL contigs carrying the 2,228
  DUF4229 proteins (3,906,510 CDS). Every DUF4229 protein was placed: 2,213 by locus tag, 13 by UniProtKB
  cross-reference, 1 by the GenBank gene qualifier (Cgl0447 on BA000036), 1 by EMBL protein ID (C. efficiens
  BAC17278.1). Window = 10 CDS either side; offsets run in DUF4229's direction of transcription.
- **Controls.** For each locus, 5 random genes from the same contig whose protein is within 25% of the DUF4229
  protein's length (DUF4229 median 102 aa) and more than 20 genes from any DUF4229, so no windows overlap
  (seed 4229). 10,221 control windows; 264 loci on short contigs had fewer than 5 eligible genes. Contig
  ends cut 8.8% of DUF4229 windows and 7.5% of control windows.
- **Men calls, windows.** All 230,421 window proteins from NCBI, searched with hmmsearch 3.4 `--cut_ga`
  (each model's own gathering threshold) using NCBIfam models from NCBI's hmm_PGAP release of 25 June 2026:
  MenA TIGR00751, NF004751, NF004752; MenE TIGR01923, NF005877, NF002966; MenB TIGR01929, NF006186, NF005637;
  MenD TIGR00173; MenC TIGR01927, TIGR01928, NF002782, NF002739; MenF TIGR00543; MenH TIGR03695; MenG
  TIGR01934, NF001240 to NF001244; plus MenJ NF045655, MqnP NF041584, MqnE TIGR03700, MqnC TIGR03699, and
  Pfam UbiA PF01040 as context. A protein passing models of two genes takes its best-scoring one.
- **Cross-checks.** The curated July MenB and MenD profiles (`mencheck/`, NOTEBOOK section 43; cut 0.97 x
  the weaker positive control: MenB 400.5 bits, MenD 419.7 bits) call 692 and 741 window proteins, every one
  of them also an NCBIfam call (the NCBIfam calls add 370 MenB and 172 MenD). Product names agree with the
  HMM calls for MenA (1,772 both, 29 HMM only, 5 name only), MenB (1,057, 5, 4) and MenD (872, 41, 13).
- **Genome level.** UniProtKB proteins of reference proteomes that InterPro matches to NCBIfam MenD
  TIGR00173, MenB TIGR01929/NF006186/NF005637, MenA TIGR00751, MenE TIGR01923/NF005877/NF002966 and
  MqnP NF041584 (InterProScan runs each model at its own cutoff). Classical = at least one MenD and one MenB.
- **Statistics.** One-sided Fisher exact test on windows; an exact matched-set permutation test (each DUF4229
  window against its own control windows from the same contig); the genus-level repeat with one random locus
  per genus. p-values are computed in log space (pure Python; checked against brute-force enumeration).

## Limits

- Genomic neighborhood shows linkage across genomes. It says nothing about contact or mechanism.
- The arrangement differs between branches. With MqnE, DUF4229 is on the same strand (341 of 342 pairs, the
  paper's operon). With *menA* it is on the opposite strand and convergent, so the two are unlikely to share
  one transcript; the conserved feature is adjacency.
- MenE at the gathering threshold is a floor. 374 more window proteins are named o-succinylbenzoate-CoA
  ligase, and 354 of them score 250 to 345 bits on NF005877, just under its threshold of 345; by product
  name, MenE is within 10 genes at 949 loci (42.6%). MenH cannot be called with NCBIfam models in these genomes.
- Loci from related genomes are not independent. The genus-level test covers part of this; the counts are
  still dominated by Micrococcales and Mycobacteriales (1,360 of 2,228 loci).
- "Without MqnE" rests on UniProtKB's TIGR03700 match (see the 55 futalosine genomes above).

## Files (`_data/men_neighborhood/`)

`README.md` lists every file and the run order. Main tables: `summary_by_gene.tsv`,
`summary_by_genome_class.tsv`, `summary_by_order.tsv`, `duf4229_loci.tsv` (one row per locus, offset and
strand of each gene), `men_hits_in_windows.tsv` (every call in every window). Scripts:
`_scripts/men_neighborhood/` (fetch_duf_contigs, fetch_genome_men, build_windows, fetch_window_proteins,
call_men_genes, analyze_windows).
