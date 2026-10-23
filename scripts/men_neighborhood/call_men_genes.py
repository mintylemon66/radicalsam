#!/usr/bin/env python3
"""Call classical menaquinone (men) genes in the window proteins with profile HMMs.

    python3 _scripts/men_neighborhood/call_men_genes.py          # after fetch_window_proteins.py
    python3 _scripts/men_neighborhood/call_men_genes.py --futalosine   # genomes WITH MqnE, in FUTALOSINE_DIR

Every call is a sequence call. The script downloads the NCBIfam models in
MEN_MODELS (NCBI's hmm_PGAP release) into HMM_DIR and runs hmmsearch with
--cut_ga, so a protein only counts when it passes the model's own gathering
threshold. A protein that passes models of two genes (for example MenA and the
futalosine-branch prenyltransferase MqnP, both UbiA family) takes the gene of
its best-scoring model. The Pfam UbiA family model (PF01040) is searched the
same way as context, outside the men calls.

The curated MenB and MenD profiles built in July (mencheck/, NOTEBOOK.md
section 43) have no gathering threshold, so they are applied with the rule from
mencheck/menbd_check.py: a protein passes when it scores at least
CURATED_TOLERANCE times the weaker of the two positive controls (MenB or MenD of
Mycobacterium tuberculosis H37Rv and Corynebacterium glutamicum ATCC 13032).
They serve as a cross-check on the NCBIfam MenB and MenD calls.

With --futalosine every file, the downloaded models included, is read from and
written to FUTALOSINE_DIR (the DUF4229 genomes that carry MqnE).
"""

import argparse
import os
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict

OUT_DIR = "_data/men_neighborhood/"
PROTEINS_FAA = OUT_DIR + "window_proteins.faa"
HMM_DIR = OUT_DIR + "hmm/"
MEN_HMM = HMM_DIR + "men_models.hmm"
NCBIFAM_TBL = OUT_DIR + "hmmsearch_men_models.tbl"
CALLS_TSV = OUT_DIR + "men_calls.tsv"
MODELS_TSV = OUT_DIR + "men_models.tsv"
FUTALOSINE_DIR = OUT_DIR + "futalosine/"      # --futalosine: the same files for the DUF4229 genomes with MqnE

NCBIFAM_URL = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/"
UBIA_HMM = "_hmm/cat/PF01040.hmm"
CURATED_HMMS = {"MenB": "mencheck/menB.hmm", "MenD": "mencheck/menD.hmm"}
POSITIVE_CONTROLS = "mencheck/positive_control_menABD.faa"
CURATED_TOLERANCE = 0.97
CURATED_REPORT_EVALUE = "1e-5"

# (gene, NCBIfam model, what NCBIfam calls it); every model whose taxonomic range covers Actinomycetota
MEN_MODELS = [
    ("MenA", "TIGR00751.1", "1,4-dihydroxy-2-naphthoate octaprenyltransferase"),
    ("MenA", "NF004751.0", "1,4-dihydroxy-2-naphthoate polyprenyltransferase (PRK06080)"),
    ("MenA", "NF004752.0", "1,4-dihydroxy-2-naphthoate polyprenyltransferase (PRK06080, Bacteria)"),
    ("MenE", "TIGR01923.1", "o-succinylbenzoate--CoA ligase"),
    ("MenE", "NF005877.1", "o-succinylbenzoate--CoA ligase (PRK07824)"),
    ("MenE", "NF002966.0", "o-succinylbenzoate--CoA ligase (PRK03640)"),
    ("MenB", "TIGR01929.1", "1,4-dihydroxy-2-naphthoyl-CoA synthase"),
    ("MenB", "NF006186.0", "1,4-dihydroxy-2-naphthoyl-CoA synthase (PRK08321)"),
    ("MenB", "NF005637.0", "1,4-dihydroxy-2-naphthoyl-CoA synthase (PRK07396)"),
    ("MenD", "TIGR00173.1", "SEPHCHC synthase"),
    ("MenC", "TIGR01927.1", "o-succinylbenzoate synthase (gamma-proteobacteria and Gram-positive type)"),
    ("MenC", "TIGR01928.1", "o-succinylbenzoate synthase (low GC and archaeal type)"),
    ("MenC", "NF002782.0", "o-succinylbenzoate synthase (PRK02901)"),
    ("MenC", "NF002739.0", "o-succinylbenzoate synthase (PRK02714)"),
    ("MenF", "TIGR00543.1", "isochorismate synthase"),
    ("MenH", "TIGR03695.1", "SHCHC synthase"),
    ("MenG", "TIGR01934.1", "ubiquinone/menaquinone biosynthesis methyltransferase"),
    ("MenG", "NF001240.0", "UbiE/MenG methyltransferase (PRK00216.1-1)"),
    ("MenG", "NF001241.0", "demethylmenaquinone methyltransferase (PRK00216.1-2)"),
    ("MenG", "NF001242.0", "UbiE/MenG methyltransferase (PRK00216.1-3)"),
    ("MenG", "NF001243.0", "demethylmenaquinone methyltransferase (PRK00216.1-4)"),
    ("MenG", "NF001244.0", "UbiE/MenG methyltransferase (PRK00216.1-5)"),
    ("MenJ", "NF045655.1", "menaquinone reductase MenJ (Mycobacteriaceae)"),
    # futalosine-branch genes, to see whether a window the genome annotation calls MqnE-free still holds one
    ("MqnP", "NF041584.1", "futalosine-branch prenyltransferase MqnP (Actinomycetes)"),
    ("MqnE", "TIGR03700.1", "aminofutalosine synthase MqnE"),
    ("MqnC", "TIGR03699.1", "cyclic dehypoxanthinyl futalosine synthase MqnC"),
]


def log(msg):
    print(msg, flush=True)


def in_folder(path, futalosine):
    """path itself, or the file of the same name in FUTALOSINE_DIR for the --futalosine run."""
    return FUTALOSINE_DIR + path[len(OUT_DIR):] if futalosine else path


def download_models(hmm_dir, men_hmm):
    """Fetch each NCBIfam HMM once into hmm_dir and concatenate them, plus Pfam UbiA, into men_hmm."""
    os.makedirs(hmm_dir, exist_ok=True)
    paths = []
    for gene, model, _ in MEN_MODELS:
        path = hmm_dir + model + ".HMM"
        if not os.path.exists(path):
            with urllib.request.urlopen(NCBIFAM_URL + model + ".HMM", timeout=300) as response:
                data = response.read()
            with open(path + ".part", "wb") as out:
                out.write(data)
            os.rename(path + ".part", path)
            log(f"  downloaded {model}")
        paths.append(path)
    paths.append(UBIA_HMM)
    with open(men_hmm, "w") as out:
        for path in paths:
            out.write(open(path).read())
    return paths


def model_cutoffs(path):
    """NAME, ACC and the sequence GA threshold of every model in an HMM file."""
    models, current = [], {}
    for line in open(path):
        if line.startswith("NAME "):
            current = {"name": line.split()[1]}
        elif line.startswith("ACC "):
            current["acc"] = line.split()[1]
        elif line.startswith("GA "):
            current["ga"] = float(line.split()[1].rstrip(";"))
        elif line.startswith("//"):
            models.append(current)
    return models


def run_hmmsearch(hmm, fasta, tbl, threshold_args, cpus):
    command = ["hmmsearch", "--noali", "--cpu", str(cpus), "--tblout", tbl] + threshold_args + [hmm, fasta]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)


def read_tblout(path):
    """(target, query_name, query_acc, evalue, score, description) for every line of a --tblout file."""
    hits = []
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split(None, 18)
        hits.append((p[0], p[2], p[3], float(p[4]), float(p[5]), p[18] if len(p) > 18 else ""))
    return hits


def curated_cut(gene, cpus, folder):
    """Score cut for a curated July profile: CURATED_TOLERANCE x the weaker positive control."""
    tbl = folder + f"hmmsearch_curated_{gene}_controls.tbl"
    run_hmmsearch(CURATED_HMMS[gene], POSITIVE_CONTROLS, tbl, ["-E", CURATED_REPORT_EVALUE], cpus)
    controls = [score for target, _, _, _, score, desc in read_tblout(tbl) if f"GN={gene[0].lower()}{gene[1:]} " in desc + " "]
    if len(controls) != 2:
        sys.exit(f"expected 2 positive controls for {gene}, found {len(controls)}")
    return CURATED_TOLERANCE * min(controls), controls


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--futalosine", action="store_true",
                        help=f"the DUF4229 genomes that carry MqnE; files are read from and written to {FUTALOSINE_DIR}")
    parser.add_argument("--proteins", help=f"default {PROTEINS_FAA}")
    parser.add_argument("--output", help=f"default {CALLS_TSV}")
    parser.add_argument("--cpus", type=int, default=os.cpu_count() or 4)
    args = parser.parse_args()
    if shutil.which("hmmsearch") is None:
        sys.exit("hmmsearch not on PATH; install HMMER 3")
    folder = in_folder(OUT_DIR, args.futalosine)
    proteins_faa = args.proteins or in_folder(PROTEINS_FAA, args.futalosine)
    calls_tsv = args.output or in_folder(CALLS_TSV, args.futalosine)
    men_hmm = in_folder(MEN_HMM, args.futalosine)
    ncbifam_tbl = in_folder(NCBIFAM_TBL, args.futalosine)

    log("NCBIfam models")
    download_models(in_folder(HMM_DIR, args.futalosine), men_hmm)
    gene_of = {model.split(".")[0]: gene for gene, model, _ in MEN_MODELS}
    gene_of["PF01040"] = "UbiA"
    cutoffs = model_cutoffs(men_hmm)
    with open(in_folder(MODELS_TSV, args.futalosine), "w") as out:
        out.write("gene\tmodel\tname\tgathering_threshold\tdescription\n")
        described = {model.split(".")[0]: what for _, model, what in MEN_MODELS}
        described["PF01040"] = "Pfam UbiA prenyltransferase family (context, not a men call)"
        for m in cutoffs:
            acc = m["acc"].split(".")[0]
            out.write(f"{gene_of[acc]}\t{m['acc']}\t{m['name']}\t{m['ga']}\t{described[acc]}\n")
    log(f"  {len(cutoffs)} models in {men_hmm}, every one with a GA line: {all('ga' in m for m in cutoffs)}")

    log("hmmsearch --cut_ga")
    run_hmmsearch(men_hmm, proteins_faa, ncbifam_tbl, ["--cut_ga"], args.cpus)
    hits = defaultdict(list)            # protein -> [(score, gene, model, evalue)]
    for target, name, acc, evalue, score, _ in read_tblout(ncbifam_tbl):
        hits[target].append((score, gene_of[acc.split(".")[0]], acc, evalue))

    log("curated July MenB and MenD profiles")
    curated = {}
    for gene in CURATED_HMMS:
        cut, controls = curated_cut(gene, args.cpus, folder)
        tbl = folder + f"hmmsearch_curated_{gene}.tbl"
        run_hmmsearch(CURATED_HMMS[gene], proteins_faa, tbl, ["-E", CURATED_REPORT_EVALUE], args.cpus)
        curated[gene] = {target: score for target, _, _, _, score, _ in read_tblout(tbl)}
        curated[gene + "_cut"] = cut
        log(f"  {gene}: controls {controls[0]:.1f} and {controls[1]:.1f} bits, cut {cut:.1f} bits; "
            f"{sum(s >= cut for s in curated[gene].values()):,} window proteins pass")

    proteins = set(hits) | set(curated["MenB"]) | set(curated["MenD"])
    with open(calls_tsv, "w") as out:
        out.write("protein_id\tgene\tmodel\tscore\tevalue\tall_men_hits\tubia_pfam_score\t"
                  "curated_menB_score\tcurated_menB_pass\tcurated_menD_score\tcurated_menD_pass\n")
        n_called = defaultdict(int)
        for pid in sorted(proteins):
            men = sorted((h for h in hits.get(pid, []) if h[1] != "UbiA"), reverse=True)
            ubia = max((h[0] for h in hits.get(pid, []) if h[1] == "UbiA"), default="")
            best = men[0] if men else ("", "", "", "")
            if men:
                n_called[best[1]] += 1
            listing = ";".join(f"{g}:{m}:{s:.1f}" for s, g, m, e in men)
            b, d = curated["MenB"].get(pid, ""), curated["MenD"].get(pid, "")
            out.write(f"{pid}\t{best[1]}\t{best[2]}\t{best[0]}\t{best[3]}\t{listing}\t{ubia}\t"
                      f"{b}\t{int(b != '' and b >= curated['MenB_cut'])}\t"
                      f"{d}\t{int(d != '' and d >= curated['MenD_cut'])}\n")
    log(f"wrote {calls_tsv}: " + ", ".join(f"{g} {n:,}" for g, n in sorted(n_called.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
