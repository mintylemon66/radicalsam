#!/bin/bash
# Runs the whole ensemble pipeline end to end. Every stage is resumable, so if this is interrupted
# you can run it again and it picks up where it stopped.
#   caffeinate -i nohup bash _scripts/ensemble/run_all.sh > _data/ensemble/run.log 2>&1 &
set -e
cd "$(dirname "$0")/../.."
export PYTHONPATH=_scripts/ensemble:_scripts/evoweaver/genomewide
stamp() { date "+%Y-%m-%d %H:%M:%S"; }

echo "[$(stamp)] stage 1: fetch families from UniProt"
python3 _scripts/ensemble/fetch_systems.py

echo "[$(stamp)] stage 2: fetch contig feature tables from NCBI"
python3 _scripts/ensemble/fetch_system_contigs.py

echo "[$(stamp)] stage 3: build EvoWeaver labels and the pair table"
python3 _scripts/ensemble/build_training_pairs.py

echo "[$(stamp)] stage 4: score, fit the combiner, evaluate, score the held-out system"
Rscript _scripts/ensemble/train_ensemble.R

echo "[$(stamp)] finished"
