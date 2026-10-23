#!/bin/bash
# Retrains the combiner on neighborhood negatives, end to end. Every stage is resumable.
#   caffeinate -i nohup bash _scripts/ensemble/run_neighbor_model.sh > _data/ensemble/neighbor_model.log 2>&1 &
set -e
cd "$(dirname "$0")/../.."
stamp() { date "+%Y-%m-%d %H:%M:%S"; }

echo "[$(stamp)] stage 1: neighborhood families for the eight training systems"
python3 _scripts/ensemble/pick_system_neighbors.py

echo "[$(stamp)] stage 2: labels and the pair table"
python3 _scripts/ensemble/build_neighbor_pairs.py

echo "[$(stamp)] stage 3: score, retrain, evaluate, apply once to futalosine"
Rscript _scripts/ensemble/train_neighbor_model.R

echo "[$(stamp)] finished"
