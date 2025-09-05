#!/bin/bash
# ============================================
# Train GCNN model using Chemprop
#
# NOTE:
#   The --data-path argument specifies the training data file (CSV format).
#   You can change this to any of your own processed datasets, for example:
#     --data-path "data/processed/data_cmnpd_after2000.csv"
#     --data-path "data/processed/trainingset.csv"
#
# Usage:
#   1) Run with the default data path:
#        ./train.sh
#
#   2) Run with a custom data path (pass it as the first argument):
#        ./train.sh data/processed/trainingset.csv
# ============================================

DATA_PATH="../data/processed/data_cmnpd_after2000.csv"

if [ ! -z "$1" ]; then
  DATA_PATH="$1"
fi

chemprop train \
    --data-path "$DATA_PATH" \
    --task-type multiclass \
    --multiclass-num-classes 3 \
    --epochs 200 \
    --batch-size 32 \
    --output-dir "model_after2000" \
    --num-workers 0 \
    --ffn-hidden-dim 1100 \
    --ffn-num-layers 3 \
    --message-hidden-dim 1100 \
    --depth 3 \
    --num-replicates 5\