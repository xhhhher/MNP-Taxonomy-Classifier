#!/bin/bash
# ============================================
# Predict with trained GCN model
#
# NOTE:
#   This script runs Chemprop prediction.
#   You can specify custom paths for:
#     1) Test set CSV file (with SMILES)
#     2) Output predictions file
#     3) Trained model checkpoint (.ckpt)
#
# Usage:
#   1) Run with defaults:
#        ./predict_chemprop.sh
#
#   2) Run with custom arguments:
#        ./predict_chemprop.sh <test_path> <preds_path> <model_path>
#
# Example:
#   ./predict_chemprop.sh ../data/processed/final_test_set_smiles.csv \
#                         ../results/predictions/final_test_set_predict.csv \
#                         ../results/models/model_GCN_cleaned_finetuned.ckpt
# ============================================

TEST_PATH="../data/processed/CMNPD2.0_test_set_smiles.csv"
PREDS_PATH="../results/predictions/CMNPD2.0_test_set_predict.csv"
MODEL_PATH="../results/models/model_GCN_cleaned_finetuned.ckpt"

if [ ! -z "$1" ]; then
  TEST_PATH="$1"
fi
if [ ! -z "$2" ]; then
  PREDS_PATH="$2"
fi
if [ ! -z "$3" ]; then
  MODEL_PATH="$3"
fi

echo "Running prediction..."
echo "Test set: $TEST_PATH"
echo "Predictions will be saved to: $PREDS_PATH"
echo "Using model: $MODEL_PATH"

chemprop predict \
  --test-path "$TEST_PATH" \
  --preds-path "$PREDS_PATH" \
  --model-paths "$MODEL_PATH" \
  --num-workers 0