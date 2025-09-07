# 🔍 Error Data Detection and Cleaning

This directory contains scripts and notebooks for **error data detection and cleaning**, corresponding to the relevant parts of the paper.  
The purpose of these scripts is to identify mislabeled or noisy samples, evaluate their impact on model training, and generate a cleaned dataset for the development of the final model.

## 📂 Files

- **`cross_training.py`**  
  Implements multiple **cross-training** experiments.  
  - Constructs stratified train/val/test splits;  
  - Repeats splitting and training **n = 3 / 5 / 10** times;  
  - Records performance fluctuations and misclassified samples across repeats;  
  - In the paper, **Repeat = 10** was finally chosen as the basis for cleaning.

- **`detection.ipynb`**  
  Identifies potential mislabels based on the workflow **error prediction list → similarity search in the training set**:  
  1) **Input**: a list of misclassified samples from cross-training/evaluation;  
  2) **Similarity search**: computes Morgan fingerprints (RDKit) for each misclassified sample, and retrieves highly similar molecules (default **Tanimoto ≥ 0.7**, adjustable) from the training set.  

- **`error_prediction_and_cleaning.ipynb`**  
  Integrates results from repeated cross-training, using the **10-repeat models**:  
  - Aggregates samples that were **consistently misclassified** across multiple runs;  
  - Produces the cleaned dataset used in the main experiments of the paper.

## 📝 Notes

- The analysis shows that **label noise** has a significant impact on natural product origin classification;  
- Using **Repeat = 10** more reliably exposes problematic samples.