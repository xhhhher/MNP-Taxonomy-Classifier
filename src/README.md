# Source Code

This directory contains all scripts and notebooks for the **Marine Natural Product (MNP) Taxonomy Classifier** project.  
It is organized into submodules that correspond to different parts of the workflow described in the paper, from **data preprocessing**, **model training**, **error detection**, to **interpretability analysis**.

## Directory Structure

### 1. `cheminformatic_analysis/`
Scripts for **chemical space analysis**:  
- `properties_calculate.py` → Compute molecular properties (e.g., MW, LogP).  
- `properties_analysis.py` → Analyze property distributions.  
- `scaffold.py` → Scaffold overlap and visualization.  
- `README.md` → Notes on additional analysis (e.g., SOM with [DataWarrior](https://openmolecules.org/datawarrior/)).  

### 2. `train/`
Scripts for **training baseline and advanced models**:
- `train_chemprop.sh` → Train the GCN/MPNN model with Chemprop.  
- `train_last_ffn_SVM.py` → Train an SVM classifier on the last hidden layer (FFN) features.  
- `train_last_ffn_XGBoost.py` → Train an XGBoost classifier on the last hidden layer (FFN) features.  

### 3. `evaluate/`
Scripts for **evaluating trained models**:
- `evaluate_GCN.py` → Evaluate Chemprop-trained GCN.  
- `evaluate_last_ffn_SVM.py` → Evaluate SVM performance on test set.  
- `evaluate_last_ffn_XGBoost.py` → Evaluate XGBoost performance on test set.  

### 4. `error_data_detection_and_cleaning/`
Implements the **error detection and data cleaning pipeline** described in the paper:  
- `cross_training.py` → Perform repeated cross-training experiments (`n=3/5/10`) to identify noisy labels.  
- `detection.ipynb` → Use **mis-predicted samples** and search for **similar molecules** in the training set to locate potential mislabels.  
- `error_prediction_and_cleaning.ipynb` → Aggregate repeat-10 results to generate a cleaned dataset.  
- `README.md` → Documentation of methodology.  

### 5. `intepretability/`
Notebooks and scripts for **interpretability analysis**:  
- `murcko_scaffold_analysis.py` → Extract and visualize Murcko scaffolds.  
- `NPClassifier.ipynb` → Compare predictions against **[NPClassifier](https://pubs.acs.org/doi/10.1021/acs.jnatprod.1c00399)** categories.  
- `Monte_Carlo_Tree_Search.ipynb` → Monte Carlo Tree Search–based rationales (Chemprop interpret module).  
- `rationale_visualization.ipynb` → Visualize atom/substructure-level rationales.  
- `README.md` → Documentation of analysis.  

### 6. other scripts
- `preprocessing.ipynb` → Data cleaning and preprocessing pipeline.  
- `final_testset.ipynb` → Preparing test sets for final evaluation.  
- `fingerprints.ipynb` → Fingerprint generation for downstream classifiers.  
- `predict_chemprop.sh` → Run Chemprop prediction with trained models.  
- `finetune.ipynb` → Fine-tune pretrained Chemprop model.  

