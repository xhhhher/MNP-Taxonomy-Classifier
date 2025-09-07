# 🧩 Interpretability Analysis  

This directory contains scripts and notebooks used for the **interpretability analysis** described in our paper.  
The goal of this part is to understand **which chemical features or scaffolds drive the classification decisions** of our natural product origin prediction models.  

## 📂 Contents

- **`murcko_scaffold_analysis.py`**  
  Extracts Murcko scaffolds from molecules and visualizes the most frequent scaffolds per origin (Animalia, Bacteria, Fungi).  
  This corresponds to the scaffold-level analysis in the paper, where we compared the chemical backbones enriched in different biological origins.  

- **`NPClassifier.ipynb`**  
  Notebook for analyzing predictions with respect to existing classification systems (e.g., **[NPClassifier](https://pubs.acs.org/doi/10.1021/acs.jnatprod.1c00399)** categories).  
  This allows comparison between model-derived rationales and curated ontology labels.  

- **`Monte_Carlo_Tree_Search.ipynb`**  
  Implements **Monte Carlo Tree Search (MCTS)**-based interpretability.  
  - MCTS is a heuristic search algorithm widely used in decision-making problems.  
  - In our framework, MCTS was adapted to explore **substructures (rationales)** that contribute most to classification.  
  - Starting from the full molecular graph, MCTS iteratively samples and expands subgraphs, retaining those that maximize model confidence for a given class.  
  - By repeating this process, the algorithm identifies **minimal predictive fragments**—small sets of atoms/bonds sufficient to drive the same prediction.  
  - This helps highlight **key substructures** that underlie the origin-specific decision boundary.  

- **`rationale_visualization.ipynb`**   
  Highlights which atoms or substructures contribute most strongly to the prediction, providing a molecular-level explanation of model decisions.  

## 📝 Notes  

- Scaffold analysis was carried out using **RDKit**’s `MurckoScaffold` implementation.  
- Monte Carlo Tree Search leverages Chemprop’s **interpretability module** (currently available only in version 1), consistent with our pipeline in the main manuscript.  
- These analyses complement the **predictive performance results** by offering **chemical insights** into how the model distinguishes origins.  