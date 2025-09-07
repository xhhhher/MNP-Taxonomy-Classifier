import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from matplotlib_venn import venn3
import os

sns.set_palette("Set2")

def read_sdf_with_origin(filename):
    suppl = Chem.SDMolSupplier(filename)
    origin = os.path.basename(filename).split('.')[0]  
    for mol in suppl:
        if mol is not None:
            mol.SetProp('origin', origin)  
            yield mol  

def get_scaffold(mol):
    return MurckoScaffold.GetScaffoldForMol(mol)

def calculate_scaffolds_by_label(mols, origin):
    scaffolds = set()
    for mol in mols:
        if mol.GetProp("origin") == origin:  
            scaffold = Chem.MolToSmiles(get_scaffold(mol))
            if scaffold:  
                scaffolds.add(scaffold)
    return scaffolds

def generate_fingerprints(mols):
    fingerprints = []
    for mol in mols:
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 6, nBits=1024)
        fingerprints.append(fingerprint)
    return fingerprints

def pca_and_clustering(fingerprints, n_clusters=5000):
    pca = PCA(n_components=85) 
    pca_result = pca.fit_transform(fingerprints)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    clusters = kmeans.fit_predict(pca_result)
    return clusters

def main():
    fig, axes = plt.subplots(1, 2, figsize=(18, 10))
    df_file_animalia = 'animalia.sdf'
    sdf_file_bacteria = 'bacteria.sdf'
    sdf_file_fungi = 'fungi.sdf'
    
    mols_animalia = list(read_sdf_with_origin(sdf_file_animalia))
    mols_bacteria = list(read_sdf_with_origin(sdf_file_bacteria))
    mols_fungi = list(read_sdf_with_origin(sdf_file_fungi))
    
    scaffolds_animalia = calculate_scaffolds_by_label(mols_animalia, "animalia")
    scaffolds_bacteria = calculate_scaffolds_by_label(mols_bacteria, "bacteria")
    scaffolds_fungi = calculate_scaffolds_by_label(mols_fungi, "fungi")
    
    colors = sns.color_palette("Set2", n_colors=3)  

    venn3(subsets=(len(scaffolds_animalia - scaffolds_bacteria - scaffolds_fungi),
                   len(scaffolds_bacteria - scaffolds_animalia - scaffolds_fungi),
                   len(scaffolds_animalia & scaffolds_bacteria - scaffolds_fungi),
                   len(scaffolds_fungi - scaffolds_animalia - scaffolds_bacteria),
                   len(scaffolds_animalia & scaffolds_fungi - scaffolds_bacteria),
                   len(scaffolds_bacteria & scaffolds_fungi - scaffolds_animalia),
                   len(scaffolds_animalia & scaffolds_bacteria & scaffolds_fungi)),
          set_labels=('Animalia', 'Bacteria', 'Fungi'),
          ax=axes[0], 
          set_colors=colors, 
          alpha=0.5)
    axes[0].set_title('Overlap of Murcko Scaffolds')

    mols = mols_animalia + mols_bacteria + mols_fungi 
    fingerprints = generate_fingerprints(mols)
    clusters = pca_and_clustering(fingerprints)
    animalia_indices = [i for i, mol in enumerate(mols) if mol.GetProp("origin") == "animalia"]
    bacteria_indices = [i for i, mol in enumerate(mols) if mol.GetProp("origin") == "bacteria"]
    fungi_indices = [i for i, mol in enumerate(mols) if mol.GetProp("origin") == "fungi"]
    
    clusters_animalia = clusters[animalia_indices]
    clusters_bacteria = clusters[bacteria_indices]
    clusters_fungi = clusters[fungi_indices]
    
    venn3(subsets=(len(set(clusters_animalia) - set(clusters_bacteria) - set(clusters_fungi)),
                   len(set(clusters_bacteria) - set(clusters_animalia) - set(clusters_fungi)),
                   len(set(clusters_animalia) & set(clusters_bacteria) - set(clusters_fungi)),
                   len(set(clusters_fungi) - set(clusters_animalia) - set(clusters_bacteria)),
                   len(set(clusters_animalia) & set(clusters_fungi) - set(clusters_bacteria)),
                   len(set(clusters_bacteria) & set(clusters_fungi) - set(clusters_animalia)),
                   len(set(clusters_animalia) & set(clusters_bacteria) & set(clusters_fungi))),
          set_labels=('Animalia', 'Bacteria', 'Fungi'),
          ax=axes[1], 
          set_colors=colors,
          alpha=0.5)
    axes[1].set_title('Overlap of Morgan Fingerprint Clusters')
    
    plt.tight_layout()
    plt.savefig('../../results/figures/scaffold_venn.svg',format='svg')

if __name__ == "__main__":
    main()
 