import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
from collections import Counter
import os

def get_murcko_scaffold(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold, isomericSmiles=False)
    except:
        return None

def get_top_scaffolds(df, label_col="labels", smiles_col="SMILES", top_k=20):
    df["scaffold"] = df[smiles_col].apply(get_murcko_scaffold)
    df = df[df["scaffold"].notnull() & (df["scaffold"] != "")]
    label2scaffold = {}
    for label in df[label_col].unique():
        sub = df[df[label_col] == label]
        scaffold_counts = Counter(sub["scaffold"].dropna())
        top_scaffolds = scaffold_counts.most_common(top_k)
        label2scaffold[label] = top_scaffolds
    return label2scaffold

def draw_scaffold_grid(scaffold_list, label_name, save_path, mols_per_row=5):
    fig_rows = (len(scaffold_list) + mols_per_row - 1) // mols_per_row
    fig, axes = plt.subplots(fig_rows, mols_per_row, figsize=(12, 8))
    axes = axes.flatten()

    draw_options = Draw.rdMolDraw2D.MolDrawOptions()
    draw_options.fixedBondLength = 25

    for i, (scaffold, count) in enumerate(scaffold_list):
        mol = Chem.MolFromSmiles(scaffold)
        if mol:
            drawer = Draw.MolDraw2DCairo(300, 150)
            drawer.SetDrawOptions(draw_options)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            img = drawer.GetDrawingText()

            axes[i].imshow(plt.imread(io.BytesIO(img), format='png'))
            axes[i].set_title(f"{count}", fontsize=11, pad=2)
            axes[i].axis('off')

    for j in range(len(scaffold_list), len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    filename = os.path.join(save_path, f"scaffold_top20_label_{label_name}.svg")
    plt.savefig(filename, format='svg', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {filename}")

def main():
    input_csv = "../../data/processed/trainingset.csv"
    output_dir = "../../results/figures"

    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_csv)
    top_scaffold_dict = get_top_scaffolds(df)

    for label, scaffold_list in top_scaffold_dict.items():
        draw_scaffold_grid(scaffold_list, label_name=label, save_path=output_dir)

if __name__ == "__main__":
    main()

