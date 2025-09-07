import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem

def extract_properties_from_sdf(sdf_file, species_name):
    supplier = Chem.SDMolSupplier(sdf_file)
    data = {
        'Species': [],
        'Molecular Weight': [],
        'AlogP': [],
        'H-bond Donors': [],
        'H-bond Acceptors': [],
        'Rotatable Bonds': [],
        'Polar Surface Area': []
    }
    
    for mol in supplier:
        if mol is not None:
            data['Species'].append(species_name)
            data['Molecular Weight'].append(float(mol.GetProp('MOLECULAR_WEIGHT')))
            data['AlogP'].append(float(mol.GetProp('ALOGP')))
            data['H-bond Donors'].append(int(mol.GetProp('HBD')))
            data['H-bond Acceptors'].append(int(mol.GetProp('HBA')))
            data['Rotatable Bonds'].append(int(mol.GetProp('ROTATABLE_BONDS')))
            data['Polar Surface Area'].append(float(mol.GetProp('POLAR_SURFACE_AREA')))
    
    return pd.DataFrame(data)

df_animalia = extract_properties_from_sdf("animalia.sdf", "Animalia")
df_bacteria = extract_properties_from_sdf("bacteria.sdf", "Bacteria")
df_fungi = extract_properties_from_sdf("fungi.sdf", "Fungi")

df = pd.concat([df_animalia, df_bacteria, df_fungi])

fig, axes = plt.subplots(2, 3, figsize=(18, 12)) 

properties = ['Molecular Weight', 'AlogP', 'H-bond Donors', 'H-bond Acceptors', 'Rotatable Bonds', 'Polar Surface Area']
for i, prop in enumerate(properties):
    ax = axes[i//3, i%3]  
    sns.violinplot(x='Species', y=prop, data=df, palette='Set2', ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(f'{prop}')

plt.tight_layout()
plt.savefig('../../results/figures/properties_taxonomy.svg',format='svg')