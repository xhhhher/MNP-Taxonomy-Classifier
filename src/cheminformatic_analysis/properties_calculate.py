from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

def calculate_properties(mol):
    properties = {
        'MolecularWeight': Descriptors.MolWt(mol),
        'LogP': Crippen.MolLogP(mol),
        'H-bond Donors': rdMolDescriptors.CalcNumHDonors(mol),
        'H-bond Acceptors': rdMolDescriptors.CalcNumHAcceptors(mol),
        'Rotatable Bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'TPSA': rdMolDescriptors.CalcTPSA(mol)
    }
    return properties

def add_properties_to_sdf(input_sdf, output_sdf):

    supplier = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)

    for mol in supplier:
        if mol is None: 
            continue
        properties = calculate_properties(mol)
        
        for key, value in properties.items():
            mol.SetProp(key, str(value))

        writer.write(mol)
    writer.close()

add_properties_to_sdf("animalia.sdf","canimalia_with_properties.sdf")
add_properties_to_sdf("bacteria.sdf", "bacteria_with_properties.sdf")
add_properties_to_sdf("fungi.sdf", "fungi_with_properties.sdf")
