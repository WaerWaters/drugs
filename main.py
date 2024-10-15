import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdmolops, AllChem
import plotly.graph_objects as go

# Define the InChIKey of the compound
inchikey = 'MAEIEVLCKWDQJH-UHFFFAOYSA-N'

# Convert InChIKey to SMILES using PubChem
compound = pcp.get_compounds(inchikey, 'inchikey')[0]
smiles = compound.canonical_smiles

def get_molecule_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    smarts = Chem.MolToSmarts(mol)
    print(smiles, smarts)
    updated_mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(updated_mol)
    sssr = rdmolops.GetSymmSSSR(updated_mol)
    print(list(sssr[0]), list(sssr[1]))
    substructure_dict = {}

    # Gather atom information
    atoms = []
    for atom in updated_mol.GetAtoms():
        pos = updated_mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atom_info = {
            'index': atom.GetIdx(),
            'element': atom.GetSymbol(),
            'degree': atom.GetDegree(),
            'is_aromatic': atom.GetIsAromatic(),
            'x': pos.x,
            'y': pos.y,
        }
        atoms.append(atom_info)

    # Gather bond information
    bonds = []
    for bond in updated_mol.GetBonds():
        bond_info = {
            'start_atom_idx': bond.GetBeginAtomIdx(),
            'end_atom_idx': bond.GetEndAtomIdx(),
            'bond_type': str(bond.GetBondType()),
            'is_aromatic': bond.GetIsAromatic()
        }
        bonds.append(bond_info)

    # Print extracted atom and bond information
    print("\nAtoms:")
    for atom in atoms:
        print(atom)

    print("\nBonds:")
    for bond in bonds:
        print(bond)
    
    return atoms, bonds




def visualize_2d_molecule(atoms, bonds):
    # Visualization using Plotly
    fig = go.Figure()

    # Draw the bonds
    for bond in bonds:
        start_atom = atoms[bond['start_atom_idx']]
        end_atom = atoms[bond['end_atom_idx']]
        bond_type = bond['bond_type']
        hovertext = f"Bond between {start_atom['element']} (index {start_atom['index']}) and {end_atom['element']} (index {end_atom['index']})<br>Type: {bond_type}"

        fig.add_trace(go.Scatter(
            x=[start_atom['x'], end_atom['x']],
            y=[start_atom['y'], end_atom['y']],
            mode='lines',
            line=dict(color='black'),
            hoverinfo='text',
            hovertext=hovertext
        ))

    # Draw the atoms
    x_coords = [atom['x'] for atom in atoms]
    y_coords = [atom['y'] for atom in atoms]
    indices = [atom['index'] for atom in atoms]
    hovertexts = [f"Element: {atom['element']}<br>Index: {atom['index']}<br>Degree: {atom['degree']}<br>Aromatic: {atom['is_aromatic']}" for atom in atoms]

    fig.add_trace(go.Scatter(
        x=x_coords,
        y=y_coords,
        mode='markers+text',
        marker=dict(color='lightgreen', size=20),
        text=indices,
        textposition='middle center',
        hoverinfo='text',
        hovertext=hovertexts
    ))

    # Add labels and title
    fig.update_layout(
        title='2D Visualization of Molecule',
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        showlegend=False
    )

    # Show the plot
    fig.show()
    
#visualize_2d_molecule(get_molecule_properties(smiles)[0], get_molecule_properties(smiles)[1])

def pattern_search():
    smiles = 'CCCCNC1=C(C(=CC(=C1)C(=O)O)S(=O)(=O)N)OC2=CC=CC=C2'
    smarts = '[#6]-[#6]-[#6]-[#6]-[#7]-[#6]1:[#6](:[#6](:[#6]:[#6](:[#6]:1)-[#6](=[#8])-[#8])-[#16](=[#8])(=[#8])-[#7])-[#8]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'
    smarts_sub = '[#6]-[#6]-[#7]-[#6]-[#8]'

    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmiles(smarts_sub)
    hit_ats = list(mol.GetSubstructMatch(pattern))
    hit_bonds = []
    for bond in pattern.GetBonds():
        aid1 = hit_ats[bond.GetBeginAtomIdx()]
        aid2 = hit_ats[bond.GetEndAtomIdx()]
        hit_bonds.append(mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())
    d = Draw.rdMolDraw2D.MolDraw2DCairo(500, 500)

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

""" # This is the code to generate the 2D structure of the molecule
print(Chem.MolToMolBlock(mol))
"""
""" # This adds hydrogen atoms to the molecule
updated_mol = Chem.AddHs(mol)
print(Chem.MolToMolBlock(updated_mol))
"""
""" # this saves the molecule structure to a file
updated_mol = Chem.AddHs(mol)
print(Chem.MolToMolBlock(updated_mol),file=open('foo.mol','w+'))
"""
""" # Get the neighbor atoms of a specific atom
atom = mol.GetAtomWithIdx(22)
print([x.GetAtomicNum() for x in atom.GetNeighbors()])
"""
""" # This gets the rings in the molecule
ssr = Chem.GetSymmSSSR(mol)
print(len(ssr))
print(list(ssr[0]))
"""
""" # This provides an image of the molecule
img = Draw.MolToImage(mol)
img.save('mol.png')
"""
"""
compound_series = c.to_series(properties=['atoms', 'bonds', 'inchi', 'molecular_formula'])
"print(compound_series)"
"""

# An obvious use is to show atoms and bonds that have matched a substructure query
# Substructure Searching
# Generating Similarity Maps Using Fingerprints


