import pubchempy as pcp
import pprint
import numpy as np
import plotly.graph_objects as go

results = pcp.get_compounds('MAEIEVLCKWDQJH-UHFFFAOYSA-N', 'inchikey')
print(results)

"""
results = c.record
"""

def molecular_visualization():
    # Get the compound by InChIKey
    results = pcp.get_compounds('MJIHNNLFOKEZEW-UHFFFAOYSA-N', 'inchikey')
    c = pcp.Compound.from_cid(3883)

    # Convert the compound to a dictionary
    compound_dict = c.to_dict(properties=['atoms', 'bonds', 'inchi', 'elements', 'molecular_formula', 'complexity', 'synonyms', 'charge', 'molecular_weight', 'canonical_smiles', 'isomeric_smiles', 'iupac_name'])
    
    pprint.pprint(compound_dict)
    
    # Extract 2D coordinates and element symbols if available
    atoms = compound_dict['atoms']
    coords_2d = [(atom['x'], atom['y'], atom['element'], atom['aid']) for atom in atoms if 'x' in atom and 'y' in atom and 'element' in atom]

    # Extract bond information
    bonds = compound_dict['bonds']

    # Create a dictionary to map atom indices to coordinates
    atom_coords = {atom['aid']: (atom['x'], atom['y']) for atom in atoms}

    # Define a color map for all elements on the periodic table with readable colors
    color_map = {
        'H': 'lightgray', 'He': 'cyan', 'Li': 'purple', 'Be': 'darkgreen', 'B': 'salmon', 'C': 'magenta', 'N': 'blue', 'O': 'red', 'F': 'green', 'Ne': 'orange',
        'Na': 'blue', 'Mg': 'green', 'Al': 'lightgray', 'Si': 'brown', 'P': 'orange', 'S': 'yellow', 'Cl': 'green', 'Ar': 'cyan',
        'K': 'purple', 'Ca': 'darkgreen', 'Sc': 'lightgray', 'Ti': 'lightgray', 'V': 'lightgray', 'Cr': 'lightgray', 'Mn': 'lightgray', 'Fe': 'brown', 'Co': 'lightgray', 'Ni': 'lightgray', 'Cu': 'orange', 'Zn': 'lightgray',
        'Ga': 'lightgray', 'Ge': 'lightgray', 'As': 'lightgray', 'Se': 'lightgray', 'Br': 'brown', 'Kr': 'cyan',
        'Rb': 'purple', 'Sr': 'darkgreen', 'Y': 'lightgray', 'Zr': 'lightgray', 'Nb': 'lightgray', 'Mo': 'lightgray', 'Tc': 'lightgray', 'Ru': 'lightgray', 'Rh': 'lightgray', 'Pd': 'lightgray', 'Ag': 'lightgray', 'Cd': 'lightgray',
        'In': 'lightgray', 'Sn': 'lightgray', 'Sb': 'lightgray', 'Te': 'lightgray', 'I': 'purple', 'Xe': 'cyan',
        'Cs': 'purple', 'Ba': 'darkgreen', 'La': 'lightgray', 'Ce': 'lightgray', 'Pr': 'lightgray', 'Nd': 'lightgray', 'Pm': 'lightgray', 'Sm': 'lightgray', 'Eu': 'lightgray', 'Gd': 'lightgray', 'Tb': 'lightgray', 'Dy': 'lightgray',
        'Ho': 'lightgray', 'Er': 'lightgray', 'Tm': 'lightgray', 'Yb': 'lightgray', 'Lu': 'lightgray',
        'Hf': 'lightgray', 'Ta': 'lightgray', 'W': 'lightgray', 'Re': 'lightgray', 'Os': 'lightgray', 'Ir': 'lightgray', 'Pt': 'lightgray', 'Au': 'orange', 'Hg': 'lightgray',
        'Tl': 'lightgray', 'Pb': 'lightgray', 'Bi': 'lightgray', 'Po': 'lightgray', 'At': 'lightgray', 'Rn': 'cyan',
        'Fr': 'purple', 'Ra': 'darkgreen', 'Ac': 'lightgray', 'Th': 'lightgray', 'Pa': 'lightgray', 'U': 'lightgray', 'Np': 'lightgray', 'Pu': 'lightgray', 'Am': 'lightgray', 'Cm': 'lightgray', 'Bk': 'lightgray', 'Cf': 'lightgray',
        'Es': 'lightgray', 'Fm': 'lightgray', 'Md': 'lightgray', 'No': 'lightgray', 'Lr': 'lightgray',
        'Rf': 'lightgray', 'Db': 'lightgray', 'Sg': 'lightgray', 'Bh': 'lightgray', 'Hs': 'lightgray', 'Mt': 'lightgray', 'Ds': 'lightgray', 'Rg': 'lightgray', 'Cn': 'lightgray', 'Nh': 'lightgray', 'Fl': 'lightgray', 'Mc': 'lightgray', 'Lv': 'lightgray', 'Ts': 'lightgray', 'Og': 'lightgray'
    }

    # Create a plotly figure
    fig = go.Figure()

    # Draw the bonds first without offsets
    for bond in bonds:
        start_atom = bond['aid1']
        end_atom = bond['aid2']
        x_start, y_start = atom_coords[start_atom]
        x_end, y_end = atom_coords[end_atom]
        bond_order = bond.get('order', 1)  # Default to single bond if order is not specified

        hovertext = f"Bond between AID {start_atom} and AID {end_atom}<br>Order: {bond_order}"

        if bond_order == 1:
            fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='black'), hoverinfo='text', hovertext=hovertext))
        elif bond_order == 2:
            fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='blue', dash='dash'), hoverinfo='text', hovertext=hovertext))
        elif bond_order == 3:
            fig.add_trace(go.Scatter(x=[x_start, x_end], y=[y_start, y_end], mode='lines', line=dict(color='red', dash='dot'), hoverinfo='text', hovertext=hovertext))

    # Plot the 2D coordinates with different colors and larger size
    x_coords, y_coords, elements, aids = zip(*coords_2d)
    colors = [color_map.get(element, 'lightgray') for element in elements]  # Default to 'lightgray' if element not in color_map

    fig.add_trace(go.Scatter(
        x=x_coords,
        y=y_coords,
        mode='markers+text',
        marker=dict(color=colors, size=20),
        text=elements,
        textposition='middle center',
        hoverinfo='text',
        hovertext=[f"Element: {element}<br>AID: {aid}<br>X: {x}<br>Y: {y}" for x, y, element, aid in coords_2d]
    ))

    # Add labels and title
    fig.update_layout(
        title='2D Coordinates of Compound',
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        showlegend=False
    )

    # Show the plot
    fig.show()

molecular_visualization()