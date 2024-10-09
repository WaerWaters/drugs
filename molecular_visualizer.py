import pubchempy as pcp
import pprint
import matplotlib.pyplot as plt
import numpy as np

results = pcp.get_compounds('MAEIEVLCKWDQJH-UHFFFAOYSA-N', 'inchikey')
print(results)
"""
results = c.record
"""

def molecular_visualization():
    c = pcp.Compound.from_cid(962)

    # Convert the compound to a dictionary
    compound_dict = c.to_dict(properties=['atoms', 'bonds', 'inchi', 'elements', 'molecular_formula', 'complexity'])

    # Extract 2D coordinates
    atoms = compound_dict['atoms']
    coords_2d = [(atom['x'], atom['y'], atom['element']) for atom in atoms if 'x' in atom and 'y' in atom and 'element' in atom]

    # Extract bond information
    bonds = compound_dict['bonds']

    # Create a dictionary to map atom indices to coordinates
    atom_coords = {atom['aid']: (atom['x'], atom['y']) for atom in atoms}

    # Print the 2D coordinates and atomic letter
    pprint.pprint(coords_2d)

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

    # Draw bonds
    for bond in bonds:
        start_atom = bond['aid1']
        end_atom = bond['aid2']
        x_start, y_start = atom_coords[start_atom]
        x_end, y_end = atom_coords[end_atom]
        bond_order = bond.get('order', 1)  # Default to single bond if order is not specified

        # Calculate the perpendicular offset for double and triple bonds
        dx = x_end - x_start
        dy = y_end - y_start
        length = np.sqrt(dx**2 + dy**2)
        offset_x = -dy / length * 0.05
        offset_y = dx / length * 0.05

        if bond_order == 1:
            plt.plot([x_start, x_end], [y_start, y_end], 'k-', zorder=1)
        elif bond_order == 2:
            plt.plot([x_start, x_end], [y_start, y_end], 'k-', zorder=1)
            plt.plot([x_start + offset_x, x_end + offset_x], [y_start + offset_y, y_end + offset_y], 'k-', zorder=1)
        elif bond_order == 3:
            plt.plot([x_start, x_end], [y_start, y_end], 'k-', zorder=1)
            plt.plot([x_start + offset_x, x_end + offset_x], [y_start + offset_y, y_end + offset_y], 'k-', zorder=1)
            plt.plot([x_start - offset_x, x_end - offset_x], [y_start - offset_y, y_end - offset_y], 'k-', zorder=1)

    # Plot the 2D coordinates with different colors and larger size
    x_coords, y_coords, elements = zip(*coords_2d)
    colors = [color_map[element] for element in elements]
    plt.scatter(x_coords, y_coords, c=colors, s=200, zorder=2)  # s=100 sets the size of the points

    # Add labels and title
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('2D Coordinates of Compound')

    # Annotate each point with the element symbol
    for (x, y, element) in coords_2d:
        plt.annotate(element, (x, y), textcoords="offset points", xytext=(0, 0), ha='center', va='center', zorder=3)

    # Show the plot
    plt.show()
    
molecular_visualization()