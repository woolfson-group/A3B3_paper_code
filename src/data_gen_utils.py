
# internal modules
import subprocess

# external modules
import ampal
import numpy as np

def calculate_z_shift_between_atoms(pdb_ampal: ampal.assembly.Assembly, atom1_coords: np.ndarray, atom2_coords: np.ndarray) -> float:
    """
    Calculate the z-shift between two atoms in a PDB file.
    
    Parameters:
    - pdb_ampal (ampal.assembly.Assembly): The PDB file as an ampal assaembly.
    - atom1_coords (np.ndarray): The coordinates of the first atom.
    - atom2_coords (np.ndarray): The coordinates of the second atom.
    
    Returns:
    - float: The z-shift, which is the Euclidean distance between the foot coordinates of the two atoms.
    """
    
    # Calculate the reference axis from the chains in the PDB file
    ref_ax = ampal.analyse_protein.reference_axis_from_chains(pdb_ampal)
    
    # Get the start and end points of superhelical axis for vector
    start = ref_ax.coordinates[0]
    end = ref_ax.coordinates[-1]
    
    # Call function to find the foot coordinates for both atoms
    atom1_foot = ampal.analyse_protein.find_foot(start, end, atom1_coords)
    atom2_foot = ampal.analyse_protein.find_foot(start, end, atom2_coords)
    
    # Calculate the z-shift as the Euclidean distance between the foot coordinates of the two atoms
    z_shift = np.linalg.norm(atom1_foot - atom2_foot)
    
    # Return the calculated z-shift
    return z_shift

def p_or_ap(ampal_pdb: ampal.Assembly) -> str:
    """
    Determines if the chains in a protein structure are all parallel ('p') or anti-parallel ('ap') 
    relative to a reference chain (chain A).

    The function compares the vectors formed by the alpha carbon (CA) atoms of the first 
    and last residues of each chain in the structure. It calculates the dot product between 
    the reference chain (chain A) and every other chain to assess the relative orientation.
    
    Args:
        ampal_pdb (ampal.Assembly): An AMPAL object representing a protein structure 
                                    containing multiple chains.

    Returns:
        str: 'p' if all chains are parallel to the reference chain ('A'), 'ap' if at least 
                one chain is anti-parallel to the reference chain.
    """
    
    # Extract the coordinates of the alpha carbon (CA) atoms from the first chain (chain A)
    first_chain_coords = [residue['CA'].array for residue in ampal_pdb['A']]

    # Initialize a dictionary to store CA coordinates of all chains
    all_chain_coords = {}
    for chain in ampal_pdb:
        # Store the CA coordinates for each chain
        all_chain_coords[chain.id] = [residue['CA'].array for residue in chain]

    # Calculate the vector between the first and last residues of chain A
    first_chain_vector = first_chain_coords[0] - first_chain_coords[-1]

    # Initialize a list to store the orientation results (1 for parallel, 0 for anti-parallel)
    oris = []
    
    # Loop over all chains to compare their vectors with the reference chain (chain A)
    for chain, coords in all_chain_coords.items():
        # Skip the reference chain (chain A)
        if chain != "A":
            # Calculate the vector between the first and last residues of the current chain
            chain_vector = coords[0] - coords[-1]
            # Compute the dot product between the reference chain and the current chain
            dot_prod = np.dot(first_chain_vector, chain_vector)
            # If the dot product is positive, chains are parallel; otherwise, anti-parallel
            if dot_prod > 0:
                oris.append(1)
            else:
                oris.append(0)
    
    # If all chains are parallel to the reference chain, return 'p'; otherwise, return 'ap'
    if sum(oris) == len(ampal_pdb) - 1:
        return 'p'
    else:
        return 'ap'
    
# Function to run USalign with specified parameters
def run_usalign(p1, p2, path_to_USalign='../../../repos/USalign/USalign'):
    """
    Executes the USalign command-line tool to align two protein structures and capture the output.
    
    Args:
        p1 (str): Path to the first protein structure file.
        p2 (str): Path to the second protein structure file.
        path_to_USalign (str): Path to the US align executable
        
    Returns:
        str: The alignment result, specifically the second line of the USalign output.
    """
    # Run the USalign command with specified options and capture the output
    process = subprocess.run(
        [path_to_USalign, p1, p2, '-mm', '1', '-ter', '0', '-outfmt', '2'],
        capture_output=True, check=True
    )
    # Decode the output, split by newline, and return the second line
    return process.stdout.decode().split('\n')[1]
