"""
Script to calculate the extinction coefficient at 214 nm for a given polypeptide sequence.
Associated paper: https://pubs.acs.org/doi/10.1021/jf070337l

The authors calculate the 214 nm extinction coefficient of peptides/proteins by 
summing the extinction coefficients of peptide bonds (923 M-1 cm-1) and free amino acids, 
with an additional 2675 M-1 cm-1 accounted for Proline when not at the N-terminus.

The function mimics equation 2 in the associated paper.

The script takes a sequence as input and calculates the 214 nm extinction coefficient. 
It includes an option for specifying whether the sequence has an acetyl N-terminal cap, 
which influences the calculation.

The script requires Python 3.x and has been tested with Python 3.9.
Dependencies: argparse, collections.Counter, colorama

Usage: 
python calc_e214.py 'GKIAKSLKKIAQSLKKIAKSLFG' --acetyl

Author: Joel J. Chubb
Date: 08.07.2023
"""

import argparse
from collections import Counter

# Set up argument parsing
parser = argparse.ArgumentParser(description='Calculate the extinction coefficient at 214 nm of a given sequence.')
parser.add_argument('sequence', type=str, help='The sequence for which the coefficient is to be calculated.')
parser.add_argument('--acetyl', dest='acetyl', action='store_true', help='Specify this flag if the sequence has an acetyl N-terminal cap.')
parser.set_defaults(acetyl=False)

args = parser.parse_args()

# Define a dictionary of amino acid absorption coefficients
Abs = {
    "A":32,
    "R":102,
    "N":136,
    "D":58,
    "C":225,
    "Q":142,
    "E":78,
    "G":21,
    "H":5125,
    "I":45,
    "L":45,
    "K":41,
    "M":980,
    "F":5200,
    "P":30,
    "S":34,
    "T":41,
    "W":29050,
    "Y":5375,
    "V":43
}

# Create a counter of the amino acids in the sequence
counter = Counter(args.sequence)

# Calculate the free amino acid absorption
free_AA_abs = sum([Abs[char]*count for char, count in counter.items()])

# Calculate the number of peptide bonds
num_bonds = len(args.sequence)

# Define the absorption of a peptide bond at 214 nm
peptide_bond_214 = 923

# Calculate the extinction coefficient at 214 nm based on whether the sequence has an acetyl cap and starts with 'P'
if not args.acetyl and args.sequence[0] == 'P': 
    ex_coeff214 = (peptide_bond_214 * num_bonds) + free_AA_abs - peptide_bond_214 + 2645
elif not args.acetyl: 
    ex_coeff214 = (peptide_bond_214 * num_bonds) + free_AA_abs - peptide_bond_214
else: 
    ex_coeff214 = (peptide_bond_214 * num_bonds) + free_AA_abs

# Print the results with colorama
print("\n############### RESULTS ###############")
print(f"Sequence: {args.sequence}")
print(f"Acetyl N-cap: {args.acetyl}")
print(f"\u03B5 (214 nm): {ex_coeff214} M-1 cm-1")
print("#######################################\n")