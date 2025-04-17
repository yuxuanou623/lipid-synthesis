import csv
from rdkit import Chem
import sys
from os import path

# Import SA scorer
sys.path.append(path.join('/scratch/yo279/miniconda3/envs/dogae_py3.7_pt1.4/share/RDKit/Contrib/', 'SA_Score'))
import sascorer

def calculate_sa_score(smiles):
    """Calculate the SA score for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # Invalid SMILES
    try:
        sa_score = sascorer.calculateScore(mol)  # Negative for standard convention
        return sa_score
    except Exception as e:
        print(f"Error calculating SA score for SMILES '{smiles}': {e}")
        return None

def process_smiles_file(file_path):
    """Read SMILES from a file, check validity, and calculate SA scores."""
    results = []
    with open(file_path, 'r') as file:
        smiles_list = file.readlines()
    
    for smiles in smiles_list:
        smiles = smiles.strip()  # Remove any leading/trailing whitespace
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            sa_score = calculate_sa_score(smiles)
            results.append((smiles, True, sa_score))  # Valid SMILES
        else:
            results.append((smiles, False, None))  # Invalid SMILES
    
    return results

def save_results_to_csv(results, output_file):
    """Save results to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write header
        csvwriter.writerow(['SMILES', 'Validity', 'SA Score'])
        # Write rows
        for smiles, is_valid, sa_score in results:
            csvwriter.writerow([smiles, is_valid, sa_score if sa_score is not None else 'N/A'])

# File path to the .txt document containing SMILES
input_file_path = '/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dogae/prior_sampling/published.txt'  # Replace with your actual file path
output_file_path = '/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dogae/prior_sampling/published_results.csv'  # Replace with your desired output file path

# Process the file and save results to CSV
results = process_smiles_file(input_file_path)
save_results_to_csv(results, output_file_path)

print(f"Results saved to {output_file_path}")
