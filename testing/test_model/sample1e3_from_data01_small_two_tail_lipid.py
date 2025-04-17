import requests
import pytest
import multiset
import csv
import pandas as pd
from syn_dags.model import reaction_predictors
from multiset import Multiset
from tqdm import tqdm
# this is to evaulate the lipid rate of randomly select three tails, whether the product is lipid.
def is_multiset(element):
    # Check if the element is an instance of Multiset
    return isinstance(element, Multiset)

def lipid_construction(list_of_reactants):
    nmt_pred = reaction_predictors.OpenNMTServerPredictor(server_address = 'http://localhost:6000/translator/translate')
    print("After nmt_pred = reaction_predictors.OpenNMTServerPredictor()")

    all_out = []
    all_scores = []
    batch_size = 100
    # Process in batches
    for start_idx in tqdm(range(0, len(list_of_reactants), batch_size)):
        end_idx = min(start_idx + batch_size, len(list_of_reactants))
        batch = list_of_reactants[start_idx:end_idx]
        try:
            out= nmt_pred._run_list_of_reactant_sets(batch)
            all_out.extend(out)
            
        except requests.exceptions.ConnectionError as e:
            print(f"Connection error for batch {start_idx//batch_size}: {str(e)}")
            # Retry logic or handling code can go here

    return all_out

def main(input_head_path, input_tail_path, output_lipid_path, n):
    # Read head and tail data
    head_df = pd.read_csv(input_head_path)
    tail_df = pd.read_csv(input_tail_path)

    # Sample entire rows instead of just a column
    selected_head = head_df
    selected_tail = tail_df.sample(n, replace=True).reset_index(drop=True)

    list_of_reactants = [0 for _ in range(n)]
    for i in tqdm(range(n)):
        # Storing multiset of selected head and tail
        list_of_reactants[i] = multiset.FrozenMultiset([selected_head.iloc[i]['ip1'], selected_tail.iloc[i]['SMILES']])

    print('Start reaction prediction')
    out = lipid_construction(list_of_reactants)  # This function must be defined elsewhere
    print('Finish reaction prediction')

    # Write successfully constructed product to file
    with open(output_lipid_path, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(out)):
            sublist = out[i]

            # Filter out unsuccessful reaction
            if not is_multiset(sublist) and sublist[0] not in list(list_of_reactants[i]):
    
                # Include all original data from selected_head and selected_tail
                head_data = selected_head.iloc[i].tolist()
                tail_data = selected_tail.iloc[i].tolist()
                # Create the row combining all necessary information
                row = head_data + tail_data + sublist 
                writer.writerow(row)

if __name__ == '__main__':
    input_head_path = '/scratch/yo279/Lipid_reaction_dataset_creation/SyntheticData/smiles_index_new_onetail.csv'
    input_tail_path = '/scratch/yo279/Lipid_reaction_dataset_creation/SyntheticData/smiles_index_new_tail.csv'
    output_lipid_path = '/scratch/yo279/Lipid_reaction_dataset_creation/SyntheticData/smiles_index_new_twotail.csv'
    n = 23094
    main(input_head_path, input_tail_path, output_lipid_path, n)
