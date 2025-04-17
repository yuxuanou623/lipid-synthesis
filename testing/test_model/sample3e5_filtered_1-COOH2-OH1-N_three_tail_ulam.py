import requests
import pytest
import multiset
import csv
import pandas as pd
from syn_dags.model import reaction_predictors
from multiset import Multiset
from tqdm import tqdm

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
            out, scores = nmt_pred._run_list_of_reactant_sets(batch)
            all_out.extend(out)
            all_scores.extend(scores)
        except requests.exceptions.ConnectionError as e:
            print(f"Connection error for batch {start_idx//batch_size}: {str(e)}")
            # Retry logic or handling code can go here

    return all_out, all_scores

def main(input_head_path, input_tail_path, output_lipid_path, n):
    # Read head and tail data
    head_df = pd.read_csv(input_head_path)
    tail_df = pd.read_csv(input_tail_path)

    # Sample entire rows instead of just a column
    selected_head = head_df.sample(n, replace=True).reset_index(drop=True)
    selected_tail = tail_df.sample(n, replace=True).reset_index(drop=True)

    list_of_reactants = [0 for _ in range(n)]
    for i in tqdm(range(n)):
        # Storing multiset of selected head and tail
        list_of_reactants[i] = multiset.FrozenMultiset([selected_head.iloc[i]['6'], selected_tail.iloc[i]['Tails']])

    print('Start reaction prediction')
    out, scores = lipid_construction(list_of_reactants)  # This function must be defined elsewhere
    print('Finish reaction prediction')

    # Write successfully constructed product to file
    with open(output_lipid_path, 'w', newline='') as file:
        writer = csv.writer(file)
        for i in range(len(out)):
            sublist = out[i]
            this_score = scores[i]
            # Filter out unsuccessful reaction
            if not is_multiset(sublist) and sublist[0] not in list(list_of_reactants[i]):
                if this_score > -0.2:
                    # Include all original data from selected_head and selected_tail
                    head_data = selected_head.iloc[i].tolist()
                    tail_data = selected_tail.iloc[i].tolist()
                    # Create the row combining all necessary information
                    row = head_data + tail_data + sublist + [this_score]
                    writer.writerow(row)

if __name__ == '__main__':
    input_head_path = '/scratch/yo279/Lipid_reaction_dataset_creation/SyntheticData/twotaillipid/sample3e5_1-COOH2-OH1-N_two_tails_lipids.csv'
    input_tail_path = '/scratch/yo279/Lipid_reaction_dataset_creation/Data/Data01/purchasable_tails.csv'
    output_lipid_path = '/scratch/yo279/Lipid_reaction_dataset_creation/SyntheticData/sample3e5_1-COOH2-OH1-N_three_tails_lipid.csv'
    n = 300000
    main(input_head_path, input_tail_path, output_lipid_path, n)
