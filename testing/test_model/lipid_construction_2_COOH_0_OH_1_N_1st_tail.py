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
    nmt_pred = reaction_predictors.OpenNMTServerPredictor()
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
    head_df = pd.read_csv(input_head_path)
    selected_head = head_df['SMILES'].sample(n, replace=True).tolist()

    tail_df = pd.read_csv(input_tail_path)
    selected_tail = tail_df['Tails'].sample(n, replace=True).tolist()

    list_of_reactants = [0 for _ in range(n)]
    for i in tqdm(range(n)):
        list_of_reactants[i] = multiset.FrozenMultiset([selected_head[i], selected_tail[i]])

    print('Start reaction prediction')
    out, scores = lipid_construction(list_of_reactants)
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
                    row = list(list_of_reactants[i]) + sublist + [this_score]
                    writer.writerow(row)
    return

if __name__ == '__main__':
    input_head_path = '/home/yo279/rds/hpc-work/project/Github_lipid_dataset_creation/Lipid_reaction_dataset_creation/Data/Data01/filtered_2-COOH0-OH1-N.csv'
    input_tail_path = '/home/yo279/rds/hpc-work/project/Github_lipid_dataset_creation/Lipid_reaction_dataset_creation/Data/Data01/purchasable_tails.csv'
    output_lipid_path = '/home/yo279/rds/hpc-work/project/Github_lipid_dataset_creation/Lipid_reaction_dataset_creation/SyntheticData/filtered_2-COOH0-OH1-N_one_tail_purchasable.csv'
    n = 20000
    main(input_head_path, input_tail_path, output_lipid_path, n)
