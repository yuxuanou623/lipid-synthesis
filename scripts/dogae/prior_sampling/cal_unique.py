def count_unique_lipids(file_path):
    # Create a set to store unique lipid names
    unique_lipids = set()
    
    # Open the text file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Strip whitespace from the beginning and end of the line
            lipid = line.strip()
            # Add the cleaned lipid name to the set
            unique_lipids.add(lipid)
    
    # The length of the set is the number of unique lipids
    return len(unique_lipids)

# Specify the path to your text file
file_path = '/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dogae/prior_sampling/samples/doggen_sampling_on_weights_epoch-20_time-24-07-19_linear_data2.pth_run_at_24-07-19_19:34:03_lineararch_chemformer_smiles_20epoch.txt'
# Call the function and print the result
num_unique_lipids = count_unique_lipids(file_path)
print("Number of unique lipids:", num_unique_lipids)
