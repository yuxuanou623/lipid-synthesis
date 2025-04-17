import pickle

# Define the path to your pickle file
file_path = '/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/syn_dags/script_utils/data_2_small_numbers_list.pkl'

# Open the file in binary read mode
with open(file_path, 'rb') as file:
    # Load the data from the file
    numbers_list = pickle.load(file)

# Print the data to confirm it has been loaded correctly
print(numbers_list)
print(len(numbers_list))