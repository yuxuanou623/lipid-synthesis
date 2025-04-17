import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Function to load data from a pickle file
def load_data_from_pickle(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data

# Load data from specified pickle files
mtpickle_file_path = '/home/yo279/rds/hpc-work/project/Github_lipid_dataset_creation/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dog_gen/round_pickle/moleculartransformer20epoch_lipid_with_ph.pkl'
mtdata = load_data_from_pickle(mtpickle_file_path)

chem_file_path = '/home/yo279/rds/hpc-work/project/Github_lipid_dataset_creation/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dog_gen/round_pickle/output_chemformer_10epoch_lipid_with_ph.pkl'
chemdata = load_data_from_pickle(chem_file_path)

linear_file_path = '/home/yo279/rds/hpc-work/project/Github_lipid_dataset_creation/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dog_gen/round_pickle/output_linear_20epoch_lipid_with_ph.pkl'
lineardata = load_data_from_pickle(linear_file_path)

# Since the datasets might be lists and could be of different lengths, we will use a longer form DataFrame
data = {
    'Model': ['DAG+Chem'] * len(chemdata) + ['DAG+MT'] * len(mtdata) + ['List+Chem'] * len(lineardata),
    'Density': chemdata + mtdata + lineardata
}

df = pd.DataFrame(data)

# Create a violin plot for the combined dataset
plt.figure(figsize=(10, 6))
sns.violinplot(x='Model', y='Density', data=df)
plt.title('Violin Plot of Nearest Neighbor Distances')
plt.xlabel('Models')
plt.ylabel('Density')

# Optionally, adjust y-axis limits if necessary
plt.ylim(0, 1.2)  # Example y-axis limits

# Save and show the plot
plt.savefig('output_3vios.png', dpi=300, bbox_inches='tight')
plt.show()
