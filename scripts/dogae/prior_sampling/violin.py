import pickle
import matplotlib.pyplot as plt
import seaborn as sns

# Specify the path to your pickle file
pickle_file_path = '/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dogae/prior_sampling/moleculartransformer20epoch_lipid_with_ph.pkl'

# Function to load data from a pickle file
def load_data_from_pickle(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data

# Load data from the pickle file
data = load_data_from_pickle(pickle_file_path)

# Check if the data is a list and contains numerical data
if isinstance(data, list) and all(isinstance(x, (int, float)) for x in data):
    # Create a violin plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=data)
    plt.title('Violin Plot of Nearest Neighbor Distances')
    plt.xlabel('Reaction Predictor: Molecular Transformer, Synthesis DAG')
    plt.ylabel('Density')
    plt.ylim(0, 1.2)  # Fixing the y-axis between 0 to 1.2
    
    # Save the figure
    plt.savefig('/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dogae/prior_sampling/output_MT_20epoch.png', dpi=300, bbox_inches='tight')  # Specify your desired path and file name here
    plt.show()
else:
    print("Data is not a list of numerical values. Please check the data format.")
