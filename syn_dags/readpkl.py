import pickle


# Path to your pickle file
filename = '/scratch/yo279/Lipid_reaction_dataset_creation/synthesis-dags/scripts/dog_gen/out_samples/doggen_sampling_on_weights_epoch-30_time-24-06-13_20:10:41.pth_run_at_24-06-13_21:10:22.pick'  # Ensure the path and the extension are correct
def nx_to_tuple_tree(nx_tree, root_smi):
        def recursive_build(child_node):
            return [(elem[0][0], recursive_build(elem[0])) for elem in nx_tree.in_edges(child_node)]
        return (root_smi, recursive_build((root_smi,)))
# Open the file in binary read mode
with open(filename, 'rb') as file:
    # Load data from the file
    data = pickle.load(file)

# Now you can use 'data' as a normal Python object
all_syn_trees = data['all_syn_trees']
for syn_tree in all_syn_trees:  # Adjust the slice as needed
    print(syn_tree)
    print("Root SMILES:", syn_tree.root_smi)
    print("Number of Nodes:", syn_tree.num_nodes)
    print("Unique Molecules:", syn_tree.unique_molecules_in_tree)
    print("Order for Construction:", syn_tree.order_for_construction)
    tuple_tree = nx_to_tuple_tree(syn_tree.tree, syn_tree.root_smi)
    print("tuple_tree")
    print(tuple_tree)
    print("-----")
