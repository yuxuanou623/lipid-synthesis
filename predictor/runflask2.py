from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
from custom2removeunused import CUSTOM_REACTIONS
import random
app = Flask(__name__)



def reaction(head, tail_smiles1):
    all_indices = list(range(14))  # 0 to 13 inclusive
    tried_indices = set()
    smiles_list = []
    template_index = -1

    while len(tried_indices) < len(all_indices):
        reaction_index = random.choice([i for i in all_indices if i not in tried_indices])
        tried_indices.add(reaction_index)

        reaction = CUSTOM_REACTIONS[reaction_index]
        template_index = 0
        smiles_list = []

        # First direction: [head, tail]
        mols1 = [Chem.MolFromSmiles(head), Chem.MolFromSmiles(tail_smiles1)]
        products = reaction.run_reactants(mols1)
        for product_tuple in products:
            mol = product_tuple[0]
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
            smiles_list.append(smiles)

        # Second direction: [tail, head] if no product
        if len(smiles_list) == 0:
            template_index = 1
            mols2 = [Chem.MolFromSmiles(tail_smiles1), Chem.MolFromSmiles(head)]
            products = reaction.run_reactants(mols2)
            for product_tuple in products:
                mol = product_tuple[0]
                smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
                smiles_list.append(smiles)

        # If exactly one product found, return it
        if len(smiles_list) !=0:
            return list(set(smiles_list))[0], template_index, reaction_index

    # If no reaction yields exactly one product, return all unique products from the last attempt
    return head, template_index, reaction_index

@app.route('/reaction', methods=['POST'])
def run_reaction():
    data = request.json
    if not isinstance(data, list):
        return jsonify({"error": "Input must be a list of dictionaries."}), 400

    results = []

    for i, item in enumerate(data):
      
        head = item.get("head")
        tail_smiles1 = item.get("tail_smiles1")
  

        if not all([head, tail_smiles1]):
            print("i, item",i, item)
            results.append({
                "index": i,
                "error": "Missing one or more required fields.",
            })
            continue

        try:
            products, template_index, reaction_index = reaction(head, tail_smiles1)
            results.append({
                "index": i,
                "products": products,
                "head":head,
                "tail": tail_smiles1,
                "reaction_fg1":reaction_index,
                "template_index":template_index, 
            })
        except Exception as e:
            results.append({
                "index": i,
                "error": str(e)
            })
    return jsonify(results)

if __name__ == '__main__':
    app.run(debug=False, port=8080)
