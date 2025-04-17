from flask import Flask, request, jsonify
from rdkit import Chem
from rdkit.Chem import AllChem
from custom2removeunused import CUSTOM_REACTIONS
head = 'CC#CCCCCCCCCCN'
tail_smiles1 = 'COCCN(CC(=O)O)C[C@H]1C[C@@H](CO)CO1'
reaction_fg1 =  int(str(1))

reaction = CUSTOM_REACTIONS[reaction_fg1]
mols1 = [Chem.MolFromSmiles(head), Chem.MolFromSmiles(tail_smiles1)]
products = reaction.run_reactants(mols1)
print(products)
for product in products:
    print(Chem.MolToSmiles(Chem.RemoveHs(product[0])))