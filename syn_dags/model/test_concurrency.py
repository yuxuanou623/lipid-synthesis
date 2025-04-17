import concurrent.futures
import logging
import json
import requests
import typing
import multiset

import abc
import typing
import re
import collections
import json
import logging

import requests
import multiset

#from ..chem_ops import rdkit_general_ops
#from ..utils import misc
#from ..utils import settings
import os
import requests
import json
import csv
import random
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdchem
def are_molecules_identical(smiles1: str, smiles2: str) -> bool:
    """
    Compares two SMILES strings to determine if they represent identical molecules.

    Args:
    smiles1 (str): The SMILES representation of the first molecule.
    smiles2 (str): The SMILES representation of the second molecule.

    Returns:
    bool: True if the molecules are identical, False otherwise.
    """
    try:
        # Convert SMILES to molecule objects
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        if mol1 is None or mol2 is None:
            raise ValueError("One of the SMILES strings could not be converted to a molecule.")

        # Generate canonical SMILES from molecules for accurate comparison
        can_smiles1 = Chem.MolToSmiles(mol1, canonical=True)
        can_smiles2 = Chem.MolToSmiles(mol2, canonical=True)

        # Compare canonical SMILES
        return can_smiles1 == can_smiles2
    except rdchem.KekulizeException:
        print("Failed to kekulize the molecules.")
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False


class MolTransformerTokenizer:
    # from: https://github.com/pschwllr/MolecularTransformer
    THE_REGEX = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    REGEX = re.compile(THE_REGEX)

    @classmethod
    def to_tokens(cls, smiles_str_in: str) -> str:
        return ' '.join(cls.REGEX.findall(smiles_str_in))

    @classmethod
    def from_tokens(cls, tokenized_str_in: str) -> str:
        return ''.join(tokenized_str_in.split())


class AbstractReactionPredictor(metaclass=abc.ABCMeta):
    def __init__(self, size_of_cache=10000):

        # Will store an LRU cache of previous results
        self.cached_results = collections.OrderedDict()
        self.size_of_cache = size_of_cache

    @abc.abstractmethod
    def _run_list_of_reactant_sets(self, list_of_reactant_sets: typing.List[multiset.FrozenMultiset]) -> typing.List[multiset.Multiset]:
        raise NotImplementedError

    def __call__(self, list_of_reactant_sets: typing.List[multiset.FrozenMultiset]) -> typing.List[multiset.Multiset]:
        """
        :param list_of_reactant_sets: List of sets. Each set contains SMILES strings of the reactants in the reaction
        :return: list of product sets: list of multisets with the associated product SMILES strings in the product
        location.
        """

        # Work out which ones we need to query the reaction predictor on
        set_of_already_have = set(self.cached_results.keys())
        reactant_sets_new_needed = list(set(list_of_reactant_sets) - set_of_already_have)
        if len(reactant_sets_new_needed):
            reactants_products_new = list(zip(reactant_sets_new_needed,
                                              self._run_list_of_reactant_sets(reactant_sets_new_needed)))
        else:
            reactants_products_new = []

        # Work out which ones we can just remove out from the cache
        # NB note we POP them out so that when we reinsert them below we put them in at the end.
        reactants_products_new.extend([(item, self.cached_results.pop(item)) for item in
                                        set(list_of_reactant_sets) & set_of_already_have])

        # We collect both the old and new ones together and form the output list
        reactant_products_map = collections.OrderedDict(reactants_products_new)
        out = [reactant_products_map[item] for item in list_of_reactant_sets]

        # We also add these new results to end of our cache and remove anything at the beginning if necessary so that
        # our memory does not blow up.
        self.cached_results.update(reactant_products_map)
        for _ in range(max(0, len(self.cached_results) - self.size_of_cache)):
            self.cached_results.popitem(last=False)

        return out


class OpenNMTServerPredictor(AbstractReactionPredictor):
    """
    Runs a request against OpenNMT's RESTful API
    """
    logger = logging.getLogger('transformer-react-predict')

    def __init__(self, server_address, model_id: int=0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.server_address = settings.get_config().get('Transformer', 'address')
        self.server_address = server_address
        print(f"OpenNMTServerPredictor: server_address is {self.server_address}")
        self.model_id = model_id

    def _run_list_of_reactant_sets(self, list_of_reactant_sets: typing.List[multiset.FrozenMultiset]) -> typing.List[multiset.Multiset]:
        def input_mapper(reactant_multiset_in):
            smiles_str_in = '.'.join(reactant_multiset_in)
            tokenstr = MolTransformerTokenizer.to_tokens(smiles_str_in)
            out_ = {"src": tokenstr, "id": self.model_id}
            return out_
        json_to_send = json.dumps([input_mapper(elem) for elem in list_of_reactant_sets])

        def request_func():
            print("in reaction predictions")
            print("self.server_address")
            print(self.server_address)
            r = requests.post(self.server_address, data=json_to_send, timeout=180)
            r.raise_for_status()
            return r.json()
        return_list = misc.retry_n_times(request_func, 3, Exception, interval=0.5,
                           on_exception=lambda ex: print(ex))
        # response back should look something like:
        # [[{"n_best":1,"pred_score":-0.0002288818359375,"src":"C [S-] . [Mg+] c 1 c c c ( Cl ) c c 1","tgt":"C S c 1 c c c ( Cl ) c c 1"}]]

        def output_mapper(dict_in):
            pred_score = dict_in['pred_score']
            prediction_tokenized = dict_in['tgt']
            prediction = MolTransformerTokenizer.from_tokens(prediction_tokenized)
            prediction_split = prediction.split('.')
            pred_res = []
            for smi in prediction_split:
                try:
                    pred_res.append(rdkit_general_ops.canconicalize(smi))
                except Exception:
                    pass

            pred_res = pred_res[:1]

            out = multiset.Multiset(pred_res)
            return out, pred_score

        try:
            op_back = return_list[0]
        except KeyError as ex:
            print(return_list)
            raise ex

        output = []
        scores = []
        for elem in op_back:
            this_output, pred_score = output_mapper(elem)
            output.append(this_output)
            scores.append(pred_score)
        

        for in_, out_ in zip(list_of_reactant_sets, output):
            self.logger.debug(f"{'.'.join(in_)}>>{'.'.join(out_)}")

        # If the Transformer returned an empty reaction then randomly select one of the inputs.
        def select_input_if_op_none(input_, output_):
            output_ = [e for e in output_ if e != '']
            if len(output_):
                return output_
            else:
                return multiset.Multiset(list(input_)[:1])

        output = [select_input_if_op_none(input_elem, op_elem)
                  for input_elem, op_elem in zip(list_of_reactant_sets, output)]

        
        # return output, scores
        return output



class ChemformerModelPredictor(AbstractReactionPredictor):
    logger = logging.getLogger('Chemformer-react-predict')
    
    def __init__(self, model_id: int = 1, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.chemformer_server_address = 'http://localhost:3000/predict'
        self.grt_server_address = 'http://localhost:6006/predict'

    def _run_list_of_reactant_sets(self, list_of_reactant_sets: typing.List[multiset.FrozenMultiset]) -> typing.List[multiset.Multiset]:
        n = len(list_of_reactant_sets)
       
        list_of_reactants = [0 for _ in range(n)]
      
        for i in range(n):
            smiles_strings = list(list_of_reactant_sets[i]) 
            list_of_reactants[i] = multiset.FrozenMultiset([smiles_strings[0], smiles_strings[1]])

        def serialize_multisets(list_of_multisets):
            def multiset_to_dict(multiset):
                smiles_list = list(multiset)
                return {
                    'head_smiles': smiles_list[0],
                    'tail_smiles': smiles_list[1] 
                }

            return [multiset_to_dict(mset) for mset in list_of_multisets]

        json_to_send = json.dumps(serialize_multisets(list_of_reactants))

        def request_func(server_address):
            try:
                headers = {'Content-Type': 'application/json'}
                response = requests.post(server_address, data=json_to_send, headers=headers, timeout=180)
                response.raise_for_status()
                return response.json()
            except requests.exceptions.HTTPError as ex:
                self.logger.error(f"HTTP error occurred: {ex}")
                return None
            except requests.exceptions.RequestException as ex:
                self.logger.error(f"Request failed: {ex}")
                return None

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            future_chemformer = executor.submit(request_func, self.chemformer_server_address)
            future_grt = executor.submit(request_func, self.grt_server_address)
            chemformer_result = future_chemformer.result()
            grt_result = future_grt.result()
        

        if chemformer_result is None or grt_result is None:
            self.logger.warning("Skipping batch due to repeated request failures.")
            return None

        try:
            output = []
            for idx in range(len(chemformer_result['predictions'])):
                print(chemformer_result['predictions'][idx], grt_result['predictions'][idx])
                if are_molecules_identical(chemformer_result['predictions'][idx], grt_result['predictions'][idx]):
                    
                    output.append([chemformer_result['predictions'][idx]])
                else:
                    output.append([chemformer_result['head_smiles'][idx]])
            return output
        except KeyError as ex:
            self.logger.error(f"Error processing server response: {ex}")
            return None

# Example usage:
# list_of_reactant_sets = [
#             multiset.FrozenMultiset(['Nc1ccc(S(=O)(=O)N2CCC[C@@H](CNS(N)(=O)=O)C2)c(F)c1', 'CCCCCCCCCCC=CCCBr']),
#             multiset.FrozenMultiset(['Cc1oc2ncnc(O)c2c1C(=O)NC[C@]1(O)CCN(CCO)C1', 'CC(C)CCCCCCCCCCCN'])
#         ]
# predictor = ChemformerModelPredictor()
# result = predictor._run_list_of_reactant_sets(list_of_reactant_sets)
# print("result",result)