
import abc
import typing
import re
import collections
import json
import logging

import requests
import multiset

from ..chem_ops import rdkit_general_ops
from ..utils import misc
from ..utils import settings


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
     
       
        a = self._run_list_of_reactant_sets(list_of_reactant_sets)
        print("a",a)
        print(type(a))
        reactants_products_new = list(zip(list_of_reactant_sets,
                                              self._run_list_of_reactant_sets(list_of_reactant_sets)))
        # reactants_products_new = list(zip(list_of_reactant_sets,
        #                                       self._run_list_of_reactant_sets(list_of_reactant_sets)))
     

        # Work out which ones we can just remove out from the cache
        # NB note we POP them out so that when we reinsert them below we put them in at the end.
  

        # We collect both the old and new ones together and form the output list
        reactant_products_map = collections.OrderedDict(reactants_products_new)
        out = [reactant_products_map[item] for item in list_of_reactant_sets]

        # We also add these new results to end of our cache and remove anything at the beginning if necessary so that
        # our memory does not blow up.
      

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
        self.server_address = 'http://localhost:9090/predict'

    def _run_list_of_reactant_sets(self, list_of_reactant_sets: typing.List[multiset.FrozenMultiset]) -> typing.List[multiset.Multiset]:
        n = len(list_of_reactant_sets)
       
        
        list_of_reactants = [0 for _ in range(n)]
      
        
        for i in range(n):
            
            
            smiles_strings = list(list_of_reactant_sets[i]) 
            
            
            list_of_reactants[i] = multiset.FrozenMultiset([smiles_strings[0], smiles_strings[1]])

        def serialize_multisets(list_of_multisets):
            def multiset_to_dict(multiset):
                # Extract SMILES strings from the multiset and assign them to specific keys
                smiles_list = list(multiset)  # Convert FrozenMultiset to a list to access items
                return {
                    'head_smiles': smiles_list[0],
                    'tail_smiles': smiles_list[1] 
                }

            # This function takes a list of FrozenMultiset objects and converts each to a dictionary
            return [multiset_to_dict(mset) for mset in list_of_multisets]
        json_to_send = json.dumps(serialize_multisets(list_of_reactants))

        def request_func():
            try:
                headers = {'Content-Type': 'application/json'}
                response = requests.post(self.server_address, data=json_to_send, headers=headers, timeout=180)
                
                response.raise_for_status()  # Raises HTTPError for bad requests
                return response.json()
            except requests.exceptions.HTTPError as ex:
                self.logger.error(f"HTTP error occurred: {ex}")
                return None  # Return None to signify the request failed
            except requests.exceptions.RequestException as ex:
                self.logger.error(f"Request failed: {ex}")
                return None

        # Attempt the request up to 3 times
        result = misc.retry_n_times(request_func, 3, Exception, interval=0.5,
                                    on_exception=lambda ex: self.logger.error(f"Request error: {ex}"))
                                 

        if result is None:
            self.logger.warning("Skipping batch due to repeated request failures.")
            return None  # You can return an appropriate default or None to indicate skipping

        # Process the result
        try:
            # Assuming result processing code here
            products = result['predictions']
            probabilities = result['probabilities']
            # Further processing and mapping...
            output = [[p] for p in products]
           
            return output
        except KeyError as ex:
            self.logger.error(f"Error processing server response: {ex}")
            return None  # Handle missing keys or other issues gracefully


class TemplatePredictor(AbstractReactionPredictor):
    logger = logging.getLogger('Template-react-predict')
    
    def __init__(self, model_id: int = 1, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.server_address = 'http://localhost:8080/reaction'

    def _run_list_of_reactant_sets(self, list_of_reactant_sets: typing.List[multiset.FrozenMultiset]) -> typing.List[multiset.Multiset]:
        print("in _run_list_of_reactant_sets temp")
        n = len(list_of_reactant_sets)
       
        
        list_of_reactants = [0 for _ in range(n)]
      
        
        for i in range(n):
            
            
            smiles_strings = list(list_of_reactant_sets[i]) 
            
            
            list_of_reactants[i] = multiset.FrozenMultiset([smiles_strings[0], smiles_strings[1]])

        def serialize_multisets(list_of_multisets):
            def multiset_to_dict(multiset):
                # Convert multiset to a sorted list to ensure consistent head/tail order
                smiles_list = sorted(list(multiset))  # optional: sort if unordered
                if len(smiles_list) < 2:
                    raise ValueError(f"Multiset must have at least 2 elements, got: {smiles_list}")
              
                return {
                    'head': smiles_list[0],
                    'tail_smiles1': smiles_list[1]
                }

            # Zip the reactant multisets with corresponding templates
         
            return [multiset_to_dict(mset) for mset in list_of_multisets]
        json_to_send = json.dumps(serialize_multisets(list_of_reactants))

        def request_func():
            try:
                headers = {'Content-Type': 'application/json'}
                response = requests.post(self.server_address, data=json_to_send, headers=headers, timeout=180)
                
                response.raise_for_status()  # Raises HTTPError for bad requests
                return response.json()
            except requests.exceptions.HTTPError as ex:
                print("here1")
                self.logger.error(f"HTTP error occurred: {ex}")
                return None  # Return None to signify the request failed
            except requests.exceptions.RequestException as ex:
                print("here2")
                self.logger.error(f"Request failed: {ex}")
                return None

        # Attempt the request up to 3 times
        print("here")
        results = misc.retry_n_times(request_func, 3, Exception, interval=0.5,
                                    on_exception=lambda ex: self.logger.error(f"Request error: {ex}"))
        print("result",results)
                                 

        if results is None:
            print("resut", type(results))
            print("result is None")
            self.logger.warning("Skipping batch due to repeated request failures.")
            return None  # You can return an appropriate default or None to indicate skipping

        # Process the result
        try:
            # Assuming result processing code here
            
        
            # Further processing and mapping...
            
          
              
            
            print("111")
            print("type(results", type(results))
            for result in results:
                print(type(result))
                print(result['products'])
            output = [[result['products']] for result in results]
            print("output", output)
           
            return output
        except KeyError as ex:
            print("keyerr", ex)

            self.logger.error(f"Error processing server response: {ex}")
            return None  # Handle missing keys or other issues gracefully
        except Exception as ex:
            print(f"Unexpected error while processing server response: {ex}")
            self.logger.error(f"Unexpected error while processing server response: {ex}")
            return None  # Handle all other exceptions
