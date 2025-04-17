"""Sample from trained model

Usage:
  sample_from_model.py <weight_path>

"""
from os import path
from time import strftime, gmtime
import pickle

import numpy as np
from tqdm import tqdm
from docopt import docopt

from syn_dags.script_utils import doggen_utils
from syn_dags.utils import settings

OUT_DIR = "out_samples"

class Params:
    def __init__(self, weight_path: str):
        self.device = settings.torch_device()
        self.weight_path = weight_path
        self.dataversion = 'data3_80k'

        time_run = strftime("%y-%m-%d_%H:%M:%S", gmtime())
        f_name_weights = path.splitext(path.basename(self.weight_path))[0]
        self.run_name = f"doggen_sampling_on_weights_{f_name_weights}_run_at_{time_run}"
        print(f"Run name is {self.run_name}\n\n")

        self.batch_size = 20
        self.num_batches = 5

        self.log_for_reaction_predictor_path = path.join("logs", f"reactions-{self.run_name}.log")


def main(params: Params):
    # # Model (from chkpt)
    model, collate_func, model_other_parts = doggen_utils.load_doggen_model(device=params.device, server_address='123',log_path_for_react_predict=params.log_for_reaction_predictor_path,data_version=params.dataversion,weight_path=params.weight_path)

    model.to(params.device)
    # # Sample
    all_syn_trees = []
    all_log_probs = []
    for _ in tqdm(range(params.num_batches)):
        syn_trees = model.sample(params.batch_size)
        all_syn_trees.extend(syn_trees)
        # all_log_probs.append(log_probs.detach().cpu().numpy().T)
    # all_log_probs = np.concatenate(all_log_probs)

    # # Write out!
    with open(path.join(OUT_DIR, f"{params.run_name}_data3_chemforer_10epoch.pick"), 'wb') as fo:
        pickle.dump(dict( all_syn_trees=all_syn_trees), fo)

    smiles_only = [elem.root_smi for elem in all_syn_trees]
    with open(path.join(OUT_DIR,f"{params.run_name}_data3_chemforer_10epoch.txt"), 'w') as fo:
        fo.writelines('\n'.join(smiles_only))


if __name__ == '__main__':
    arguments = docopt(__doc__)
    weight_path = arguments['<weight_path>']
    weight_path = '/mnt/data/synthesis-dags/scripts/dog_gen/chkpts/epoch-26_time-25-03-26_23:03:58.pth.pick'
    main(Params(weight_path))
    print("Done!")

