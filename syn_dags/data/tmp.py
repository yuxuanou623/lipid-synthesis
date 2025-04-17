import numpy as np
import torch
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence

# Mock data: simulate a scenario where each element of sequence_choices has different lengths
sequence_choices_data = [
    np.array([0, 9, 1,8]),          # Sequence of length 3
    np.array([0,7,0, 5]),             # Sequence of length 2
    np.array([0,7,0,8,1,9])                 # Sequence of length 1
]

# Assuming p.sequence_choices for each p in pred_out would be each element in sequence_choices_data
# Initialize an empty list to collect padded sequences
sequence_choices_list = []
seq_sizes = np.array([sequence_choices.size for sequence_choices in sequence_choices_data])
array_containing_original_indcs = np.argsort(seq_sizes)[::-1]
print("array_containing_original_indcs")
print(array_containing_original_indcs)
# We need to determine the maximum sequence length for padding purposes
#seq_sizes = np.array([len(seq) for seq in sequence_choices_data])
seq_size_with_padding = max(seq_sizes)  # Maximum sequence length
print("seq_size_with_padding",seq_size_with_padding)
new_seq_sizes = seq_sizes[array_containing_original_indcs]
print("new_seq_sizes")
print(new_seq_sizes)

for new_idx, old_idx in enumerate(array_containing_original_indcs):
    print("new_idx",new_idx)
    print("old_idx",old_idx)
    p_seq_size = sequence_choices_data[old_idx].size
    sequence_choices_list.append(np.pad(sequence_choices_data[old_idx],
                                           (0, seq_size_with_padding-p_seq_size),
                                           'constant', constant_values=-100))
print("sequence_choices_list")
print(sequence_choices_list)
sequence_choices_list = torch.tensor(np.stack(sequence_choices_list), dtype=torch.int64)
print("sequence_choices_list")
print(sequence_choices_list)
seq_sizes = torch.tensor(new_seq_sizes)
sequence_choices = pack_padded_sequence(sequence_choices_list, seq_sizes, batch_first=True, enforce_sorted=False)
print("sequence_choices")
print(sequence_choices)
print(sequence_choices.data)
