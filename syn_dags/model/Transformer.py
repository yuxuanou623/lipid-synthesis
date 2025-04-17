import torch
import torch.nn as nn
import math

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super(PositionalEncoding, self).__init__()
        # Dropout module to improve training stability
        self.dropout = nn.Dropout(p=0.1)

        # Create a position array with shape [max_len, 1]
        position = torch.arange(max_len).unsqueeze(1)
        
        # Compute the div_term for the positional encoding formula
        # Using the formula from the "Attention is All You Need" paper
        div_term = torch.exp(torch.arange(0, d_model, 2) * -(math.log(10000.0) / d_model))
        
        # Initialize a positional encoding matrix with zeros
        pe = torch.zeros(max_len, d_model)
        
        # Apply the sine function to even indices in the positional encoding matrix
        pe[:, 0::2] = torch.sin(position * div_term)
        
        # Apply the cosine function to odd indices in the positional encoding matrix
        pe[:, 1::2] = torch.cos(position * div_term)
        
        # Add an extra dimension to pe to facilitate batch processing
        # pe shape is now [1, max_len, d_model] making it easy to add to any batch of embeddings
        self.register_buffer('pe', pe.unsqueeze(0))

    def forward(self, x):
        # x is the input embedding with shape [batch_size, seq_len, d_model]
        
        # Add positional encodings to the input embedding tensor
        # Using slicing to adjust the positional encoding to match the input sequence length
        x = x + self.pe[:, :x.size(1)]
        
        # Apply dropout to the combined embeddings to promote model generalization
        return self.dropout(x)

class AutoregressiveTransformer(nn.Module):
    def __init__(self, d_model, n_head, num_layers, d_ffn, max_seq_len, vocab_size):
        super(AutoregressiveTransformer, self).__init__()
        self.pos_encoder = PositionalEncoding(d_model, max_seq_len)
        decoder_layer = nn.TransformerDecoderLayer(d_model, n_head, d_ffn)
        self.decoder = nn.TransformerDecoder(decoder_layer, num_layers)
        self.out = nn.Linear(d_model, vocab_size)
        self.head = n_head

    def forward(self, src, src_mask):
        src = self.pos_encoder(src)
        print("src.size()",src.size())
        print("src_mask",src_mask.size())
        output = self.decoder(src, memory=None, tgt_mask=src_mask)
        #output = self.decoder(src, memory=None)
        return self.out(output)

    def generate_square_subsequent_mask_self( sz):
        mask = (torch.triu(torch.ones(sz, sz)) == 1).transpose(0, 1)
        mask = mask.float().masked_fill(mask == 0, float('-inf')).masked_fill(mask == 1, float(0.0))
        return mask
