import torch
from transformers import AutoModelForMaskedLM
model = AutoModelForMaskedLM.from_pretrained('/groups/clairemcwhite/models/ESMplusplus_large', trust_remote_code=True)
model.eval()
tokenizer = model.tokenizer

sequences = ['MPRTEIN', 'MSEQWENCE']
tokenized = tokenizer(sequences, padding=True, return_tensors='pt')

# tokenized['labels'] = tokenized['input_ids'].clone() # correctly mask input_ids and set unmasked instances of labels to -100 for MLM training

print(tokenizer.convert_ids_to_tokens(range(30)))

output = model(**tokenized) # get all hidden states with output_hidden_states=True
print(output.logits.shape) # language modeling logits, (batch_size, seq_len, vocab_size), (2, 11, 64)
print(output.logits)
print(output.logits.argmax(-1))
print(output.last_hidden_state.shape) # last hidden state of the model, (batch_size, seq_len, hidden_size), (2, 11, 1152)
print(output.loss) # language modeling loss if you passed labels
#print(output.hidden_states) # all hidden states if you passed output_hidden_states=True (in tuple)


model.eval()

with torch.no_grad():
    outputs = model(**tokenized)

logits = outputs.logits
pred_ids = logits.argmax(dim=-1)
print(pred_ids)


import torch

model.eval()
tokenized = tokenizer(sequences, padding=True, return_tensors="pt")

mask_id = tokenizer.mask_token_id

# Mask position 3 in first sequence
tokenized["input_ids"][0, 3] = mask_id

with torch.no_grad():
    output = model(**tokenized)

pred = output.logits.argmax(dim=-1)

print(tokenizer.convert_ids_to_tokens(pred[0]))

mask_token_id = tokenizer.mask_token_id

# Mask one position per sequence
tokenized["input_ids"][0, 3] = mask_token_id
tokenized["input_ids"][1, 4] = mask_token_id

with torch.no_grad():
    outputs = model(**tokenized)

logits = outputs.logits

# Get predictions only at masked positions
masked_positions = tokenized["input_ids"] == mask_token_id
pred_ids = logits.argmax(dim=-1)

print(pred_ids[masked_positions])
print(tokenizer.convert_ids_to_tokens(pred_ids[masked_positions].tolist()))

