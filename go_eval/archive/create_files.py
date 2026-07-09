# File paths
positive_ids_file = "temporal_positive_paac.txt.ids"
negative_ids_file = "temporal_negative_paac.txt.ids"

# Indices files
indices_files = {
    "positive_train": "temporal_positive_train_indices.txt",
    "positive_test": "temporal_positive_test_indices.txt",
    "positive_val": "temporal_positive_validation_indices.txt",
    "negative_train": "temporal_negative_train_indices.txt",
    "negative_test": "temporal_negative_test_indices.txt",
    "negative_val": "temporal_negative_validation_indices.txt",
}

# Output files
output_files = {
    "positive_train": "temporal_positive_train_ids.list",
    "positive_test": "temporal_positive_test_ids.list",
    "positive_val": "temporal_positive_validation_ids.list",
    "negative_train": "temporal_negative_train_ids.list",
    "negative_test": "temporal_negative_test_ids.list",
    "negative_val": "temporal_negative_validation_ids.list",
}

# Helper to read lines
def read_lines(file):
    with open(file) as f:
        return [line.strip() for line in f if line.strip()]

# Helper to read indices as integers
def read_indices(file):
    with open(file) as f:
        return [int(line.strip()) for line in f if line.strip()]

# Load all IDs
ids_data = {
    "positive": read_lines(positive_ids_file),
    "negative": read_lines(negative_ids_file),
}

# Function to map indices to IDs
def indices_to_ids(indices, id_list):
    return [id_list[i] for i in indices]

# Process all six files
for key, idx_file in indices_files.items():
    group = "positive" if "positive" in key else "negative"
    id_list = ids_data[group]
    indices = read_indices(idx_file)
    mapped_ids = indices_to_ids(indices, id_list)
    
    # Write output
    with open(output_files[key], "w") as f:
        f.write("\n".join(mapped_ids))

print("All ID files generated successfully.")
