import pickle
import numpy as np
import joblib

obj = joblib.load("reference_population/global_pca_v1.pkl")
scaler = obj["scaler"]
pca = obj["pca"]

# Fuse StandardScaler + PCA: X @ pcamatrix.T + bias == pca.transform(scaler.transform(X))
pcamatrix = pca.components_ / scaler.scale_[np.newaxis, :]
bias = ((-scaler.mean_ / scaler.scale_) - pca.mean_) @ pca.components_.T

with open("reference_population/global_pca_v1.converted.pkl", "wb") as f:
    pickle.dump({"pcamatrix": pcamatrix, "bias": bias}, f, protocol=4)

print(f"pcamatrix: {pcamatrix.shape}, bias: {bias.shape}")
