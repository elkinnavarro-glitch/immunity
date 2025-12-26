import scanpy as sc
import os

path = r"C:\gemini\SkinTissuePrograms\Atlas\03_processed_data\2021_skin\2021_skin_raw.h5ad"

if os.path.exists(path):
    adata = sc.read_h5ad(path)
    print(f"Shape: {adata.shape}")
    print(f"Obs columns: {adata.obs.columns.tolist()}")
    print(f"Counts X[0,0]: {adata.X[0,0]}")
    print(f"Processing Log: {adata.uns.get('processing_log', 'Missing')}")
else:
    print("File not found.")
