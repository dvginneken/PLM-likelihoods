import matplotlib
import torch
from evolocity.tools.fb_model import FBModel
import argparse
import pandas as pd
import anndata
import evolocity as evo
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--model", help="Choose a model: esm, ablang, protbert, sapiens")
parser.add_argument("-s", "--sample_id")
parser.add_argument("-i", "--input_source", help="cdr3_only or full_VDJ")
parser.add_argument("-d", "--dataset")
parser.add_argument("-c", "--color", help="Choose from: SHM_count or c_gene")

args = parser.parse_args()
model = args.model
sample_id = args.sample_id
input_source = args.input_source
dataset = args.dataset
color = args.color

#Get embeddings, substring the barcode, only keep heavy chain, and remove NA values from dimensions and barcodes
embedding_file = os.path.join("..","data",dataset,"VDJ",sample_id,"embeddings",input_source,f"embeddings_{model}.csv.gzip")
embeddings_df = pd.read_csv(embedding_file, compression="gzip")
embeddings_df['barcode'] = embeddings_df['barcode'].str[:-2] 
embeddings_df = embeddings_df[embeddings_df['chain'] == "IGH"]
embeddings_df = embeddings_df.dropna(subset=["dim_0"])
embeddings_df = embeddings_df.dropna(subset=["barcode"])

#Get VDJ information for a specific sample
vdj_file = os.path.join("..","data",dataset,"vdj_evolike_combine.csv")
vdj_df = pd.read_csv(vdj_file)
vdj_df = vdj_df[vdj_df['sample_id'] == sample_id]

#Match embeddings and vdj metadata
df_merged = pd.merge(embeddings_df, vdj_df, on='barcode', how='inner')

#Get the sequences
if input_source == "full_VDJ":
    sequences = [str(gene) for gene in df_merged['VDJ_sequence_aa_trimmed']]
if input_source == "cdr3_only":
    sequences = [str(gene) for gene in df_merged['VDJ_cdr3_aa']]

#Get other metadata to plot
vdj_cgene = [str(gene) for gene in df_merged['c_gene']]
shm_count = [float(count) for count in df_merged["SHM_count"]]

#Get embedding columns and metadata columns
embedding_cols = [col for col in df_merged.columns if col.startswith('dim')]
metadata_cols = [col for col in df_merged.columns if col not in embedding_cols]

#Create an AnnData object
adata = anndata.AnnData(df_merged[embedding_cols])

#Add sequence and metadata information to .obs slot of AnnData object
adata.obs['seq'] = sequences
adata.obs['c_gene'] = vdj_cgene
#adata.obs['sample_id'] = sample_ids
adata.obs["SHM_count"] = shm_count

#Construct sequence similarity network
evo.pp.neighbors(adata)
sc.tl.umap(adata)
basis = "umap"
if model  == "esm":
    evo.tl.velocity_graph(adata)
else:
    evo.tl.velocity_graph(adata, model_name = model)
if 'model' in adata.uns:
    del adata.uns['model']

#Embed network and velocities in two-dimensions
evo.tl.velocity_embedding(adata, basis = basis)

#Save the processed AnnData object to a file
output_file = os.path.join("..","data",dataset,"evo-velocity",f"adata_{sample_id}_{input_source}_{model}.h5ad")
adata.write(output_file)
print(f"Processed data saved to {output_file}")

# Plot your data using the color mapping (jgene or SHM_count), make sure data is categorical
ax = evo.pl.velocity_embedding_stream(adata, color=color, legend_loc="right margin", show=False)
save_path = os.path.join("..","data",dataset,"evo-velocity",f"embedding_{color}_{sample_id}_{input_source}_{model}.png")
plt.savefig(save_path, bbox_inches="tight")
