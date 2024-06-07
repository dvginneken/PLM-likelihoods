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
parser.add_argument("-d", "--dataset")

args = parser.parse_args()
model = args.model
dataset = args.dataset

#Get embeddings, substring the barcode, only keep heavy chain, and remove NA values from dimensions and barcodes
embedding_file = os.path.join("..","data",dataset,"all_germline_embeddings",f"all_germline_embeddings_{model}.csv.gzip")
embeddings_df = pd.read_csv(embedding_file, compression="gzip")
embeddings_df = embeddings_df[embeddings_df['VDJ_chain'] == "IGH"]
embeddings_df = embeddings_df.dropna(subset=["dim_0"])
embeddings_df = embeddings_df.dropna(subset=["barcode"])

#Get VDJ information
vdj_file = os.path.join("..","data",dataset,"vdj_evolike_combine.csv")
vdj_df = pd.read_csv(vdj_file)

#Match embeddings and vdj metadata
df_merged = pd.merge(embeddings_df, vdj_df, on='barcode', how='inner')

#Get the sequences
sequences = [str(gene) for gene in df_merged['VDJ_germline_aa_trimmed_x']]

#Get other metadata to plot
IgG_subtypes = ["IGHG", "IGHG1","IGHG2","IGHG2B","IGHG2C","IGHG3","IGHG4"]
IgA_subtypes = ["IGHA","IGHA1","IGHA2"]
df_merged["c_gene"] = df_merged["VDJ_cgene_x"].replace(IgG_subtypes,"IgG").replace(IgA_subtypes,"IgA").replace("IGHM","IgM").replace("IGHD","IgD").replace("IGHE","IgE")
isotype = [str(gene) for gene in df_merged['c_gene']]
sample = [str(s) for s in df_merged['sample_id_x']]
clonotype = [str(cl) for cl in df_merged['clonotype_id_x']]
df_merged['v_gene'] = df_merged['VDJ_vgene_x'].apply(lambda x: x.split('-')[0])
v_gene_family = [str(gene) for gene in df_merged['v_gene']]
  

#Get embedding columns and metadata columns
embedding_cols = [col for col in df_merged.columns if col.startswith('dim')]
metadata_cols = [col for col in df_merged.columns if col not in embedding_cols]

#Create an AnnData object
adata = anndata.AnnData(df_merged[embedding_cols])

#Add sequence and metadata information to .obs slot of AnnData object
adata.obs['seq'] = sequences
adata.obs['isotype'] = isotype
adata.obs["sample"] = sample
adata.obs["clonotype"] = clonotype
adata.obs["v_gene_family"] = v_gene_family

#Remove duplicate observations
adata = adata[~adata.to_df().duplicated(), :]

#Construct sequence similarity network
evo.pp.neighbors(adata)
sc.tl.umap(adata)
basis = "umap"
print(adata.to_df().duplicated().sum())
if model  == "esm":
    evo.tl.velocity_graph(adata)
else:
    evo.tl.velocity_graph(adata, model_name = model)
if 'model' in adata.uns:
    del adata.uns['model']

#Embed network and velocities in two-dimensions
evo.tl.velocity_embedding(adata, basis = basis)

#Save the processed AnnData object to a file
output_file = os.path.join("..","data",dataset,"all_germline_embeddings",f"adata_all_germline_embeddings_{model}.h5ad")
adata.write(output_file)
print(f"Processed data saved to {output_file}")
