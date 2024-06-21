#This script is used to generate lineage trees using AntibodyForests

library(Platypus)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

load("/hpc/dla_lti/dvanginneken/PLM-likelihoods/data/OVA_V7/VariantTree/Variant_tree.csv") #replicated_df

source("/hpc/dla_lti/dvanginneken/Platypus/AntibodyForests.R")
af <- AntibodyForests(VDJ = replicated_df,
                            sequence.columns = "VDJ_sequence_aa",
                            germline.columns = "VDJ_germline_aa",
                            node.features = c("),
                            construction.method = "phylo.network.default",
                            parallel = F)
save(af, file = paste0("../data/",dataset,"/AF_",dataset,"_default_HC.RData"))
