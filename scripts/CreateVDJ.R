#This script is used to build the VDJ dataframe for a given dataset and count somatic hypermutation

library(Platypus)
library(Biostrings)
library(stringdist)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

#List the sample
samples <- list.dirs(path = paste0("../data/",dataset,"/VDJ"), full.names = TRUE, recursive = FALSE)

#Build the VDJ dataframe
source("/hpc/dla_lti/dvanginneken/Platypus/VDJ_build.R")
vdj <- VDJ_build(VDJ.sample.list = samples,
                 remove.divergent.cells = T,
                 complete.cells.only = T,
                 trim.germlines = T,
                 parallel = T,
                 num.cores = 5)

#Count somatic hypermutation (hamming distance, ignore gaps)
source("SHM_functions.R")
vdj <- SHM_calculator(vdj)

save(vdj, file = paste0("../data/",dataset,"/VDJ_",dataset,".RData"))
