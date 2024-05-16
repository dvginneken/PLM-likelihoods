library(viridis)
library(igraph)

load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/AF_OVA_V7_MP_HC.RData") #af
af_ova <- af
rm(af)
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/AF_horns2020a__VDJ_RAW_MP_HC.RData")
af_horns <- af
rm(af)
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/AF_Bruhn_MP_HC.RData")
af_bruhn <- af
rm(af)
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/AF_Bieberich_MP_HC.RData")
af_bieberich <- af
rm(af)

plot_tree <- function(af, sample, clonotype, color.by){
  tree <- af[[sample]][[clonotype]][["igraph"]]
  nodes <- af[[sample]][[clonotype]][["nodes"]]
  
  igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", "G", 
                                  ifelse(startsWith(igraph::V(tree)$name, "node"), gsub(pattern = "node", replacement = "", igraph::V(tree)$name), igraph::V(tree)$name))
  nodes <- lapply(nodes, function(x){
    if(is.null(x[["size"]])){x[["size"]] <- 1};return(x)
  })
  igraph::V(tree)$size <- 30 + as.numeric(lapply(nodes,function(x){x[["size"]]})[names(igraph::V(tree))])
  igraph::V(tree)$size <- igraph::V(tree)$size * 0.3
  layout <- igraph::layout_as_tree(tree, root = "germline")
  color_values <- lapply(nodes,function(x){unique(x[[color.by]])})[names(igraph::V(tree))]
  evo <- c("germline" = NA, unlist(color_values))
  V(tree)$color <- scales::dscale(evo %>% cut(length(evo)), viridis_pal())
  igraph::plot.igraph(tree, layout = layout,
                      vertex.label = igraph::V(tree)$label,
                      edge.arrow.size = 0.1)
}


plot_tree(af_ova,
          sample = "S1",
          clonotype = "clonotype4",
          color.by = "IGH_evo_likelihood_esm_full_VDJ")
plot_tree(af_horns,
          sample = "Influenza.vac.11.12.human.S1",
          clonotype = "clonotype3",
          color.by = "IGH_evo_likelihood_esm_full_VDJ")
plot_tree(af_bruhn,
          sample = "Bruhn",
          clonotype = "clonotype4",
          color.by = "IGH_evo_likelihood_esm_full_VDJ")
plot_tree(af_bieberich,
          sample = "S2",
          clonotype = "clonotype1",
          color.by = "IGH_evo_likelihood_esm_full_VDJ")

# source("~/Documents/GitHub/Platypus/R/AntibodyForests_plot.R")
# AntibodyForests_plot(af_ova,
#                      sample = "S1",
#                      clonotype = "clonotype2",
#                      color.by = "IGH_evo_likelihood_esm_full_VDJ")
