library(viridis)
library(igraph)
library(tidyr)
library(dplyr)
library(ggplot2)

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

plot_tree <- function(af, sample, clonotype, color.by, scale.size = 0.4){
  tree <- af[[sample]][[clonotype]][["igraph"]]
  nodes <- af[[sample]][[clonotype]][["nodes"]]
  

  nodes <- lapply(nodes, function(x){
    if(is.null(x[["size"]])){x[["size"]] <- 1};return(x)
  })
  igraph::V(tree)$label <- ifelse(igraph::V(tree)$name == "germline", "G", 
                                  ifelse(startsWith(igraph::V(tree)$name, "node"), lapply(nodes,function(x){x[["size"]]})[names(igraph::V(tree))], igraph::V(tree)$name))
  
  igraph::V(tree)$size <- 10 + (as.numeric(lapply(nodes,function(x){x[["size"]]})[names(igraph::V(tree))]) * scale.size)
  layout <- igraph::layout_as_tree(tree, root = "germline")
  color_values <- lapply(nodes,function(x){unique(x[[color.by]])})[names(igraph::V(tree))]
  evo <- c("germline" = NA, unlist(color_values))[V(tree)$name]
  V(tree)$color <- scales::dscale(evo %>% cut(length(evo)-1), viridis_pal())
  igraph::plot.igraph(tree, layout = layout,
                      vertex.label = igraph::V(tree)$label,
                      edge.arrow.size = 0.1)
}

legend_image <- as.raster(viridis(100, direction = -1), ncol=1)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pseudolikelihood', cex.main = 2)
text(x=1, y = seq(0,1), labels = c("Low","High"), cex=2)
rasterImage(legend_image, 0.2, 0.05, 1.5,0.95)

plot_tree(af_ova,
          sample = "S1",
          clonotype = "clonotype4",
          color.by = "IGH_evo_likelihood_esm_full_VDJ",
          scale.size = 0.7)
# plot_tree(af_ova,
#           sample = "S2",
#           clonotype = "clonotype4",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ",
#           scale.size = 0.7)
# plot_tree(af_ova,
#           sample = "S4",
#           clonotype = "clonotype7",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ",
#           scale.size = 0.7)
# plot_tree(af_ova,
#           sample = "S5",
#           clonotype = "clonotype4",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ",
#           scale.size = 0.7)
# plot_tree(af_horns,
#           sample = "Influenza.vac.11.12.human.S1",
#           clonotype = "clonotype6",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ")
plot_tree(af_bruhn,
          sample = "Bruhn",
          clonotype = "clonotype1",
          color.by = "IGH_evo_likelihood_esm_full_VDJ")
# plot_tree(af_bruhn,
#           sample = "Bruhn",
#           clonotype = "clonotype2",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ")
# plot_tree(af_bruhn,
#           sample = "Bruhn",
#           clonotype = "clonotype3",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ")
# plot_tree(af_bruhn,
#           sample = "Bruhn",
#           clonotype = "clonotype4",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ")
# plot_tree(af_bieberich,
#           sample = "S7",
#           clonotype = "clonotype2",
#           color.by = "IGH_evo_likelihood_esm_full_VDJ")

# source("~/Documents/GitHub/Platypus/R/AntibodyForests_plot.R")
# AntibodyForests_plot(af_ova,
#                      sample = "S1",
#                      clonotype = "clonotype4",
#                      color.by = "IGH_evo_likelihood_esm_full_VDJ")

#Distance to germline

#Function to get a dataframe with the distance to the germline for each node and the pseudolikelihoods
distance_germline <- function(input){
  #Go over each tree in the AntibodyForests object and create a metric dataframe
  metric_df <- lapply(seq_along(input), function(sample){
    lapply(seq_along(input[[sample]]), function(clonotype){
      sample_name <- names(input)[sample]
      clonotype_name <- names(input[[sample_name]])[clonotype]
      
      tree = input[[sample_name]][[clonotype_name]][["igraph"]]
      nodes = igraph::V(tree)[names(igraph::V(tree)) != "germline"]
      node_features = input[[sample_name]][[clonotype_name]][["nodes"]][names(input[[sample_name]][[clonotype_name]][["nodes"]]) != "germline"]
      
      #Get the total length of shortest paths between each node and the germline
      distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                                    weights = edge_attr(tree)$edge.length)
      distance = t(as.data.frame(distance))
      
      
      #Evo likelihood
      esm = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_esm_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_esm_full_VDJ"]]))]}))[rownames(distance)]
      protbert = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_protbert_full_VDJ"]]))]}))[rownames(distance)]
      ablang = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_ablang_full_VDJ"]]))]}))[rownames(distance)]
      sapiens = unlist(lapply(node_features,function(x){unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]])[!is.na(unique(x[["IGH_evo_likelihood_sapiens_full_VDJ"]]))]}))[rownames(distance)]
      
      
      distance = data.frame("distance" = distance, "sample" = sample_name, "clonotype" = clonotype_name, 
                            "node" = rownames(distance), "esm" = esm, "protbert" = protbert, 
                            "ablang" = ablang, "sapiens" = sapiens)
      return(distance)
    })
  })
  #Transform list into dataframe
  dfs <- lapply(metric_df, function(x){do.call(rbind, x)})
  df <- do.call(rbind, dfs)
  rownames(df) <- paste0(df$sample, "_", df$clonotype, "_", df$node)
  
  return(df)
}

#Get distance to germline for all samples
distance_ova <- distance_germline(af_ova)
distance_horns <- distance_germline(af_horns)
distance_bruhn <- distance_germline(af_bruhn)
distance_bieberich <- distance_germline(af_bieberich)

#change sample names
distance_ova$sample <- case_match(distance_ova$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")

distance_horns$sample <- case_match(distance_horns$sample,
                              "Influenza.vac.11.12.human.S1" ~ "Individual1",
                              "Influenza.vac.11.12.human.S2" ~ "Human2",
                              "Influenza.vac.11.12.human.S3" ~ "Human3",
                              "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
distance_horns <- distance_horns[distance_horns$sample == "Individual1",]

distance_bruhn$sample <- "Individual2"

distance_bieberich$sample <- case_match(distance_bieberich$sample,
                                  "S2" ~ "Individual3",
                                  "S4" ~ "Individual4",
                                  "S7" ~ "Individual5")

distance_all <- rbind(distance_ova, distance_horns, distance_bruhn, distance_bieberich)
distance_human <- rbind(distance_horns, distance_bruhn, distance_bieberich)

#Plot the distance to the germline against the pseudolikelihoods
ggplot(distance_ova, aes(x = germline, y = esm, colour = sample)) + 
  geom_point() + 
  theme_minimal()

ggplot(distance_horns, aes(x = germline, y = esm, colour = sample)) + 
  geom_point() + 
  theme_minimal()

ggplot(distance_bruhn, aes(x = germline, y = esm, colour = sample)) + 
  geom_point() + 
  theme_minimal()

ggplot(distance_bieberich, aes(x = germline, y = esm, colour = sample)) + 
  geom_point() + 
  theme_minimal()

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

#esm
cor_human <- cor.test(distance_human$germline, distance_human$esm)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$esm)$estimate
ggplot(distance_all, aes(x = germline, y = esm, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Distance to germline") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#protbert
cor_human <- cor.test(distance_human$germline, distance_human$protbert)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$protbert)$estimate
ggplot(distance_all, aes(x = germline, y = protbert, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Distance to germline") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#sapiens
cor_human <- cor.test(distance_human$germline, distance_human$sapiens)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$sapiens)$estimate
ggplot(distance_all, aes(x = germline, y = sapiens, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Distance to germline") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(distance_human$germline, distance_human$ablang)$estimate
cor_mouse <- cor.test(distance_ova$germline, distance_ova$ablang)$estimate
ggplot(distance_all, aes(x = germline, y = ablang, colour = sample)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values = c("Individual1" = "#99d8c9",
                                "Individual2" = "#66c2a4",
                                "Individual3" = "#41ae76",
                                "Individual4" = "#238b45",
                                "Individual5" = "#005824",
                                "Mouse1" = "#fcc5c0",
                                "Mouse2" = "#fa9fb5",
                                "Mouse3" = "#f768a1",
                                "Mouse4" = "#c51b8a",
                                "Mouse5" = "#7a0177"),
                     name = "Sample") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Distance to germline") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))


