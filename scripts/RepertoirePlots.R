library(ggplot2)
library(ggpubr)
#library(Platypus)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
#Read VDJ dataframes
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/VDJ_PLL_OVA_V7.RData")
vdj_ova <- vdj
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/VDJ_PLL_horns2020a__VDJ_RAW.RData")
vdj_horns <- vdj
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/VDJ_PLL_Bruhn.RData")
vdj_bruhn <- vdj
rm(vdj)
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/VDJ_PLL_Bieberich.RData")
vdj_bieberich <- vdj
rm(vdj)

#Change sample names
vdj_ova$sample <- case_match(vdj_ova$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")
vdj_ova$sample_id <- vdj_ova$sample

vdj_horns$sample <- case_match(vdj_horns$sample,
                              "Influenza.vac.11.12.human.S1" ~ "Individual1",
                              "Influenza.vac.11.12.human.S2" ~ "Human2",
                              "Influenza.vac.11.12.human.S3" ~ "Human3",
                              "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
vdj_horns <- vdj_horns[vdj_horns$sample == "Individual1",]
vdj_horns$sample_id <- vdj_horns$sample

vdj_bruhn$sample <- "Individual2"
vdj_bruhn$sample_id <- vdj_bruhn$sample

vdj_bieberich$sample <- case_match(vdj_bieberich$sample,
                                  "S2" ~ "Individual3",
                                  "S4" ~ "Individual4",
                                  "S7" ~ "Individual5")
vdj_bieberich$sample_id <- vdj_bieberich$sample

#Set v-gene families
vdj_ova$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_ova$VDJ_vgene)
vdj_bruhn$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_bruhn$VDJ_vgene)
vdj_horns$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_horns$VDJ_vgene)
vdj_bieberich$v_gene_family <- gsub(pattern = "-.*", replacement = "", x = vdj_bieberich$VDJ_vgene)



#Remove NA cgene
vdj_ova <- vdj_ova[!is.na(vdj_ova$VDJ_cgene),]
vdj_bieberich <- vdj_bieberich[!is.na(vdj_bieberich$VDJ_cgene),]
vdj_horns <- vdj_horns[!is.na(vdj_horns$VDJ_cgene),]
vdj_bruhn <- vdj_bruhn[!is.na(vdj_bruhn$VDJ_cgene),]

#combine human samples
vdj_human <- rbind(vdj_horns, vdj_bruhn, vdj_bieberich)
vdj_all <- rbind(vdj_ova, vdj_human)

#clonal expansion
source("~/OneDrive - UMC Utrecht/Documenten/Platypus/VDJ_clonal_expansion_Daphne.R")
VDJ_clonal_expansion_Daphne(vdj_bruhn, clones = 30, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)
VDJ_clonal_expansion_Daphne(vdj_horns, clones = 30, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)
VDJ_clonal_expansion_Daphne(vdj_ova, clones = 50, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)
VDJ_clonal_expansion_Daphne(vdj_bieberich, clones = 50, group.by = "sample_id", color.by = "VDJ_cgene", text.size=20)


source("~/OneDrive - UMC Utrecht/Documenten/Platypus/VDJ_clonal_barplot.R")
vdj_human$sample_id <- gsub(pattern = "Individual", replacement = "I", x = vdj_human$sample_id)
vdj_ova$sample_id <- gsub(pattern = "Mouse", replacement = "M", x = vdj_ova$sample_id)
out1 <- VDJ_clonal_barplot(vdj_human, counts.to.use = "clonotype_id", group.by = "sample_id", expanded.colors = brewer.pal(n=5, "BuGn")[2:5])
out2 <- VDJ_clonal_barplot(vdj_ova, counts.to.use = "clonotype_id", group.by = "sample_id", expanded.colors = brewer.pal(n=5, "RdPu")[2:5])
do.call("grid.arrange", c(out1, out2, ncol = 10))


nrow(vdj_all)

sub_colors <- c("IGHA" = "red",
               "IGHA1" = "#fb6a4a",
               "IGHA2" = "#de2d26",
               "IGHE" = "purple",
               "IGHG1" = "#a1d99b",
               "IGHG2" = "#74c476",
               "IGHG2B" = "#41ab5d",
               "IGHG2C" = "#238b45",
               "IGHG3" = "#005a32",
               "IGHM" = "black",
               "IGHD" = "blue",
               "NA" = "grey")

isotype_colors <- c("IgA" = "red",
               "IgE" = "purple",
               "IgG" = "green",
               "IgM" = "black",
               "IgD" = "blue",
               "NA" = "grey")

vdj_ova$isotype <- case_match(vdj_ova$isotype,
                              "IgG1" ~ "IgG",
                              "IgG2" ~ "IgG",
                              "IgG3" ~ "IgG",
                              "IGHA" ~ "IgA",
                              "IGHD" ~ "IgD",
                              "IGHE" ~ "IgE",
                              "IGHM" ~ "IgM")
vdj_human$isotype <- case_match(vdj_human$isotype,
                              "IgG1" ~ "IgG",
                              "IgG2" ~ "IgG",
                              "IgG3" ~ "IgG",
                              "IgG4" ~ "IgG",
                              "IGHA" ~ "IgA",
                              "IGHD" ~ "IgD",
                              "IGHE" ~ "IgE",
                              "IGHM" ~ "IgM",
                              "IgA1" ~ "IgA",
                              "IgA2" ~ "IgA")

#Isotype
#Mouse foundational PLMs
a <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Full VDJ")

b <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_cdr3_from_VDJ, y=IGH_evo_likelihood_protbert_cdr3_from_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("CDR3 from VDJ")

c <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_cdr3_only, y=IGH_evo_likelihood_protbert_cdr3_only, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("CDR3 only")

ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#Mouse antibody-specific PLMs
a <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=IGH_evo_likelihood_ablang_full_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("Full VDJ")

b <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_sapiens_cdr3_from_VDJ, y=IGH_evo_likelihood_ablang_cdr3_from_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("CDR3 from VDJ")

c <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_sapiens_cdr3_only, y=IGH_evo_likelihood_ablang_cdr3_only, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("CDR3 only")

ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#Human foundational PLMs
a <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Full VDJ")

b <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_cdr3_from_VDJ, y=IGH_evo_likelihood_protbert_cdr3_from_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("CDR3 from VDJ")

c <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_cdr3_only, y=IGH_evo_likelihood_protbert_cdr3_only, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("CDR3 only")

ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#Human antibody-specific PLMs
a <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=IGH_evo_likelihood_ablang_full_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("Full VDJ")

b <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_sapiens_cdr3_from_VDJ, y=IGH_evo_likelihood_ablang_cdr3_from_VDJ, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("CDR3 from VDJ")

c <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_sapiens_cdr3_only, y=IGH_evo_likelihood_ablang_cdr3_only, color=isotype)) +
  geom_point(size = 1) +
  scale_color_manual(name = "Isotype", values = isotype_colors) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("CDR3 only")

ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#V-gene family
vgene_colors <- c("IGHV1" = "coral1",
                  "IGHV2" = "cadetblue",
                  "IGHV3" = "darkolivegreen3",
                  "IGHV4" = "darkorchid2",
                  "IGHV5" = "darkgoldenrod4",
                  "IGHV6" = "snow4",
                  "IGHV7" = "gold3",
                  "IGHV8" = "hotpink",
                  "IGHV9" = "dodgerblue3",
                  "IGHV10" = "violetred4",
                  "IGHV11" = "aquamarine4",
                  "IGHV12" = "deeppink",
                  "IGHV13" = "darkgreen",
                  "IGHV14" = "orangered",
                  "IGHV15" = "mediumpurple")
#Human foundational PLMs
a1 <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=v_gene_family)) +
  geom_point(size = 1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon()+
  theme(text = element_text(size = 18)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Human")

#Human antibody-specific PLMs
b1 <- ggplot(vdj_human, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=IGH_evo_likelihood_ablang_full_VDJ, color=v_gene_family)) +
  geom_point(size = 1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("Human")

#Mouse foundational PLMs
a2 <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=IGH_evo_likelihood_protbert_full_VDJ, color=v_gene_family)) +
  geom_point(size = 1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon()+
  theme(text = element_text(size = 18)) +
  xlab("ESM-1b Pseudolikelihood") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle("Mouse")

#Mouse antibody-specific PLMs
b2 <- ggplot(vdj_ova, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=IGH_evo_likelihood_ablang_full_VDJ, color=v_gene_family)) +
  geom_point(size = 1) +
  scale_color_manual(name = "V-gene family", values = vgene_colors,
                     breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
                                "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Sapiens Pseudolikelihood") + ylab("Ablang Pseudolikelihood") +
  ggtitle("Mouse")

ggarrange(a1, a2, ncol = 2, nrow = 1, common.legend = FALSE, legend = "right")
ggarrange(b1, b2, ncol = 2, nrow = 1, common.legend = FALSE, legend = "right")

#clonal expansion

#normalize by sample size
vdj_ova %>% group_by(sample_id) %>% summarise(sample_size = n()) -> sample_size
vdj_ova$samplesize <- sample_size$sample_size[match(vdj_ova$sample_id, sample_size$sample_id)]
vdj_ova %>% mutate(normalized_expansion = clonotype_frequency / samplesize) -> vdj_ova
vdj_human %>% group_by(sample_id) %>% summarise(sample_size = n()) -> sample_size
vdj_human$samplesize <- sample_size$sample_size[match(vdj_human$sample_id, sample_size$sample_id)]
vdj_human %>% mutate(normalized_expansion = clonotype_frequency / samplesize) -> vdj_human

vdj_all <- rbind(vdj_ova, vdj_human)

#esm
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_esm_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_esm_full_VDJ)$estimate
ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_esm_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Clonal Expansion") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#protbert
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_protbert_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_protbert_full_VDJ)$estimate
ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Clonal Expansion") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#sapiens
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Clonal Expansion") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(vdj_human$normalized_expansion, vdj_human$IGH_evo_likelihood_ablang_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$normalized_expansion, vdj_ova$IGH_evo_likelihood_ablang_full_VDJ)$estimate
ggplot(vdj_all, aes(x=normalized_expansion, y=IGH_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  scale_x_continuous(trans='log10') +
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("Clonal Expansion") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))


# ggplot(vdj_human, aes(x=SHM_count, y=normalized_expansion, color=v_gene_family)) +
#   geom_point(size = 1) +
#   scale_color_manual(name = "V-gene family", values = vgene_colors,
#                      breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
#                                 "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme_steropodon() +
#   theme(text = element_text(size = 18)) +
#   xlab("SHM count") + ylab("Clonal Expansion") +
#   ggtitle("Human")
# 
# ggplot(vdj_ova, aes(x=SHM_count, y=normalized_expansion, color=v_gene_family)) +
#   geom_point(size = 1) +
#   scale_color_manual(name = "V-gene family", values = vgene_colors,
#                      breaks = c("IGHV1", "IGHV2", "IGHV3","IGHV4","IGHV5","IGHV6", "IGHV7","IGHV8",
#                                 "IGHV9", "IGHV10", "IGHV11", "IGHV12", "IGHV13", "IGHV14", "IGHV15")) +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme_steropodon() +
#   theme(text = element_text(size = 18)) +
#   xlab("SHM count") + ylab("Clonal Expansion") +
#   ggtitle("Mouse")

#SHM
#esm
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_esm_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_esm_full_VDJ)$estimate
ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_esm_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("SHM count") + ylab("ESM-1b Pseudolikelihood") +
  ggtitle(paste0("R\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#protbert
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_protbert_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_protbert_full_VDJ)$estimate
ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("SHM count") + ylab("ProtBERT Pseudolikelihood") +
  ggtitle(paste0("ProtBERT\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#Sapiens
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ)$estimate
ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("SHM count") + ylab("Sapiens Pseudolikelihood") +
  ggtitle(paste0("Sapiens\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))

#ablang
cor_human <- cor.test(vdj_human$SHM_count, vdj_human$IGH_evo_likelihood_ablang_full_VDJ)$estimate
cor_mouse <- cor.test(vdj_ova$SHM_count, vdj_ova$IGH_evo_likelihood_ablang_full_VDJ)$estimate
ggplot(vdj_all, aes(x=SHM_count, y=IGH_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
  geom_point(size =1) +
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
  geom_smooth(method = "lm", se=F) +
  theme_steropodon() +
  theme(text = element_text(size = 18)) +
  xlab("SHM count") + ylab("Ablang Pseudolikelihood") +
  ggtitle(paste0("Ablang\nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mouse, digits = 2)))



#HC LC correlation
source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")
#Ablang
vdj_human$LC_evo_likelihood_ablang_full_VDJ <- vdj_human$IGL_evo_likelihood_ablang_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_ablang_full_VDJ)),"LC_evo_likelihood_ablang_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_ablang_full_VDJ)),"IGK_evo_likelihood_ablang_full_VDJ"]

vdj_ova$LC_evo_likelihood_ablang_full_VDJ <- vdj_ova$IGL_evo_likelihood_ablang_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_ablang_full_VDJ)),"LC_evo_likelihood_ablang_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_ablang_full_VDJ)),"IGK_evo_likelihood_ablang_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_ablang_full_VDJ, vdj_human$LC_evo_likelihood_ablang_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_ablang_full_VDJ, vdj_ova$LC_evo_likelihood_ablang_full_VDJ)$estimate
ggplot(vdj_all, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=LC_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  geom_smooth(method="lm", se = F) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("Ablang \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))

#Sapiens
vdj_human$LC_evo_likelihood_sapiens_full_VDJ <- vdj_human$IGL_evo_likelihood_sapiens_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_sapiens_full_VDJ)),"LC_evo_likelihood_sapiens_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_sapiens_full_VDJ)),"IGK_evo_likelihood_sapiens_full_VDJ"]

vdj_ova$LC_evo_likelihood_sapiens_full_VDJ <- vdj_ova$IGL_evo_likelihood_sapiens_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)),"LC_evo_likelihood_sapiens_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)),"IGK_evo_likelihood_sapiens_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_sapiens_full_VDJ, vdj_human$LC_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ, vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)$estimate

ggplot(vdj_all, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=LC_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  geom_smooth(method="lm", se = F) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("Sapiens \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))


#esm
vdj_human$LC_evo_likelihood_esm_full_VDJ <- vdj_human$IGL_evo_likelihood_esm_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_esm_full_VDJ)),"LC_evo_likelihood_esm_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_esm_full_VDJ)),"IGK_evo_likelihood_esm_full_VDJ"]

vdj_ova$LC_evo_likelihood_esm_full_VDJ <- vdj_ova$IGL_evo_likelihood_esm_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_esm_full_VDJ)),"LC_evo_likelihood_esm_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_esm_full_VDJ)),"IGK_evo_likelihood_esm_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_esm_full_VDJ, vdj_human$LC_evo_likelihood_esm_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_esm_full_VDJ, vdj_ova$LC_evo_likelihood_esm_full_VDJ)$estimate

ggplot(vdj_all, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=LC_evo_likelihood_esm_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  geom_smooth(method="lm", se = F) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("ESM-1b \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))
  

#protbert
vdj_human$LC_evo_likelihood_protbert_full_VDJ <- vdj_human$IGL_evo_likelihood_protbert_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_protbert_full_VDJ)),"LC_evo_likelihood_protbert_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_protbert_full_VDJ)),"IGK_evo_likelihood_protbert_full_VDJ"]

vdj_ova$LC_evo_likelihood_protbert_full_VDJ <- vdj_ova$IGL_evo_likelihood_protbert_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_protbert_full_VDJ)),"LC_evo_likelihood_protbert_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_protbert_full_VDJ)),"IGK_evo_likelihood_protbert_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_protbert_full_VDJ, vdj_human$LC_evo_likelihood_protbert_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_protbert_full_VDJ, vdj_ova$LC_evo_likelihood_protbert_full_VDJ)$estimate

ggplot(vdj_all, aes(x=IGH_evo_likelihood_protbert_full_VDJ, y=LC_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
  geom_point(size = 1) +
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
  geom_smooth(method="lm", se = F) +
  theme_steropodon() +
  theme(text = element_text(size = 15)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("ProtBERT \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))
