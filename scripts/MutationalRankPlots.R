#Mutational rank plots

library(ggpubr)
library(dplyr)
library(rstatix)

#Mutations
mut_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_table_OVA_V7_default.csv")
mut_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_table_horns2020a__VDJ_RAW_default.csv")
mut_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_table_Bruhn_default.csv")
mut_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/mutational_rank_table_Kim_default.csv")

#change sample names
mut_ova$sample <- case_match(mut_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")

mut_horns$sample <- case_match(mut_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
mut_horns <- mut_horns[mut_horns$sample == "Individual1",]

mut_bruhn$sample <- "Individual2"

mut_kim$sample <- case_match(mut_kim$sample,
                                        "SRR17729703" ~ "Individual3",
                                        "SRR17729692" ~ "Individual4",
                                        "SRR17729726" ~ "Individual5")
mut_kim <- mut_kim[mut_kim$sample %in% paste0("Individual",2:5),]

mut_human <- rbind(mut_horns, mut_bruhn, mut_kim)
mut_all <- rbind(mut_ova, mut_human)

mut_all1 <- mut_all[mut_all$n_subs == 1,]
mut_all_more <- mut_all[mut_all$n_subs > 1,]

#separate on model
mut_all_esm <- mut_all[mut_all$model == "esm",]
mut_all_protbert <- mut_all[mut_all$model == "protbert",]
mut_all_sapiens <- mut_all[mut_all$model == "sapiens",]
mut_all_ablang <- mut_all[mut_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a <- ggplot(mut_all_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Substitution")
  #ggtitle("ESM-1b")

b <- ggplot(mut_all_protbert, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("ProtBERT")

c <- ggplot(mut_all_sapiens, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Sapiens")
  
d <- ggplot(mut_all_ablang, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Ablang")

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")

#separated on nr mutations 
mut_all1 <- mut_all[mut_all$n_subs == 1,]
mut_all1 <- mut_all1[!(is.na(mut_all1$mean_sub_rank)),]
mut_all_more <- mut_all[mut_all$n_subs > 1,]
mut_all_more <- mut_all_more[!(is.na(mut_all_more$mean_sub_rank)),]

mut_all1_esm <- mut_all1[mut_all1$model == "esm",]
mut_all_more_esm <- mut_all_more[mut_all_more$model == "esm",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a <- ggplot(mut_all1_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Substitution")

b3 <- ggplot(mut_all_more_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlab("Substitution Rank") + ylab("Number of edges") +
  ggtitle("Multiple mutations")

#Figure main
ggarrange(a, a2, b3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

#Average rank plot
mut_all_esm$group <- "all_edges"
mut_all1_esm$group <- "single_mutations"
mut_all_more_esm$group <- "multiple_mutations"
mut_all_grouped <- rbind(mut_all_esm, mut_all1_esm, mut_all_more_esm)

stat.test <- mut_all_grouped %>%
  t_test(mean_sub_rank ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

ggplot(mut_all_grouped, aes(y = as.numeric(mean_sub_rank), x = group)) +
  geom_boxplot() +
  theme_steropodon() +
  theme(text = element_text(size = 14)) +
  xlab("Mutation type") + ylab("Substitution rank") +
  ggtitle("ESM-1b") +
  geom_signif(comparisons=list(c("single_mutations","multiple_mutations")),map_signif_level = TRUE,annotations="***")



#Original residue
original_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_default.csv")
original_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_reversed_table_horns2020a__VDJ_RAW_default.csv")
original_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_default.csv")
original_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/mutational_rank_reversed_table_Kim_default.csv")

#change sample names
original_ova$sample <- case_match(original_ova$sample,
                             "S1" ~ "Mouse1",
                             "S2" ~ "Mouse2",
                             "S3" ~ "Mouse3",
                             "S4" ~ "Mouse4",
                             "S5" ~ "Mouse5")

original_horns$sample <- case_match(original_horns$sample,
                               "Influenza.vac.11.12.human.S1" ~ "Individual1",
                               "Influenza.vac.11.12.human.S2" ~ "Human2",
                               "Influenza.vac.11.12.human.S3" ~ "Human3",
                               "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
original_horns <- original_horns[original_horns$sample == "Individual1",]

original_bruhn$sample <- "Individual2"

original_kim$sample <- case_match(original_kim$sample,
                             "SRR17729703" ~ "Individual3",
                             "SRR17729692" ~ "Individual4",
                             "SRR17729726" ~ "Individual5")
original_kim <- original_kim[original_kim$sample %in% paste0("Individual",2:5),]


original_human <- rbind(original_horns, original_bruhn, original_kim)
original_all <- rbind(original_ova, original_human)

#separate on model
original_all_esm <- original_all[original_all$model == "esm",]
original_all_protbert <- original_all[original_all$model == "protbert",]
original_all_sapiens <- original_all[original_all$model == "sapiens",]
original_all_ablang <- original_all[original_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a2 <- ggplot(original_all_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Original Residue Rank") + ylab("Number of edges") +
  ggtitle("Original")

b <- ggplot(original_all_protbert, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Original Residue Rank") + ylab("Number of edges") +
  ggtitle("ProtBERT")

c <- ggplot(original_all_sapiens, aes(mean_original_residue_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Original Residue Rank") + ylab("Number of edges") +
  ggtitle("Sapiens")

d <- ggplot(original_all_ablang, aes(mean_original_residue_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Original Residue Rank") + ylab("Number of edges") +
  ggtitle("Ablang")

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")

#separated on nr mutations 
original_all1 <- original_all[original_all$n_subs == 1,]
original_all1 <- original_all1[!(is.na(original_all1$mean_sub_rank)),]
original_all_more <- original_all[original_all$n_subs > 1,]
original_all_more <- original_all_more[!(is.na(original_all_more$mean_sub_rank)),]

original_all1_esm <- original_all1[original_all1$model == "esm",]
original_all_more_esm <- original_all_more[original_all_more$model == "esm",]

a1 <- ggplot(original_all1_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Original Residue Rank") + ylab("Number of edges") +
  ggtitle("Original")

b <- ggplot(original_all_more_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Original Residue Rank") + ylab("Number of edges") +
  ggtitle("Multiple mutations ESM-1b")

ggarrange(a, b, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
ggarrange(a, a1, b, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")

a <- ggplot(original_all1_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original residue Likelihood") + ylab("Number of edges") +
  ggtitle("Single mutation ESM-1b")

b <- ggplot(original_all_more_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original residue Likelihood") + ylab("Number of edges") +
  ggtitle("Multiple mutations ESM-1b")

ggarrange(a, b, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")

#Conserved (unmutated) residues
conserved_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_default.csv")
conserved_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/conserved_rank_table_horns2020a__VDJ_RAW_default.csv")
conserved_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/conserved_rank_table_Bruhn_default.csv")
conserved_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/conserved_rank_table_Kim_default.csv")

#change sample names
conserved_ova$sample <- case_match(conserved_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")

conserved_horns$sample <- case_match(conserved_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
conserved_horns <- conserved_horns[conserved_horns$sample == "Individual1",]

conserved_bruhn$sample <- "Individual2"

conserved_kim$sample <- case_match(conserved_kim$sample,
                                  "SRR17729703" ~ "Individual3",
                                  "SRR17729692" ~ "Individual4",
                                  "SRR17729726" ~ "Individual5")
conserved_kim <- conserved_kim[conserved_kim$sample %in% paste0("Individual",2:5),]

conserved_human <- rbind(conserved_horns, conserved_bruhn, conserved_kim)
conserved_all <- rbind(conserved_ova, conserved_human)

#separate on model
conserved_all_esm <- conserved_all[conserved_all$model == "esm",]
conserved_all_protbert <- conserved_all[conserved_all$model == "protbert",]
conserved_all_sapiens <- conserved_all[conserved_all$model == "sapiens",]
conserved_all_ablang <- conserved_all[conserved_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a3 <- ggplot(conserved_all_esm, aes(mean_sub_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Conserved Residue Rank") + ylab("Number of edges") +
  ggtitle("Conserved")

b <- ggplot(conserved_all_protbert, aes(mean_residue_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Conserved Residue Rank") + ylab("Number of edges") +
  ggtitle("ProtBERT")

c <- ggplot(conserved_all_sapiens, aes(mean_residue_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Conserved Residue Rank") + ylab("Number of edges") +
  ggtitle("Sapiens")

d <- ggplot(conserved_all_ablang, aes(mean_residue_rank)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,20) +
  xlab("Conserved Residue Rank") + ylab("Number of edges") +
  ggtitle("Ablang")

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")

png(file = "~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/Figures/RankPlot.png")
ggarrange(p_mut, p_orig, p_cons, ncol = 1, nrow = 3)
dev.off()

mean(mut_all_esm$mean_sub_rank)
mean(original_all_esm$mean_original_residue_rank)
mean(conserved_all_esm$mean_residue_rank)

mean(mut_all_ablang$mean_sub_rank)
mean(original_all_ablang$mean_original_residue_rank)
mean(conserved_all_ablang$mean_residue_rank)

####Probability in stead of rank
#Mutations
mut_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_table_OVA_V7_default.csv")
mut_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_table_horns2020a__VDJ_RAW_default.csv")
mut_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_table_Bruhn_default.csv")
mut_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/mutational_rank_table_Kim_default.csv")

#change sample names
mut_ova$sample <- case_match(mut_ova$sample,
                             "S1" ~ "Mouse1",
                             "S2" ~ "Mouse2",
                             "S3" ~ "Mouse3",
                             "S4" ~ "Mouse4",
                             "S5" ~ "Mouse5")

mut_horns$sample <- case_match(mut_horns$sample,
                               "Influenza.vac.11.12.human.S1" ~ "Individual1",
                               "Influenza.vac.11.12.human.S2" ~ "Human2",
                               "Influenza.vac.11.12.human.S3" ~ "Human3",
                               "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
mut_horns <- mut_horns[mut_horns$sample == "Individual1",]

mut_bruhn$sample <- "Individual2"

mut_kim$sample <- case_match(mut_kim$sample,
                             "SRR17729703" ~ "Individual3",
                             "SRR17729692" ~ "Individual4",
                             "SRR17729726" ~ "Individual5")
mut_kim <- mut_kim[mut_kim$sample %in% paste0("Individual",2:5),]

mut_human <- rbind(mut_horns, mut_bruhn, mut_kim)
mut_all <- rbind(mut_ova, mut_human)

#separate on model
mut_all_esm <- mut_all[mut_all$model == "esm",]
mut_all_protbert <- mut_all[mut_all$model == "protbert",]
mut_all_sapiens <- mut_all[mut_all$model == "sapiens",]
mut_all_ablang <- mut_all[mut_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a <- ggplot(mut_all_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Substitution Likelihood") + ylab("Number of edges") +
  ggtitle("Substitution")

b <- ggplot(mut_all_protbert, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Substitution Likelihood") + ylab("Number of edges") +
  ggtitle("ProtBERT")

c <- ggplot(mut_all_sapiens, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Substitution Likelihood") + ylab("Number of edges") +
  ggtitle("Sapiens")

d <- ggplot(mut_all_ablang, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Substitution Likelihood") + ylab("Number of edges") +
  ggtitle("Ablang")

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")


#Original residue
original_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_default.csv")
original_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_reversed_table_horns2020a__VDJ_RAW_default.csv")
original_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_default.csv")
original_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/mutational_rank_reversed_table_Kim_default.csv")

#change sample names
original_ova$sample <- case_match(original_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")

original_horns$sample <- case_match(original_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
original_horns <- original_horns[original_horns$sample == "Individual1",]

original_bruhn$sample <- "Individual2"

original_kim$sample <- case_match(original_kim$sample,
                                  "SRR17729703" ~ "Individual3",
                                  "SRR17729692" ~ "Individual4",
                                  "SRR17729726" ~ "Individual5")
original_kim <- original_kim[original_kim$sample %in% paste0("Individual",2:5),]


original_human <- rbind(original_horns, original_bruhn, original_kim)
original_all <- rbind(original_ova, original_human)

#separate on model
original_all_esm <- original_all[original_all$model == "esm",]
original_all_protbert <- original_all[original_all$model == "protbert",]
original_all_sapiens <- original_all[original_all$model == "sapiens",]
original_all_ablang <- original_all[original_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a2 <- ggplot(original_all_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Original")

b <- ggplot(original_all_protbert, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("ProtBERT")

c <- ggplot(original_all_sapiens, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Sapiens")

d <- ggplot(original_all_ablang, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Ablang")

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")


#Conserved (unmutated) residues
conserved_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_default.csv")
conserved_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/conserved_rank_table_horns2020a__VDJ_RAW_default.csv")
conserved_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/conserved_rank_table_Bruhn_default.csv")
conserved_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/conserved_rank_table_Kim_default.csv")

#change sample names
conserved_ova$sample <- case_match(conserved_ova$sample,
                                   "S1" ~ "Mouse1",
                                   "S2" ~ "Mouse2",
                                   "S3" ~ "Mouse3",
                                   "S4" ~ "Mouse4",
                                   "S5" ~ "Mouse5")

conserved_horns$sample <- case_match(conserved_horns$sample,
                                     "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                     "Influenza.vac.11.12.human.S2" ~ "Human2",
                                     "Influenza.vac.11.12.human.S3" ~ "Human3",
                                     "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
conserved_horns <- conserved_horns[conserved_horns$sample == "Individual1",]

conserved_bruhn$sample <- "Individual2"

conserved_kim$sample <- case_match(conserved_kim$sample,
                                   "SRR17729703" ~ "Individual3",
                                   "SRR17729692" ~ "Individual4",
                                   "SRR17729726" ~ "Individual5")
conserved_kim <- conserved_kim[conserved_kim$sample %in% paste0("Individual",2:5),]

conserved_human <- rbind(conserved_horns, conserved_bruhn, conserved_kim)
conserved_all <- rbind(conserved_ova, conserved_human)

#separate on model
conserved_all_esm <- conserved_all[conserved_all$model == "esm",]
conserved_all_protbert <- conserved_all[conserved_all$model == "protbert",]
conserved_all_sapiens <- conserved_all[conserved_all$model == "sapiens",]
conserved_all_ablang <- conserved_all[conserved_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a3 <- ggplot(conserved_all_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Conserved")

b <- ggplot(conserved_all_protbert, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("ProtBERT")

c <- ggplot(conserved_all_sapiens, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Sapiens")

d <- ggplot(conserved_all_ablang, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Ablang")

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")

png(file = "~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/Figures/RankPlot.png")
ggarrange(p_mut, p_orig, p_cons, ncol = 1, nrow = 3)
dev.off()

mean(mut_all_esm$mean_sub_rank)
mean(original_all_esm$mean_original_residue_rank)
mean(conserved_all_esm$mean_residue_rank)

mean(mut_all_ablang$mean_sub_rank)
mean(original_all_ablang$mean_original_residue_rank)
mean(conserved_all_ablang$mean_residue_rank)


conserved_all$group <- "conserved"

original_all$group <- "mutating"
df <- rbind(conserved_all[,-8], original_all[-8])

stat.test <- df %>%
  t_test(mean_sub_rank ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

ggplot(df, aes(y = as.numeric(mean_sub_rank), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ESM-1b Likelihood rank") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations="****",
              color = "black",
              size = 1,
              textsize = 5)



###
#Fig C
ggarrange(a, a2, a3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
ggarrange(a, a2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
#Fig D
conserved_all$group <- "conserved"
original_all$group <- "mutating"
df <- rbind(conserved_all[,-8], original_all[-8])
df <- df[df$model == "esm",]

stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ESM-1b Likelihood") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations="****",
              color = "black",
              size = 1,
              textsize = 5)



#Original residue
original_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_default.csv")
original_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_reversed_table_horns2020a__VDJ_RAW_default.csv")
original_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_default.csv")
original_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/mutational_rank_reversed_table_Kim_default.csv")
#change sample names
original_ova$sample <- case_match(original_ova$sample,
                                  "S1" ~ "Mouse1",
                                  "S2" ~ "Mouse2",
                                  "S3" ~ "Mouse3",
                                  "S4" ~ "Mouse4",
                                  "S5" ~ "Mouse5")
original_horns$sample <- case_match(original_horns$sample,
                                    "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                    "Influenza.vac.11.12.human.S2" ~ "Human2",
                                    "Influenza.vac.11.12.human.S3" ~ "Human3",
                                    "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
original_horns <- original_horns[original_horns$sample == "Individual1",]
original_bruhn$sample <- "Individual2"
original_kim$sample <- case_match(original_kim$sample,
                                  "SRR17729703" ~ "Individual3",
                                  "SRR17729692" ~ "Individual4",
                                  "SRR17729726" ~ "Individual5")
original_kim <- original_kim[original_kim$sample %in% paste0("Individual",2:5),]
original_human <- rbind(original_horns, original_bruhn, original_kim)
original_all <- rbind(original_ova, original_human)
#separated on nr mutations 
original_all1 <- original_all[original_all$n_subs == 1,]
original_all1 <- original_all1[!(is.na(original_all1$mean_sub_rank)),]
original_all_more <- original_all[original_all$n_subs > 1,]
original_all_more <- original_all_more[!(is.na(original_all_more$mean_sub_rank)),]

original_all1_esm <- original_all1[original_all1$model == "esm",]
original_all_more_esm <- original_all_more[original_all_more$model == "esm",]

a <- ggplot(original_all1_esm, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Original Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Original")

conserved_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_default.csv")
conserved_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/conserved_rank_table_horns2020a__VDJ_RAW_default.csv")
conserved_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/conserved_rank_table_Bruhn_default.csv")
conserved_kim <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/conserved_rank_table_Kim_default.csv")

#change sample names
conserved_ova$sample <- case_match(conserved_ova$sample,
                                   "S1" ~ "Mouse1",
                                   "S2" ~ "Mouse2",
                                   "S3" ~ "Mouse3",
                                   "S4" ~ "Mouse4",
                                   "S5" ~ "Mouse5")

conserved_horns$sample <- case_match(conserved_horns$sample,
                                     "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                     "Influenza.vac.11.12.human.S2" ~ "Human2",
                                     "Influenza.vac.11.12.human.S3" ~ "Human3",
                                     "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
conserved_horns <- conserved_horns[conserved_horns$sample == "Individual1",]

conserved_bruhn$sample <- "Individual2"

conserved_kim$sample <- case_match(conserved_kim$sample,
                                   "SRR17729703" ~ "Individual3",
                                   "SRR17729692" ~ "Individual4",
                                   "SRR17729726" ~ "Individual5")
conserved_kim <- conserved_kim[conserved_kim$sample %in% paste0("Individual",2:5),]

conserved_human <- rbind(conserved_horns, conserved_bruhn, conserved_kim)
conserved_all <- rbind(conserved_ova, conserved_human)

#separate on model
conserved_all_esm <- conserved_all[conserved_all$model == "esm",]

conserved_all1 <- conserved_all_esm[conserved_all_esm$n_subs == 1,]
conserved_all1 <- conserved_all1[!(is.na(conserved_all1$mean_sub_prob)),]

b <- ggplot(conserved_all1, aes(mean_sub_prob)) +
  #geom_histogram(color = "white", fill = "black", binwidth = 1) +
  geom_freqpoly(aes(colour = sample), binwidth = 0.1, linewidth = 2) +
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
  theme(text = element_text(size = 14)) +
  xlim(0,1) +
  xlab("Conserved Residue Likelihood") + ylab("Number of edges") +
  ggtitle("Conserved")

ggarrange(a, b, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")


conserved_all1$group <- "conserved"
original_all1_esm$group <- "mutating"
df <- rbind(conserved_all1[,-8], original_all1_esm[-8])
df <- df[df$model == "esm",]

stat.test <- df %>%
  t_test(mean_sub_prob ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

ggplot(df, aes(y = as.numeric(mean_sub_prob), x = group, color=group)) +
  geom_boxplot() +
  scale_color_manual(values = c("conserved" = "orange",
                                "mutating" = "dodgerblue3")) +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("ESM-1b Likelihood") +
  geom_signif(comparisons=list(c("conserved","mutating")),
              map_signif_level = TRUE,
              annotations="****",
              color = "black",
              size = 1,
              textsize = 5) +
  ggtitle("Single mutation edges")
