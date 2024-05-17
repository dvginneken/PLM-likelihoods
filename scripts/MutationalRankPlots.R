#Mutational rank plots
library(ggpubr)

#Mutations
mut_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_table_OVA_V7_MP.csv")
mut_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_table_horns2020a__VDJ_RAW_MP.csv")
mut_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_table_Bruhn_MP.csv")
mut_bieberich <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/mutational_rank_table_Bieberich_MP.csv")

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

mut_bieberich$sample <- case_match(mut_bieberich$sample,
                                        "S2" ~ "Individual3",
                                        "S4" ~ "Individual4",
                                        "S7" ~ "Individual5")

mut_human <- rbind(mut_horns, mut_bruhn, mut_bieberich)
mut_all <- rbind(mut_ova, mut_human)

#separate on model
mut_all_esm <- mut_all[mut_all$model == "esm",]
mut_all_protbert <- mut_all[mut_all$model == "protbert",]
mut_all_sapiens <- mut_all[mut_all$model == "sapiens",]
mut_all_ablang <- mut_all[mut_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a <- ggplot(mut_all_esm, aes(mean_sub_rank)) +
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  ggtitle("ESM-1b")

b <- ggplot(mut_all_protbert, aes(mean_sub_rank)) +
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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


#Original residue
original_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/mutational_rank_reversed_table_OVA_V7_MP.csv")
original_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/mutational_rank_reversed_table_horns2020a__VDJ_RAW_MP.csv")
original_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/mutational_rank_reversed_table_Bruhn_MP.csv")
original_bieberich <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/mutational_rank_reversed_table_Bieberich_MP.csv")

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

original_bieberich$sample <- case_match(original_bieberich$sample,
                                   "S2" ~ "Individual3",
                                   "S4" ~ "Individual4",
                                   "S7" ~ "Individual5")

original_human <- rbind(original_horns, original_bruhn, original_bieberich)
original_all <- rbind(original_ova, original_human)

#separate on model
original_all_esm <- original_all[original_all$model == "esm",]
original_all_protbert <- original_all[original_all$model == "protbert",]
original_all_sapiens <- original_all[original_all$model == "sapiens",]
original_all_ablang <- original_all[original_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a <- ggplot(original_all_esm, aes(mean_original_residue_rank)) +
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  ggtitle("ESM-1b")

b <- ggplot(original_all_protbert, aes(mean_original_residue_rank)) +
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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

table(original_all_esm$mean_original_residue_rank)
471/nrow(original_all_esm)

#Conserved (unmutated) residues
conserved_ova <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/conserved_rank_table_OVA_V7_MP.csv")
conserved_horns <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/conserved_rank_table_horns2020a__VDJ_RAW_MP.csv")
conserved_bruhn <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/conserved_rank_table_Bruhn_MP.csv")
conserved_bieberich <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/conserved_rank_table_Bieberich_MP.csv")

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

conserved_bieberich$sample <- case_match(conserved_bieberich$sample,
                                        "S2" ~ "Individual3",
                                        "S4" ~ "Individual4",
                                        "S7" ~ "Individual5")

conserved_human <- rbind(conserved_horns, conserved_bruhn, conserved_bieberich)
conserved_all <- rbind(conserved_ova, conserved_human)

#separate on model
conserved_all_esm <- conserved_all[conserved_all$model == "esm",]
conserved_all_protbert <- conserved_all[conserved_all$model == "protbert",]
conserved_all_sapiens <- conserved_all[conserved_all$model == "sapiens",]
conserved_all_ablang <- conserved_all[conserved_all$model == "ablang",]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")

a <- ggplot(conserved_all_esm, aes(mean_residue_rank)) +
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  ggtitle("ESM-1b")

b <- ggplot(conserved_all_protbert, aes(mean_residue_rank)) +
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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
  geom_histogram(color = "white", fill = "black", binwidth = 1) +
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

mean(mut_all_esm$mean_sub_rank)
mean(original_all_esm$mean_original_residue_rank)
mean(conserved_all_esm$mean_residue_rank)
