library(ggplot2) 
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

#Source correlation
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/SourceCorrelation.csv",
               header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/SourceCorrelation.csv",
                   header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/SourceCorrelation.csv",
                     header = TRUE, sep = ",")
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/SourceCorrelation.csv",
                     header = TRUE, sep = ",")

df_ova$sample <- case_match(df_ova$sample,
                       "S1" ~ "Mouse1",
                       "S2" ~ "Mouse2",
                       "S3" ~ "Mouse3",
                       "S4" ~ "Mouse4",
                       "S5" ~ "Mouse5")

df_horns$sample <- case_match(df_horns$sample,
                            "Influenza.vac.11.12.human.S1" ~ "Individual1",
                            "Influenza.vac.11.12.human.S2" ~ "Human2",
                            "Influenza.vac.11.12.human.S3" ~ "Human3",
                            "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
df_horns <- df_horns[df_horns$sample == "Individual1",]

df_bruhn$sample <- "Individual2"

df_kim$sample <- case_match(df_kim$sample,
                             "SRR17729703" ~ "Individual3",
                             "SRR17729692" ~ "Individual4",
                             "SRR17729726" ~ "Individual5")
df_kim <- df_kim[df_kim$sample %in% paste0("Individual",2:5),]

df <- rbind(df_ova, df_horns, df_bruhn, df_kim)

df <- pivot_longer(df, cols = 1:3, names_to = "source", values_to = "correlation")
df$source<- case_match(df$source,
                            "full_VDJ__CDR3_only" ~ "VDJ\nCDR3",
                            "full_VDJ__CDR3_from_VDJ" ~ "VDJ\nCDR3-VDJ",
                            "CDR3_only__CDR3_from_VDJ" ~ "CDR3\nCDR3-VDJ")
df_main <- df[df$model %in% c("ESM-1b", "Ablang"),]

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")
ggplot(df_main, aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 4) +
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
  theme(text = element_text(size = 18)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model)
ggsave("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/Figures/Figure2/SourceCorrelation.png", width = 8, height = 6)

#supplementals
ggplot(df, aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~model, nrow = 1, ncol = 4)
ggsave("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/Figures/Figure2/SourceCorrelation_sup1.png", width = 12, height = 4)


#HC LC correlation
#Read VDJ dataframes
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/VDJ_PLL_OVA_V7.RData")
vdj_ova <- vdj
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/VDJ_PLL_horns2020a__VDJ_RAW.RData")
vdj_horns <- vdj
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/VDJ_PLL_Bruhn.RData")
vdj_bruhn <- vdj
rm(vdj)
load("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/VDJ_PLL_Kim.RData")
vdj_kim <- vdj
rm(vdj)

#Change sample names
vdj_ova$sample_id <- case_match(vdj_ova$sample_id,
                                "S1" ~ "Mouse1",
                                "S2" ~ "Mouse2",
                                "S3" ~ "Mouse3",
                                "S4" ~ "Mouse4",
                                "S5" ~ "Mouse5")
vdj_ova$sample <- vdj_ova$sample_id

vdj_horns$sample_id <- case_match(vdj_horns$sample_id,
                                  "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                  "Influenza.vac.11.12.human.S2" ~ "Human2",
                                  "Influenza.vac.11.12.human.S3" ~ "Human3",
                                  "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
vdj_horns <- vdj_horns[vdj_horns$sample_id == "Individual1",]
vdj_horns$sample <- vdj_horns$sample_id

vdj_bruhn$sample_id <- "Individual2"
vdj_bruhn$sample <- vdj_bruhn$sample_id

vdj_kim <- vdj_kim[vdj_kim$sample_id %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
vdj_kim$sample_id <- case_match(vdj_kim$sample_id,
                                "SRR17729703" ~ "Individual3",
                                "SRR17729692" ~ "Individual4",
                                "SRR17729726" ~ "Individual5")
vdj_kim$sample<- vdj_kim$sample_id


#combine human samples
vdj_human <- rbind(vdj_horns, vdj_bruhn, vdj_kim)
vdj_all <- rbind(vdj_ova, vdj_human)



#Ablang
vdj_human$LC_evo_likelihood_ablang_full_VDJ <- vdj_human$IGL_evo_likelihood_ablang_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_ablang_full_VDJ)),"LC_evo_likelihood_ablang_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_ablang_full_VDJ)),"IGK_evo_likelihood_ablang_full_VDJ"]

vdj_ova$LC_evo_likelihood_ablang_full_VDJ <- vdj_ova$IGL_evo_likelihood_ablang_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_ablang_full_VDJ)),"LC_evo_likelihood_ablang_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_ablang_full_VDJ)),"IGK_evo_likelihood_ablang_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_ablang_full_VDJ, vdj_human$LC_evo_likelihood_ablang_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_ablang_full_VDJ, vdj_ova$LC_evo_likelihood_ablang_full_VDJ)$estimate
a <-ggplot(vdj_all, aes(x=IGH_evo_likelihood_ablang_full_VDJ, y=LC_evo_likelihood_ablang_full_VDJ, color=sample_id)) +
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
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("Ablang \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))

#esm
vdj_human$LC_evo_likelihood_esm_full_VDJ <- vdj_human$IGL_evo_likelihood_esm_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_esm_full_VDJ)),"LC_evo_likelihood_esm_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_esm_full_VDJ)),"IGK_evo_likelihood_esm_full_VDJ"]

vdj_ova$LC_evo_likelihood_esm_full_VDJ <- vdj_ova$IGL_evo_likelihood_esm_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_esm_full_VDJ)),"LC_evo_likelihood_esm_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_esm_full_VDJ)),"IGK_evo_likelihood_esm_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_esm_full_VDJ, vdj_human$LC_evo_likelihood_esm_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_esm_full_VDJ, vdj_ova$LC_evo_likelihood_esm_full_VDJ)$estimate

b <- ggplot(vdj_all, aes(x=IGH_evo_likelihood_esm_full_VDJ, y=LC_evo_likelihood_esm_full_VDJ, color=sample_id)) +
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
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("ESM-1b \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))

ggarrange(a, b, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")

#Sapiens
vdj_human$LC_evo_likelihood_sapiens_full_VDJ <- vdj_human$IGL_evo_likelihood_sapiens_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_sapiens_full_VDJ)),"LC_evo_likelihood_sapiens_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_sapiens_full_VDJ)),"IGK_evo_likelihood_sapiens_full_VDJ"]

vdj_ova$LC_evo_likelihood_sapiens_full_VDJ <- vdj_ova$IGL_evo_likelihood_sapiens_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)),"LC_evo_likelihood_sapiens_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)),"IGK_evo_likelihood_sapiens_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_sapiens_full_VDJ, vdj_human$LC_evo_likelihood_sapiens_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_sapiens_full_VDJ, vdj_ova$LC_evo_likelihood_sapiens_full_VDJ)$estimate

c <- ggplot(vdj_all, aes(x=IGH_evo_likelihood_sapiens_full_VDJ, y=LC_evo_likelihood_sapiens_full_VDJ, color=sample_id)) +
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
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("Sapiens \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))



#protbert
vdj_human$LC_evo_likelihood_protbert_full_VDJ <- vdj_human$IGL_evo_likelihood_protbert_full_VDJ
vdj_human[which(is.na(vdj_human$LC_evo_likelihood_protbert_full_VDJ)),"LC_evo_likelihood_protbert_full_VDJ"] <- vdj_human[which(is.na(vdj_human$LC_evo_likelihood_protbert_full_VDJ)),"IGK_evo_likelihood_protbert_full_VDJ"]

vdj_ova$LC_evo_likelihood_protbert_full_VDJ <- vdj_ova$IGL_evo_likelihood_protbert_full_VDJ
vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_protbert_full_VDJ)),"LC_evo_likelihood_protbert_full_VDJ"] <- vdj_ova[which(is.na(vdj_ova$LC_evo_likelihood_protbert_full_VDJ)),"IGK_evo_likelihood_protbert_full_VDJ"]

vdj_all <- rbind(vdj_ova, vdj_human)

cor_human <- cor.test(vdj_human$IGH_evo_likelihood_protbert_full_VDJ, vdj_human$LC_evo_likelihood_protbert_full_VDJ)$estimate
cor_mice <- cor.test(vdj_ova$IGH_evo_likelihood_protbert_full_VDJ, vdj_ova$LC_evo_likelihood_protbert_full_VDJ)$estimate

d <- ggplot(vdj_all, aes(x=IGH_evo_likelihood_protbert_full_VDJ, y=LC_evo_likelihood_protbert_full_VDJ, color=sample_id)) +
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
  theme(text = element_text(size = 12)) +
  xlab("Heavy Chain Pseudolikelihood") + ylab("Light Chain Pseudolikelihood") +
  ggtitle(paste0("ProtBERT \nR\u00b2 Individuals = ", round(cor_human, digits = 2), ", R\u00b2 Mice = ", round(cor_mice, digits = 2)))

ggarrange(a, b, c, d, ncol = 4, nrow = 1, common.legend = TRUE, legend = "right")

#HC and LC separate
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/SourceCorrelation_chains.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")

df_ova$sample <- case_match(df_ova$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")

df_horns$sample <- case_match(df_horns$sample,
                              "Influenza.vac.11.12.human.S1" ~ "Individual1",
                              "Influenza.vac.11.12.human.S2" ~ "Human2",
                              "Influenza.vac.11.12.human.S3" ~ "Human3",
                              "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
df_horns <- df_horns[df_horns$sample == "Individual1",]

df_bruhn$sample <- "Individual2"

df_kim$sample <- case_match(df_kim$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")
df_kim <- df_kim[df_kim$sample %in% paste0("Individual",2:5),]

df <- rbind(df_ova, df_horns, df_bruhn, df_kim)
df <- pivot_longer(df, cols = 1:3, names_to = "source", values_to = "correlation")
df$source<- case_match(df$source,
                       "full_VDJ__CDR3_only" ~ "VDJ\nCDR3",
                       "full_VDJ__CDR3_from_VDJ" ~ "VDJ\nCDR3-VDJ",
                       "CDR3_only__CDR3_from_VDJ" ~ "CDR3\nCDR3-VDJ")

df_sup<- df[df$model %in% c("ESM-1b", "Ablang"),]

a <- ggplot(df[df$model == "ESM-1b",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("ESM-1b")

b <- ggplot(df[df$model == "ProtBERT",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("ProtBERT")

c <- ggplot(df[df$model == "Ablang",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("Ablang")

d <- ggplot(df[df$model == "Sapiens",], aes(x=source, y=correlation, col=factor(sample))) + 
  geom_point(size = 2) +
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
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("Sapiens")

ggarrange(a, b, c, d, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")

#PLM correlation
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/PLMCorrelation.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/PLMCorrelation.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/PLMCorrelation.csv",
                     header = TRUE, sep = ",")
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/PLMCorrelation.csv",
                   header = TRUE, sep = ",")

df_ova$sample <- case_match(df_ova$sample,
                            "S1" ~ "Mouse1",
                            "S2" ~ "Mouse2",
                            "S3" ~ "Mouse3",
                            "S4" ~ "Mouse4",
                            "S5" ~ "Mouse5")

df_horns$sample <- case_match(df_horns$sample,
                              "Influenza.vac.11.12.human.S1" ~ "Individual1",
                              "Influenza.vac.11.12.human.S2" ~ "Human2",
                              "Influenza.vac.11.12.human.S3" ~ "Human3",
                              "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
df_horns <- df_horns[df_horns$sample == "Individual1",]

df_bruhn$sample <- "Individual2"
df_kim <- df_kim[df_kim$sample %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
df_kim$sample <- case_match(df_kim$sample,
                            "SRR17729703" ~ "Individual3",
                            "SRR17729692" ~ "Individual4",
                            "SRR17729726" ~ "Individual5")

df <- rbind(df_ova, df_horns, df_bruhn, df_kim)


df <- pivot_longer(df, cols = 1:6, names_to = "PLM", values_to = "correlation")

ggplot(df, aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 4) +
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
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~source)


#PLM correlation
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/PLMCorrelation_chains.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/PLMCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/PLMCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/PLMCorrelation_chains.csv",
                     header = TRUE, sep = ",")

#Change sample names
df_ova$sample <- case_match(df_ova$sample,
                                "S1" ~ "Mouse1",
                                "S2" ~ "Mouse2",
                                "S3" ~ "Mouse3",
                                "S4" ~ "Mouse4",
                                "S5" ~ "Mouse5")

df_horns$sample <- case_match(df_horns$sample,
                                  "Influenza.vac.11.12.human.S1" ~ "Individual1",
                                  "Influenza.vac.11.12.human.S2" ~ "Human2",
                                  "Influenza.vac.11.12.human.S3" ~ "Human3",
                                  "Influenza.vac.11.12.human.S4" ~ "Human4")
#Only keep 1 replicate
df_horns <- df_horns[df_horns$sample == "Individual1",]

df_bruhn$sample <- "Individual2"

df_kim <- df_kim[df_kim$sample %in% c("SRR17729703", "SRR17729692", "SRR17729726"),]
df_kim$sample <- case_match(df_kim$sample,
                                "SRR17729703" ~ "Individual3",
                                "SRR17729692" ~ "Individual4",
                                "SRR17729726" ~ "Individual5")

df <- rbind(df_ova, df_horns, df_bruhn, df_kim)

df <- pivot_longer(df, cols = 1:6, names_to = "PLM", values_to = "correlation")


b <- ggplot(df[df$source == "CDR3 from VDJ",], aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  ylim(-0.25,1) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("CDR3 from VDJ")

a <- ggplot(df[df$source == "CDR3 only",], aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  ylim(-0.25,1) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("CDR3 only")

c <- ggplot(df[df$source == "Full VDJ",], aes(x=PLM, y=correlation, col=factor(sample))) + 
  geom_point(size = 3) +
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
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  ylim(-0.25,1) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain) + ggtitle("Full VDJ")

ggarrange(a, b, c, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
