library(ggplot2) 
library(tidyr)
library(dplyr)
library(RColorBrewer)


#Source correlation
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/SourceCorrelation.csv",
               header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/SourceCorrelation.csv",
                   header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/SourceCorrelation.csv",
                     header = TRUE, sep = ",")
df_bieberich <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/SourceCorrelation.csv",
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

df_bieberich$sample <- case_match(df_bieberich$sample,
                              "S2" ~ "Individual3",
                              "S4" ~ "Individual4",
                              "S7" ~ "Individual5")

df <- rbind(df_ova, df_horns, df_bruhn, df_bieberich)

df <- pivot_longer(df, cols = 1:3, names_to = "source", values_to = "correlation")

source("~/OneDrive - UMC Utrecht/Documenten/Steropodon_theme.R")
ggplot(df, aes(x=source, y=correlation, col=factor(sample), shape=factor(model))) + 
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
  scale_shape_manual(values = c("Ablang" = 15,
                                "Sapiens" = 19,
                                "ESM-1b" = 3,
                                "ProtBERT" = 4),
                     name = "PLM") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient")



#PLM correlation
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/PLMCorrelation.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/PLMCorrelation.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/PLMCorrelation.csv",
                     header = TRUE, sep = ",")
df_bieberich <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/PLMCorrelation.csv",
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

df_bieberich$sample <- case_match(df_bieberich$sample,
                                  "S2" ~ "Individual3",
                                  "S4" ~ "Individual4",
                                  "S7" ~ "Individual5")

df <- rbind(df_ova, df_horns, df_bruhn, df_bieberich)


df <- pivot_longer(df, cols = 1:6, names_to = "PLM", values_to = "correlation")

ggplot(df, aes(x=PLM, y=correlation, col=factor(sample), shape=factor(source))) + 
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
  scale_shape_manual(values = c("Full VDJ" = 15,
                                "CDR3 from VDJ" = 19,
                                "CDR3 only" = 17),
                     name = "Source") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient")

#HC and LC separate
df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/SourceCorrelation_chains.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/SourceCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bieberich <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/SourceCorrelation_chains.csv",
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

df_bieberich$sample <- case_match(df_bieberich$sample,
                                  "S2" ~ "Individual3",
                                  "S4" ~ "Individual4",
                                  "S7" ~ "Individual5")

df <- rbind(df_ova, df_horns, df_bruhn, df_bieberich)
df <- pivot_longer(df, cols = 1:3, names_to = "source", values_to = "correlation")

ggplot(df, aes(x=source, y=correlation, col=factor(sample), shape=factor(model))) + 
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
  scale_shape_manual(values = c("Ablang" = 15,
                                "Sapiens" = 19,
                                "ESM-1b" = 3,
                                "ProtBERT" = 4),
                     name = "PLM") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("Source Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain)




df_ova <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/PLMCorrelation_chains.csv",
                   header = TRUE, sep = ",")
df_horns <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Horns/PLMCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bruhn <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bruhn/PLMCorrelation_chains.csv",
                     header = TRUE, sep = ",")
df_bieberich <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Bieberich/PLMCorrelation_chains.csv",
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

df_bieberich$sample <- case_match(df_bieberich$sample,
                                  "S2" ~ "Individual3",
                                  "S4" ~ "Individual4",
                                  "S7" ~ "Individual5")

df <- rbind(df_ova, df_horns, df_bruhn, df_bieberich)

df <- pivot_longer(df, cols = 1:6, names_to = "PLM", values_to = "correlation")

ggplot(df, aes(x=PLM, y=correlation, col=factor(sample), shape=factor(source))) + 
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
  scale_shape_manual(values = c("Full VDJ" = 15,
                                "CDR3 from VDJ" = 19,
                                "CDR3 only" = 17),
                     name = "Source") +
  theme_steropodon() +
  theme(text = element_text(size = 18),
        axis.text.x = element_text(angle = 50, hjust = 1)) +
  xlab("PLM Comparison") + ylab("Correlation Coefficient") +
  facet_wrap(~chain)

