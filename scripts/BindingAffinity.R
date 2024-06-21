library("readxl")
library(dplyr)
library(ggplot2)

##OVA
#ESM
df<-read.delim("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/BinderNonBinder/evo_likelihood_esm.csv",sep=",", header = T)
df<-df[df$Bind..ELISA.signal.0.2. == "yes",]
df<-df[df$octet.affinity..nM. != "",]
df<-df[df$octet.affinity..nM. != "nd",]


##Kim
df_kim<-read.delim("~/OneDrive - UMC Utrecht/Documenten/Kim/WU368_kim_et_al_nature_2022_bcr_heavy.tsv",sep="\t", header = T, row.names = NULL)
affinity_kim <- read_xlsx("~/OneDrive - UMC Utrecht/Documenten/Kim/2022-12-05_ed_table_6.xlsx")
affinity_kim$sequence_id <- affinity_kim$h_id
affinity_kim$sequence_id %in% df_kim$sequence_id
new_df <- left_join(affinity_kim, df_kim, by = "sequence_id")

translate_DNA<- function(sequence){
  
  #Translate a nucleotide sequence into an amino acid sequence
  #Arguments:
  #- sequence: nucleotide sequence to be translated
  
  if (sequence == ""){
    return(NA)
  }
  
  #Genetic code
  genetic_code <- list(
    "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
    "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
    "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
    "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
    "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
    "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
    "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
    "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
    "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
    "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
    "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
    "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
    "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
    "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
    "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
    "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
  )
  #Split the sequence into codons
  codons <- strsplit(sequence, "(?<=.{3})", perl=TRUE)[[1]]
  
  #Translate the codons
  for (codon_id in 1:length(codons)){
    #Remove codons that are not complete
    if(nchar(codons[codon_id]) < 3){
      codons[codon_id] = ""
    }
    #Codons that contain "-" are replaced with "-"
    else if (grepl("-", codons[codon_id], fixed = TRUE)){
      codons[codon_id] = "-"
    }
    #Codons that contain "N" are replaced with "-"
    else if (grepl("N", codons[codon_id], fixed = TRUE)){
      codons[codon_id] = "-"
    }
    #Translate codons according to the genetic code
    else{
      codons[codon_id] = genetic_code[[codons[codon_id]]]
    }
  }
  
  #Paste the codons together
  sequence <- paste(codons, collapse="")
  
  #Return the sequence
  return(sequence)
}
new_df$sequence_alignment<-gsub("\\.","-",new_df$sequence_alignment)
new_df$full_sequence <- apply(new_df, 1, function(x) translate_DNA(x["sequence_alignment"]))
new_df$full_sequence <-gsub("-","",new_df$full_sequence)
write.csv(new_df, file = "~/OneDrive - UMC Utrecht/Documenten/Kim/Affinity_dataframe_HC.csv", row.names = FALSE)

#esm
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/Affinity/evo_likelihood_esm.csv", header = T)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=donor)) +
  geom_point() +
  scale_y_log10() +
  ylab("Affinity (Kd)") + xlab("ESM-1b HC pseudolikelihood") +
  ggtitle(paste0("Human Dataset\nR\u00b2 = ", round(cor, digits = 3))) +
  theme_minimal()
ggplot(df_kim, aes(x=tissue, y=evo_likelihood)) +
  geom_point()


#ablang
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/Affinity/evo_likelihood_ablang.csv", header = T)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=tissue)) +
  geom_point() +
  scale_y_log10() +
  ylab("Affinity (Kd)") + xlab("Ablang HC pseudolikelihood") +
  ggtitle(paste0("Human Dataset\nR\u00b2 = ", round(cor, digits = 3))) +
  theme_minimal()
ggplot(df_kim, aes(x=tissue, y=evo_likelihood)) +
  geom_point()

#protbert
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/Affinity/evo_likelihood_protbert.csv", header = T)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=tissue)) +
  geom_point() +
  scale_y_log10() +
  ylab("Affinity (Kd)") + xlab("ProtBERT HC pseudolikelihood") +
  ggtitle(paste0("Human Dataset\nR\u00b2 = ", round(cor, digits = 3))) +
  theme_minimal()
ggplot(df_kim, aes(x=tissue, y=evo_likelihood)) +
  geom_point()

#sapiens
df_kim <- read.csv("~/OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/Kim/Affinity/evo_likelihood_sapiens.csv", header = T)
cor <- cor.test(df_kim$evo_likelihood, df_kim$K_D_nM, method = "spearman")$estimate
ggplot(df_kim, aes(x=evo_likelihood, y=K_D_nM, color=tissue)) +
  geom_point() +
  scale_y_log10() +
  ylab("Affinity (Kd)") + xlab("Sapiens HC pseudolikelihood") +
  ggtitle(paste0("Human Dataset\nR\u00b2 = ", round(cor, digits = 3))) +
  theme_minimal()
ggplot(df_kim, aes(x=tissue, y=evo_likelihood)) +
  geom_point()
