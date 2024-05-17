#Binder Properties

#Clean the data
df <- read.csv("OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/Binder_properties.csv", sep = "\t")
df <- df[df$Bind..ELISA.signal.0.2. != "", ]
df <- df[df$to.remove != "yes",]
df <- df[,1:28]
#Remove lower bar from sequences
df <- df %>% mutate(HC_VDJ_AA = gsub("_", "", HC_VDJ_AA))
write.csv(df, file = "OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/Binder_properties_cleaned.csv",
          row.names = F)
