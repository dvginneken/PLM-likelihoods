library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

df<-read.delim("~/OneDrive - UMC Utrecht/Documenten/Kim/evo_likelihood_esm.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)

stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(
  df, x = "ELISA", y = "evo_likelihood", 
  color = "ELISA", palette = c("#00AFBB", "#E7B800")
)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "ELISA", dodge = 0.8)
bxp <- bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + ylab("ESM-1b pseudolikelihood") + ggtitle("ESM-1b")


df<-read.delim("~/OneDrive - UMC Utrecht/Documenten/Kim/evo_likelihood_ablang.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)

stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp_ablang <- ggboxplot(
  df, x = "ELISA", y = "evo_likelihood", 
  color = "ELISA", palette = c("#00AFBB", "#E7B800")
)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "ELISA", dodge = 0.8)
bxp_ablang <- bxp_ablang + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + ylab("Ablang pseudolikelihood") + ggtitle("Ablang")

df<-read.delim("~/OneDrive - UMC Utrecht/Documenten/Kim/evo_likelihood_protbert.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)

stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp_protbert <- ggboxplot(
  df, x = "ELISA", y = "evo_likelihood", 
  color = "ELISA", palette = c("#00AFBB", "#E7B800")
)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "ELISA", dodge = 0.8)
bxp_protbert <- bxp_protbert + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + ylab("ProtBERT pseudolikelihood") + ggtitle("ProtBERT")

df<-read.delim("~/OneDrive - UMC Utrecht/Documenten/Kim/evo_likelihood_sapiens.csv",sep=",", header = T)
df$ELISA <- case_match(df$elisa,
                       "True" ~ "Binder",
                       "False" ~ "Non-binder")
df$ELISA <- as.factor(df$ELISA)

stat.test <- df %>%
  t_test(evo_likelihood ~ ELISA) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp_sapiens <- ggboxplot(
  df, x = "ELISA", y = "evo_likelihood", 
  color = "ELISA", palette = c("#00AFBB", "#E7B800")
)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "ELISA", dodge = 0.8)
bxp_sapiens <- bxp_sapiens + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) + ylab("Sapiens pseudolikelihood") + ggtitle("Sapiens")

ggarrange(bxp, bxp_protbert, bxp_sapiens, bxp_ablang, ncol = 4, nrow = 1, common.legend = TRUE, legend = "top")

