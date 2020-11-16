library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

load("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/R-objects/DNA_Methylation_Groups_new.Rds")

melt.df <- melt(df, id.vars = c("id","symbol","Gene.in.bulk", "Cont.vs.Diff", "Diff.vs.Reversible", "Cont.vs.Reversible", "Group", "in.bulk", "Regulation.bulk", "Symbol_ident"), variable.name = "Treatment", value.name = "Mean_Methylation")

GOI <- c("PTPRC", "COL4A1", "CNN1", "TAGLN", "NOS3", "KDR")
melt.df <- na.omit(melt.df)
melt.df <- filter(melt.df, Treatment %in% c("mean.mean.Control", "mean.mean.Differentiated", "mean.mean.Reversible"))
melt.df$Mean_Methylation <- as.numeric(melt.df$Mean_Methylation)
X <- melt.df %>% filter(Gene.in.bulk %in% GOI)

X$Gene.in.bulk <- factor(X$Gene.in.bulk, levels = GOI)

ggplot(X, aes(x=Gene.in.bulk, y = Mean_Methylation, fill = Treatment))+
  geom_col(color = "black", position = "dodge")+
  scale_fill_manual(values = c("grey80", "grey50", "grey30"))+
  scale_y_continuous(breaks = c(0, 0.3, 0.6))+
  labs(y = "Mean Methylation across sites", x = "Gene")+
  theme(legend.position = "bottom")
