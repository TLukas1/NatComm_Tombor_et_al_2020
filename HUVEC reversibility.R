# This file contains all qPCR analysis from Tombor et al. 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyr)
library(scales)

################################## Figure 4 H and I ##############################################

#Significance calculated by Stefanie (additional file /media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/qPCR data Reversibility/Statistic reversibility.pzf)

# Reversible HUVEC experiment

# Data from Simone HUVEC assay with reversible (EndMT Reversibility n1-n3 d1-d7 an d3+FM.xlsx)
SEM <- function(x) sqrt(var(x)/length(x))


df.qPCR <- read.csv("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/qPCR data Reversibility/EndMT Reversibility_Summary.csv")

df.qPCR <- df.qPCR %>% group_by(Sample) %>% mutate(SM22.mean = mean(Fold_to_d1_FM_SM22), SM22.SEM = SEM(Fold_to_d1_FM_SM22),
                                                   CNN1.mean = mean(Fold_to_d1_FM_CNN1), CNN1.SEM = SEM(Fold_to_d1_FM_CNN1))


df.qPCR$Sample <- factor(df.qPCR$Sample, levels = c("FM d1", "FM d7", "DM+T d1", "DM+T d2", "DM+T d3", "DM+T d7", "DM+T + FM d7"))

df.qPCR <- na.omit(df.qPCR)


ggplot(filter(df.qPCR), aes(x = Sample, y = Fold_to_d1_FM_SM22))+
  stat_summary(geom = "col", fun = mean, color = "black", fill = "white")+
  geom_jitter(height = 0)+
  geom_errorbar(aes(ymax = SM22.mean + SM22.SEM, ymin = SM22.mean-SM22.SEM), width = .2)+
  theme_classic()+
  geom_vline(xintercept = 5.5, linetype = "dotted")+
  labs(y = "Fold to SM22", title = "qPCR reversibility SM22")

ggplot(filter(df.qPCR), aes(x = Sample, y = Fold_to_d1_FM_CNN1))+
  stat_summary(geom = "col", fun = mean, color = "black", fill = "white")+
  geom_jitter(height = 0)+
  geom_errorbar(aes(ymax = CNN1.mean + CNN1.SEM, ymin = CNN1.mean-CNN1.SEM), width = .2)+
  theme_classic()+
  geom_vline(xintercept = 5.5, linetype = "dotted")+
  labs(y = "Fold to CNN1", title = "qPCR reversibility CNN1")

######################## Figure Supplement 15d and e #########################################

Analysis <- read_csv("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/qPCR data Reversibility/Analysis_LT20_006.csv")

Desc <- read.table("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/qPCR data Reversibility/Sample_desc_LT_20_006.csv", sep = ";", header = T)
Desc <- Desc %>% replace_na(list(Supplement = ""))
Analysis$Sample <- factor(Analysis$Sample)

Analysis <- unique(Analysis)
Analysis <- Analysis[1:51, ]

Analysis$Sample <- factor(Analysis$Sample)
Desc$Sample <- factor(Desc$Sample)
Desc$Treatment <- factor(Desc$Treatment, levels = c("FM", "DM", "DM+FM"))
Desc$N <- factor(Desc$N)
Desc$Timepoint <- factor(Desc$Timepoint, levels = c("d0", "d3", "d7", "d3+d7"))
Desc$Supplement <- factor(Desc$Supplement, levels = c("TgfbII", "TgfbI + Il1b", "L-NAME", "TgfbII + L-NAME", "TgfbI+Il1b+L-NAME", ""))

df <- full_join(Analysis, Desc)

df <- df %>% group_by(Treatment, Timepoint, Supplement) %>% mutate(SM22_F_Mean = mean(SM22), SM22_F_SEM = SEM(SM22), CNN1_F_Mean = mean(CNN1), CNN1_F_SEM = SEM(CNN1),
                                                                   CDH5_F_Mean = mean(CDH5), CDH5_F_SEM = SEM(CDH5), F_Group = paste(Treatment, Timepoint, Supplement))

df$F_Group <- factor(df$F_Group, levels = c("FM d0 ", "FM d3 ", "FM d7 ", "DM d3 TgfbII", "DM d7 TgfbII", "DM+FM d3+d7 TgfbII", "DM d3 TgfbI + Il1b", "DM d7 TgfbI + Il1b", "DM+FM d3+d7 TgfbI + Il1b", "DM d3 L-NAME", "DM d7 L-NAME", "DM d3 TgfbII + L-NAME", "DM d7 TgfbII + L-NAME", "DM+FM d3+d7 TgfbII + L-NAME", "DM d3 TgfbI+Il1b+L-NAME", "DM d7 TgfbI+Il1b+L-NAME", "DM+FM d3+d7 TgfbI+Il1b+L-NAME")) 
ggplot(df, aes(x = F_Group, y= SM22, group = F_Group))+
  geom_col(aes(y = SM22_F_Mean, fill = F_Group), position = "dodge", color = "black", alpha = 0.4, width = .8)+
  geom_errorbar(aes(ymin = SM22_F_Mean-SM22_F_SEM, ymax = SM22_F_Mean+SM22_F_SEM), width = .1)+
  geom_jitter(height = 0, width = 0.2)+
  #stat_compare_means(label = "p.signif", method = "t.test", ref.group = "siRNA Control FM")+
  theme_minimal()+
  scale_x_discrete(name = "Group")+
  scale_fill_manual(values = c("#383131", "#383131", "#383131", "#ed3621", "#ed3621", "#faa79d", "#1c79eb", "#1c79eb", "#9bc2f2", "#46ebe8", "#46ebe8", "#ed3621", "#ed3621", "#faa79d", "#1c79eb", "#1c79eb", "#9bc2f2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank())+
  labs(title = "LT20/006 - HUVEC different conditions", 
       y= "2^ddCT normalized to d0 FM \n
       GAPDH and RPLP0")

ggplot(df, aes(x = F_Group, y= CNN1, group = F_Group))+
  geom_col(aes(y = CNN1_F_Mean, fill = F_Group), position = "dodge", color = "black", alpha = 0.4, width = .8)+
  geom_errorbar(aes(ymin = CNN1_F_Mean-CNN1_F_SEM, ymax = CNN1_F_Mean+CNN1_F_SEM), width = .1)+
  geom_jitter(height = 0, width = 0.2)+
  #stat_compare_means(label = "p.signif", method = "t.test", ref.group = "siRNA Control FM")+
  theme_minimal()+
  scale_x_discrete(name = "Group")+
  scale_fill_manual(values = c("#383131", "#383131", "#383131", "#ed3621", "#ed3621", "#faa79d", "#1c79eb", "#1c79eb", "#9bc2f2", "#46ebe8", "#46ebe8", "#ed3621", "#ed3621", "#faa79d", "#1c79eb", "#1c79eb", "#9bc2f2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank())+
  labs(title = "LT20/006 - HUVEC different conditions", 
       y= "2^ddCT normalized to d0 FM \n
       GAPDH and RPLP0")

ggplot(df, aes(x = F_Group, y= CDH5, group = F_Group))+
  geom_col(aes(y = CDH5_F_Mean, fill = F_Group), position = "dodge", color = "black", alpha = 0.4, width = .8)+
  geom_errorbar(aes(ymin = CDH5_F_Mean-CDH5_F_SEM, ymax = CDH5_F_Mean+CDH5_F_SEM), width = .1)+
  geom_jitter(height = 0, width = 0.2)+
  #stat_compare_means(label = "p.signif", method = "t.test", ref.group = "siRNA Control FM")+
  theme_minimal()+
  scale_x_discrete(name = "Group")+
  scale_fill_manual(values = c("#383131", "#383131", "#383131", "#ed3621", "#ed3621", "#faa79d", "#1c79eb", "#1c79eb", "#9bc2f2", "#46ebe8", "#46ebe8", "#ed3621", "#ed3621", "#faa79d", "#1c79eb", "#1c79eb", "#9bc2f2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank())+
  labs(title = "LT20/006 - HUVEC different conditions", 
       y= "2^ddCT normalized to d0 FM \n
       GAPDH and RPLP0")

DM_Tgfb2 <- filter(df, F_Group == "DM d7 TgfbII")
Rev_Tgfb2 <- filter(df, F_Group == "DM+FM d3+d7 TgfbII")
DM_Tgfb1 <- filter(df, F_Group == "DM d7 TgfbI + Il1b")
Rev_Tgfb1 <- filter(df, F_Group == "DM+FM d3+d7 TgfbI + Il1b")
DM_Tgfb2_L <- filter(df, F_Group == "DM d7 TgfbII + L-NAME")
Rev_Tgfb2_L <- filter(df, F_Group == "DM+FM d3+d7 TgfbII + L-NAME")
DM_Tgfb1_L <- filter(df, F_Group == "DM d7 TgfbI+Il1b+L-NAME")
Rev_Tgfb1_L <- filter(df, F_Group == "DM+FM d3+d7 TgfbI+Il1b+L-NAME")


t.test(x = DM_Tgfb2$SM22, y= Rev_Tgfb2$SM22)
t.test(x = DM_Tgfb1$SM22, y= Rev_Tgfb1$SM22)
t.test(x = DM_Tgfb2_L$SM22, y= Rev_Tgfb2_L$SM22)
t.test(x = DM_Tgfb1_L$SM22, y= Rev_Tgfb1_L$SM22)

t.test(x = DM_Tgfb2$CNN1, y= Rev_Tgfb2$CNN1)
t.test(x = DM_Tgfb1$CNN1, y= Rev_Tgfb1$CNN1)
t.test(x = DM_Tgfb2_L$CNN1, y= Rev_Tgfb2_L$CNN1)
t.test(x = DM_Tgfb1_L$CNN1, y= Rev_Tgfb1_L$CNN1)

t.test(x = DM_Tgfb2$CDH5, y= Rev_Tgfb2$CDH5)
t.test(x = DM_Tgfb1$CDH5, y= Rev_Tgfb1$CDH5)
t.test(x = DM_Tgfb2_L$CDH5, y= Rev_Tgfb2_L$CDH5)
t.test(x = DM_Tgfb1_L$CDH5, y= Rev_Tgfb1_L$CDH5)

##################### Figure Supplement 15e and f ######################################

Analysis <- read_csv("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/qPCR data Reversibility/Analysis_LT20_002.csv")

Desc <- read_csv("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/qPCR data Reversibility/Sample_desc_LT20_002.csv")

Analysis$Sample <- factor(Analysis$Sample)

Analysis <- unique(Analysis)

Desc$Sample <- factor(Desc$Sample)
Desc$Treatment <- factor(Desc$Treatment, levels = c("FM", "DM", "DM + FM"))
Desc$N <- factor(Desc$N)
Desc$Timepoint <- factor(Desc$Timepoint, levels = c("d0", "d7", "d7+d11", "d11", "d10", "d10+d14", "d14"))

df <- full_join(Analysis, Desc)

df <- df %>% group_by(Treatment, Timepoint) %>% mutate(SM22_F_Mean = mean(SM22), SM22_F_SEM = SEM(SM22), CNN1_F_Mean = mean(CNN1), CNN1_F_SEM = SEM(CNN1),
                                                       CDH5_F_Mean = mean(CDH5), CDH5_F_SEM = SEM(CDH5), F_Group = paste(Treatment, Timepoint))

df$F_Group <- factor(df$F_Group, levels = c("FM d0", "FM d11", "DM d7", "DM d11", "DM + FM d7+d11", "FM d14", "DM d10", "DM d14", "DM + FM d10+d14"))

ggplot(df, aes(x = F_Group, y= SM22, group = F_Group))+
  geom_col(aes(y = SM22_F_Mean, fill = F_Group), position = "dodge", color = "black", alpha = 0.4, width = .8)+
  geom_errorbar(aes(ymin = SM22_F_Mean-SM22_F_SEM, ymax = SM22_F_Mean+SM22_F_SEM), width = .1)+
  geom_jitter(height = 0, width = 0.2)+
  #stat_compare_means(label = "p.signif", method = "t.test", ref.group = "siRNA Control FM")+
  theme_minimal()+
  scale_x_discrete(name = "Group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank())+
  labs(title = "LT20/002 - HUVEC long term reversibility", 
       y= "2^ddCT normalized to d0 FM \n
       GAPDH and RPLP0")

ggplot(df, aes(x = F_Group, y= CNN1, group = F_Group))+
  geom_col(aes(y = CNN1_F_Mean, fill = F_Group), position = "dodge", color = "black", alpha = 0.4, width = .8)+
  geom_errorbar(aes(ymin = CNN1_F_Mean-CNN1_F_SEM, ymax = CNN1_F_Mean+CNN1_F_SEM), width = .1)+
  geom_jitter(height = 0, width = 0.2)+
  #stat_compare_means(label = "p.signif", method = "t.test", ref.group = "siRNA Control FM")+
  theme_minimal()+
  scale_x_discrete(name = "Group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank())+
  labs(title = "LT20/002 - HUVEC long term reversibility", 
       y= "2^ddCT normalized to d0 FM \n
       GAPDH and RPLP0")

ggplot(df, aes(x = F_Group, y= CDH5, group = F_Group))+
  geom_col(aes(y = CDH5_F_Mean, fill = F_Group), position = "dodge", color = "black", alpha = 0.4, width = .8)+
  geom_errorbar(aes(ymin = CDH5_F_Mean-CDH5_F_SEM, ymax = CDH5_F_Mean+CDH5_F_SEM), width = .1)+
  geom_jitter(height = 0, width = 0.2)+
  #stat_compare_means(label = "p.signif", method = "t.test", ref.group = "siRNA Control FM")+
  theme_minimal()+
  scale_x_discrete(name = "Group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "none", panel.grid.major.x = element_blank())+
  labs(title = "LT20/002 - HUVEC long term reversibility", 
       y= "2^ddCT normalized to d0 FM \n
       GAPDH and RPLP0")

d11_DM <- filter(df, F_Group == "DM d11")
d11_rev <- filter(df, F_Group == "DM + FM d7+d11")
d14_DM <- filter(df, F_Group == "DM d14")
d14_rev <- filter(df, F_Group == "DM + FM d10+d14")

t.test(x = d11_DM$SM22, y= d11_rev$SM22)
t.test(x = d14_DM$SM22, y= d14_rev$SM22)

t.test(x = d11_DM$CNN1, y= d11_rev$CNN1)
t.test(x = d14_DM$CNN1, y= d14_rev$CNN1)

t.test(x = d11_DM$CDH5, y= d11_rev$CDH5)
t.test(x = d14_DM$CDH5, y= d14_rev$CDH5)

