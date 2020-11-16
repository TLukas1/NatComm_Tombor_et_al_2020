setwd("/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/DNA Methylation/Regulatory Elements mit Nina Baumgarten")

cmp1 <- read.table("cmp_1_diffMethInfo_EpiRegio.txt", header = T)
cmp3 <- read.table("cmp_3_diffMethInfo_EpiRegio.txt", header = T)

cmp1$regulation <- ifelse(cmp1$mean.quot.log2 < 0, "methylated", "demethylated")
cmp3$regulation <- ifelse(cmp3$mean.quot.log2 < 0, "methylated", "demethylated")

join <- inner_join(cmp1, cmp3, by= c("Gene_ID", "REM_ID", "region", "Gene_symbol"),suffix = c("_cmp1", "_cmp3"))

full.join <- full_join(cmp1, cmp3, by= c("Gene_ID", "REM_ID", "region", "Gene_symbol"),suffix = c("_cmp1", "_cmp3"))

write.csv(full.join, file = "HUVEC_DNA_MET_ENHANCER.csv")

join %>% group_by(regulation_cmp1, regulation_cmp3) %>% summarize(n = n()) %>%
  mutate(Regulation = paste(regulation_cmp1, regulation_cmp3, sep = "/")) %>%
  ggplot(aes(x = Regulation, y = n)) +
  geom_col(color = "black", fill = c("blue", "red"))+
  theme_classic()+
  theme(text = element_text(size = 15, color = "black"), axis.text = element_text(size = 12, color = "black"))+
  labs(y = "Number of common enhancer elements", x = "FM - DM / DM - FM", title = "DNA methylation in reversibility")+
  geom_text(aes(y = n+15, label = paste("n=",n)), size = 6)

library(ComplexHeatmap)
library(gplots)
library(RColorBrewer)
library(circlize)


input_file <- "InputHeatmap_cmp1_cmp3_common.txt"
data <- read.delim(input_file, header = TRUE, sep = "\t", row.names = 1) #REMs x TFs and last column is activity per REM


#print(data)
rnames <- rownames(data)
data_line = data[, ncol(data)]
#print(data_line)

data_heatmap = as.matrix(data[,- ncol(data)])
print(max(data_heatmap))
#print(data_heatmap)

helper = rev(brewer.pal(11, "RdBu"))
#TODO: hie rnoch anpassen
col_fun = colorRamp2(c(0.70, 0.5, 0.25, 0.15,  0.1,  0, -0.1, -0.15, -0.25,-0.5, -0.70), c(helper[1], helper[2],helper[3], helper[4], helper[5], helper[6], helper[7], helper[8], helper[9], helper[10], helper[11]))
Heatmap(data_heatmap, name = "?",
        heatmap_legend_param = list(title = "values", at = c(0.7, 0.5, 0.25, 0.1, 0,-0.1, -0.25, -0.5, -0.7)),
        col = col_fun,
        cluster_rows = as.dendrogram(hclust(dist(data_heatmap, method = "manhattan"),method="ward.D")),
        show_row_names = F,
        row_title = "differentially methylated REMs", 
        row_title_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 1),
        row_title_side = "left",
        show_row_dend = TRUE,
        row_dend_width = unit(3, "cm"),
        
        column_title = "Conditions", 
        column_title_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 10),
        show_column_dend = FALSE,
        column_title_side = "bottom"
)

library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(forcats)
library(stringr)

hs_df = msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

genes <- join[which(join$regulation_cmp1 == "demethylated"), "Gene_symbol"]
genes <- genes[which(genes %in% hs_df$gene_symbol)]
enrich_results <- enricher(genes, TERM2GENE = hs_df, pAdjustMethod = "fdr")

TF_PASTA <- read.table("PASTAA_fdr.txt", header = T)

TF_PASTA_f <- TF_PASTA %>% filter(fdr < 0.01) %>% arrange(-p_values)

TF_PASTA_f$TF <- ifelse(TF_PASTA_f$TFs %in% c("KLF4", "GATA2", "ERG", "HOXA9"), "TF", "")

ggplot(TF_PASTA_f, aes(x = fct_reorder(TFs, p_values), y = -log(p_values), fill = TF))+
  geom_col(color = "black")+
  theme_classic()+
  scale_fill_manual(values = c("grey80", "orange"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x = "TFs", y = "-log(p-value)")

REM_TF_anno <- read.table("methylatedREMs_TFsAnnotation.txt", header = T)

short.join <- data.frame(sequence_name = join$region, Gene_symbol = join$Gene_symbol)

REM_TF_anno <- full_join(REM_TF_anno, short.join)

full.join %>% replace_na(replace = list(regulation_cmp1 = "-", regulation_cmp3 = "-"))

full.x <- full.join %>% group_by(regulation_cmp1, regulation_cmp3) %>% summarize(n = n()) %>% mutate(freq = n/sum(n)) %>% replace_na(replace = list(regulation_cmp1 = "-", regulation_cmp3 = "-"))
full.x %>% mutate(Regulation = paste(regulation_cmp1, regulation_cmp3, sep = "/")) %>% filter(regulation_cmp1 == "demethylated") %>%
  ggplot(aes(x = "Regulation", y = freq)) +
  geom_col(color = "black", fill = c("blue", "red"))+
  coord_polar("y", start = 0)+
  theme_classic()+
  theme(text = element_text(size = 15, color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank())

ge <- c("G1/G2/G3")

spit <- unlist(str_split(ge, pattern = "/"))

top_n(n=10)

Simone.DEG <- read.csv("/media/ATLAS_NGS_storage/Simone/SRA-Submission-SimoneGlaser-12.2019/RNA-SEQ/Gene_expression.tsv", sep = "\t")

KLF4.regulated <- read.delim("KLF4_associatedGenes_ids.txt", header = F)
ERG.regulated <- read.delim("ERG_associatedGenes.txt", header = F, )
GATA2.regulated <- read.delim("GATA2_associatedGenes.txt", header = F)

Simone.Filter <- Simone.DEG %>% filter(gene_id %in% c(as.character(KLF4.regulated$V1), as.character(ERG.regulated$V2), as.character(GATA2.regulated$V2)) & significant_DMT.VM == "yes")

Simone.melt <- reshape2::melt(Simone.Filter[,1:13], variable.name = "Sample")

Simone.melt$Condition <- ifelse(str_detect(Simone.melt$Sample, pattern = "VM"), "FM", "DM")
Simone.melt$Condition <- factor(Simone.melt$Condition, levels = c("FM", "DM"))

ggplot(Simone.melt, aes(x = Condition, y = value, fill = Condition))+
  geom_boxplot()+
  geom_jitter(height = 0)+
  facet_wrap(~gene_short_name, ncol = 2, scales = "free")+
  theme_classic()+
  scale_fill_manual(values = c("red", "orange"))+
  theme_classic()+
  labs(y = "FPKM")

join.filter <- dplyr::filter(join, Gene_symbol %in% as.character(Simone.Filter$gene_short_name))

Met.plot <- data.frame(Region = join.filter$region, REM = join.filter$REM_ID, Gene = join.filter$Gene_symbol, Gene_ID = join.filter$Gene_ID, Control_to_FM = join.filter$regulation_cmp1, Fold_change = round(join.filter$mean.quot.log2_cmp1,2), DM_to_Reversible =join.filter$regulation_cmp3, Fold.change = round(join.filter$mean.quot.log2_cmp3,2), TF_binding_motif = c("GATA2", "ERG", "KLF4", "KLF4"))

colnames(Met.plot) <- c("chromosome region", "REM ID", "associated gene", "ENSEMBL ID", "Control to DM", "Change (mean difference log2)", "DM to DM+FM", "Change (mean difference log2)", "Predicted TF binding")

library(kableExtra)

kableExtra::kable(Met.plot) %>% column_spec(3, italic = T) %>% column_spec(5, background = "orange", color = "white", bold = T) %>% 
  column_spec(7, background = "red", color = "white", bold = T)  %>% kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, color = "black") %>% row_spec(0:4, align = "center") %>% column_spec(9, bold = T) %>% 
  column_spec(1:4, color = "black") %>% column_spec(6, color = "black") %>% column_spec(8:9, color = "black") %>% as_image(file = "R7C table.png", )

