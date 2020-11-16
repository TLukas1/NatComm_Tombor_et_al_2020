# Initial analysis of Tracing dataset 

# This R-script provides all code from fastq to R object (MILENA.filtered)

# Please find all downstream analysis in markdown (Tombor_et_al_dataanlysis.Rdm)

# Tracing dataset AMI timecourse - (18-02-2019)

# Seurat object: TracingAMI

###### Required packages ########
# R-Version 3.4.4
# Seurat v 2.3.4
library(Seurat)
# Monocle v 2.6.4
library(monocle)
# reshape2 v 1.4.3
library(reshape2)
# ggplot2 v 3.1.0
library(ggplot2)
# scales v 0.5.0
library(scales)
# ggrepel v 0.8.0
library(ggrepel)
# tidyr v 0.8.2
library(tidyr)
# ggpubr v 0.1.8
library(ggpubr)
# biomart 2.34.2 
library(biomaRt)
library(tidyverse)

###### Load data ###### 
t <- c("/media/ATLAS_NGS_storage/Lukas/Tracer_Dec2018/103548_fastq/103548-001-001_cellrangerCount_CR3.0/001CR3/outs/filtered_feature_bc_matrix",
       "/media/ATLAS_NGS_storage/Lukas/Tracer_Dec2018/103548_fastq/103548-001-003_cellrangerCount_CR3.0/003CR3/outs/filtered_feature_bc_matrix",
       "/media/ATLAS_NGS_storage/Lukas/Tracer_Dec2018/103548_fastq/103548-001-004_cellrangerCount_CR3.0/004CR3/outs/filtered_feature_bc_matrix",
       "/media/ATLAS_NGS_storage/Lukas/Tracer_Dec2018/103548_fastq/103548-001-005_cellrangerCount_CR3.0/005CR3/outs/filtered_feature_bc_matrix", 
       "/media/ATLAS_NGS_storage/Lukas/ShallowSCSeq_Tracer_GFP/Fastq_Combined/Tracer_D1_CR3_combined/outs/filtered_feature_bc_matrix",
       "/media/ATLAS_NGS_storage/Lukas/ShallowSCSeq_Tracer_GFP/Fastq_Combined/Tracer_Hom_CR3_combined/outs/filtered_feature_bc_matrix",
       "/media/ATLAS_NGS_storage/Lukas/ShallowSCSeq_Tracer_GFP/Fastq_Combined/Tracer_D3_CR3_combined/outs/filtered_feature_bc_matrix",
       "/media/ATLAS_NGS_storage/Lukas/ShallowSCSeq_Tracer_GFP/Fastq_Combined/Tracer_D7_CR3_combined/outs/filtered_feature_bc_matrix")

anno <- c("d28", "d48", "Tam-Control", "d14", "d1", "Hom", "d3", "d7")
annos <- c("Hom", "d1", "d3", "d7", "d14", "d28", "d48", "Tam-Control")
col <- c("#686de0", "#dc143c","#bc8f8f", "#7ed6df", "#ff7979", "#f6e58d", "#badc58", "#6ab04c")
cols <- c("#f6e58d", "#ff7979", "#badc58", "#6ab04c", "#7ed6df", "#686de0", "#dc143c","#bc8f8f")
#load("/media/Storage/CompleteTracer_R/With FastqCombined/Tracing.Rds")
hash <- hash::hash(keys = t, values = anno)

for (i in t){
  x <- hash::values(hash[i])
  Matrices.10X <- Read10X(i)
  Seurats<- CreateSeuratObject(Matrices.10X, min.cells= 3, min.genes = 200)
  assign(x,Seurats)
}
rm(hash, Seurats, Matrices.10X, i, x, t)


###### Quality control and filtering #####

mito.genes <- NULL
for (i in anno){
  X <- get(i)
  mito.sample <- grep(pattern = "^mt-", x = rownames(x = X@data), value = TRUE)
  mito.genes <- c(mito.genes, mito.sample)
  rm(mito.sample, X)
}
mito.genes <- unique(mito.genes)
percent.mito = NULL
UMI = NULL
nGene = NULL
mito.p = NULL
label = NULL
for (i in anno){
  X <- get(i)
  percent.mito <- Matrix::colSums(X@raw.data[mito.genes, ])/Matrix::colSums(X@raw.data)
  X <- AddMetaData(object = X, metadata = percent.mito, col.name = "percent.mito")
  file.name <- paste("Quality_plots_",i,".svg", sep ="")
  if(!file.exists(file.name)){
    svg(filename = file.name, width = 10, height = 5.5)
    par(mfrow = c(1, 2))
    GenePlot(object = X, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    GenePlot(object = X, gene1 = "nUMI", gene2 = "nGene", cex.use = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    dev.off()
  }
  UMI <- c(UMI, X@meta.data$nUMI)
  nGene <- c(nGene, X@meta.data$nGene)
  mito.p <- c(mito.p, X@meta.data$percent.mito)
  label <- c(label, rep(i, times = length(rownames(X@meta.data))))
  assign(i, X)
  rm(X, file.name, percent.mito)
}

Quality <- data.frame(UMI, nGene, mito.p, label = factor(label, levels = annos))
rm(UMI,nGene,mito.p, label)

Quantile.low.UMI <- Quality %>% group_by(label) %>%
  summarise(UMI = list(enframe(quantile(UMI,probs = 0.01)))) %>%
  unnest

Quantile.high.UMI <- Quality %>% group_by(label) %>%
  summarise(UMI = list(enframe(quantile(UMI,probs = 0.95)))) %>%
  unnest

p1 <- ggplot(Quality, aes(x=UMI, fill = label))+
  geom_density()+
  geom_vline(data = Quantile.high.UMI, aes(xintercept = value), color = "red", size = 1.5)+
  geom_text(data = Quantile.high.UMI, aes(x=45000, y= 2e-04, label = name), color = "red", size = 8)+
  geom_text(data = Quantile.high.UMI, aes(x=40000, y= 2e-04, label = value), color = "red", size = 8)+
  geom_vline(data = Quantile.low.UMI, aes(xintercept = value), color = "blue", size = 1.5)+
  geom_text(data = Quantile.low.UMI, aes(x=35000, y= 2e-04, label = name), color = "blue", size = 8)+
  geom_text(data = Quantile.low.UMI, aes(x=30000, y= 2e-04, label = value), color = "blue", size = 8)+
  scale_fill_manual(values = cols)+
  coord_flip()+
  facet_wrap(~label, nrow = 2)+
  theme(legend.position = "none", axis.text.y = element_text(size = 24), axis.text = element_text(size = 8),axis.title = element_text(size=24, face = "bold"), axis.line = element_line(size = 1), strip.background.x = element_rect(fill = "white"),
        strip.text = element_text(size = 24))



Quantile.low.Gene <- Quality %>% group_by(label) %>%
  summarise(nGene = list(enframe(quantile(nGene,probs = 0.01)))) %>%
  unnest

Quantile.high.Gene <- Quality %>% group_by(label) %>%
  summarise(nGene = list(enframe(quantile(nGene,probs = 0.95)))) %>%
  unnest

p2 <- ggplot(Quality, aes(x=nGene, fill = label))+
  geom_density()+
  geom_vline(data = Quantile.high.Gene, aes(xintercept = value), color = "red", size = 1.5)+
  geom_text(data = Quantile.high.Gene, aes(x=5500, y= 4e-04, label = name), color = "red", size = 8)+
  geom_text(data = Quantile.high.Gene, aes(x=5000, y= 4e-04, label = value), color = "red", size = 8)+
  geom_vline(data = Quantile.low.Gene, aes(xintercept = value), color = "blue", size = 1.5)+
  geom_text(data = Quantile.low.Gene, aes(x=4500, y= 4e-04, label = name), color = "blue", size = 8)+
  geom_text(data = Quantile.low.Gene, aes(x=4000, y= 4e-04, label = value), color = "blue", size = 8)+
  scale_fill_manual(values = cols)+
  coord_flip()+
  facet_wrap(~label, nrow = 2)+
  theme(legend.position = "none", axis.text.y = element_text(size = 24), axis.text = element_text(size = 8),axis.title = element_text(size=24, face = "bold"), axis.line = element_line(size = 1), strip.background.x = element_rect(fill = "white"),
        strip.text = element_text(size = 24))

file.name = c("Filtering_Overview.svg")
if(!file.exists(file.name)){
  svg(filename = file.name, width = 10, height = 20)
  print(ggarrange(p1,p2, nrow = 2))
  dev.off()
}


# Filtering for nGene, nUMI and mitochondrial content

for (i in anno){
  A <- Quality %>% filter(label == i) %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = c(0.01, 0.95)))), 
              UMI = list(enframe(quantile(UMI,probs = c(0.01, 0.95))))) %>%
    unnest
  
  low.Gene <- A$value[1]
  low.UMI <- A$value1[1]
  high.Gene <- A$value[2]
  high.UMI <- A$value1[2]
  
  X <- get(i)
  
  X<- FilterCells(X, subset.names = c("nUMI", "nGene", "percent.mito"), high.thresholds = c(high.UMI, high.Gene, 0.1),
                  low.thresholds = c(low.UMI, low.Gene, 0))
  
  assign(i, X)
  
  rm(A,X,low.Gene,low.UMI,high.Gene,high.UMI)
}

# Filter genes which do not have more than 1 transcript per in at least 2 cells
Filtering.info <- list()
number.of.genes.raw <- NULL
for (i in anno){
  print(i)
  X <- get(i)
  more.than.2= NULL
  variance = NULL
  range = NULL
  mean = NULL
  n = length(rownames(X@data))
  number.of.genes.raw[i] <- n
  for(j in rownames(X@data)){
    x = X@raw.data[j, ]
    variance = c(variance, var(x))
    range = c(range, max(x)-min(x))
    mean = c(mean, mean(x))
    n2 = length(x[which(x > 1)])
    more.than.2 = c(more.than.2, n2)
    
    if (length(variance)%%100 == 0)
    {cat(length(variance), "of", n, "genes \n")}
  }
  Summary <- data.frame(Variance = variance, Range = range, Mean = mean, More.than.1 = more.than.2)
  rownames(Summary) <- rownames(X@data)
  
  Nondetectable <- subset(Summary, More.than.1 < 2)
  Genes.to.remove <- rownames(Nondetectable)
  Genes.to.keep <- rownames(subset(Summary,More.than.1 >= 2))
  Filtering.info[[i]] <- Summary
  
  X <- X@raw.data[Genes.to.keep,rownames(X@meta.data)]
  X<- CreateSeuratObject(X,project = "Tracer")
  assign(i, X)
  
  rm(Summary, X)
}



# Re-assign again percent mito with new matix
mito.genes = NULL
for (i in anno){
  X <- get(i)
  mito.sample <- grep(pattern = "^mt-", x = rownames(x = X@data), value = TRUE)
  percent.mito <- Matrix::colSums(X@raw.data[mito.sample, ])/Matrix::colSums(X@raw.data)
  X <- AddMetaData(object = X, metadata = percent.mito, col.name = "percent.mito")
  assign(i,X)
}
# Normalization, FindVariableGenes, ScaleData

for (i in anno){
  X <- get(i)
  X<- NormalizeData(object = X, normalization.method = "LogNormalize", 
                    scale.factor = 10000)
  X <- FindVariableGenes(object = X, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  assign(paste(i, "_copy_reg", sep = ""), X)
  X <- ScaleData(object = X, vars.to.regress = c("nUMI", "percent.mito"))
  assign(i, X)
}


# How many cells have been filtered
cell.counts <- NULL
for (i in anno){
  X <- get(i)
  
  n = length(rownames(X@meta.data))
  cell.counts <- c(cell.counts, n)
  rm(n,X)
}

names(cell.counts) <- anno

ecdf_fun <- function(x,perc) ecdf(x)(perc)
percentiles <- NULL
ranges <- NULL
means <- NULL
cells <- NULL
for (i in anno){
  table <- Filtering.info[[i]]
  mean.GFP<- table["EGFP", "Mean"]
  percentile.GFP<- ecdf_fun(table$Mean, mean.GFP)
  percentiles[i] <- percentile.GFP
  ranges[i] <- table["EGFP", "Range"]
  means[i] <- mean.GFP
  cells[i] <- table["EGFP", "More.than.1"]
  
  rm(table, mean.GFP, percentile.GFP)
}

rm(file.name, i, mito.genes, Quality, Quantile.high.Gene, Quantile.high.UMI, Quantile.low.Gene, Quantile.low.UMI, p1, p2)




##### CCA clustering #####

Objects <- list(Hom, d1, d3, d7, d14, d28, d48, `Tam-Control`)

var.genes <- NULL
for (i in anno){
  X <- get(i)
  s = head(X@var.genes, 1000)
  var.genes <- c(var.genes, s)
}

var.genes <- unique(var.genes)
var.genes <- intersect(var.genes, rownames(Hom@scale.data))
var.genes <- intersect(var.genes, rownames(d1@scale.data))
var.genes <- intersect(var.genes, rownames(d3@scale.data))
var.genes <- intersect(var.genes, rownames(d7@scale.data))
var.genes <- intersect(var.genes, rownames(d14@scale.data))
var.genes <- intersect(var.genes, rownames(d28@scale.data))
var.genes <- intersect(var.genes, rownames(d48@scale.data))
var.genes <- intersect(var.genes, rownames(`Tam-Control`@scale.data))

Tracing<- RunMultiCCA(Objects, add.cell.ids = annos, genes.use = var.genes, num.ccs = 30)
MetageneBicorPlot(Tracing, grouping.var = "Sample", dims.eval = 1:30)

DimHeatmap(object = Tracing, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)

Tracing<- AlignSubspace(Tracing, grouping.var = "orig.ident", dims.align = 1:25)

Tracing <- FindClusters(object =Tracing, reduction.type = "cca.aligned", dims.use = 1:20, 
                        resolution = 0.6, print.output = 1, save.SNN = F)
# TSNE Clustering
Tracing <- RunTSNE(Tracing, do.fast = T, reduction.use = "cca")
Tracing <- SetIdent(Tracing, ident.use = Tracing@meta.data$GFP)
TSNEPlot(Tracing, do.label = T)
FeaturePlot(Tracing, features.plot = c("EGFP", "tdTomato", "Cdh5", "Pecam1"), cols.use = c("lightskyblue", "red4"), pt.size = 0.8)

Tracing<-SetI(Tracing, "orig.ident")
TSNEPlot(Tracing, do.label = F)

###### Estimating GFP threshold by density function ##### 

A <- data.frame(GFP.normalized = Tracing@data["EGFP", ], Timepoint = factor(Tracing@meta.data$orig.ident, levels = annos), Cluster = factor(Tracing@meta.data$res.0.8, levels = clus))


GFP.threshold <- NULL

for(i in anno){
  
  B <- subset(A, Timepoint == i)
  
  test <- density(B$GFP.normalized, n= length(B$GFP.normalized))
  
  y.minimum <- min(test$y[test$x > 2 & test$x < 4])
  x.minimum <- which(test$y == y.minimum)
  
  file.name = paste(i,"_GFP.threshold.svg")
  if(!file.exists(file.name)){
    svg(filename = file.name, width = 10, height = 4)
    
    print(ggplot(B, aes(x= GFP.normalized)) + 
            geom_density() + 
            geom_vline(xintercept = test$x[x.minimum], color = "blue3", size = 1)+
            geom_label(x = test$x[x.minimum]*0.9, y= y.minimum + 0.2, label = round(test$x[x.minimum], digits = 4), color = "blue")+
            ggtitle(i))
    dev.off()
  }
  GFP.threshold[i] <- test$x[x.minimum]
  
}

distinguisher <- NULL
for (i in anno){
  X <- SubsetData(Tracing, cells.use = WhichCells(Tracing, cells.use = which(Tracing@meta.data$orig.ident == i)))
  for (j in colnames(X@data)){
    
    if(Tracing@data["EGFP",j] >= GFP.threshold[i]){
      distinguisher[j] <- "GFP.positive"
    }
    else{distinguisher[j] <- "GFP.negative"}
  }
  rm(X)
}

Tracing <- AddMetaData(Tracing, distinguisher, "GFP")

A <- data.frame(GFP.normalized = Tracing@data["EGFP", ], Timepoint = factor(Tracing@meta.data$orig.ident, levels = annos), 
                Cluster = factor(Tracing@meta.data$res.0.8, levels = clus), GFP = Tracing@meta.data$GFP)

ggplot(A, aes(x=Timepoint, y= GFP.normalized, color = GFP))+ 
  geom_jitter(height = 0, width = 0.3, alpha = 0.6)+ 
  scale_color_manual(values= c("grey60", "darkgreen"))


Tracing <- SetI(Tracing, "GFP")
TSNEPlot(Tracing, colors.use = c("grey60", "darkgreen"))