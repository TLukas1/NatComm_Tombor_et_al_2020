# Innitial analysis Rosenthal data
# This R-script provides all code from fastq to R object (MILENA.filtered)

# Please find all downstream analysis in markdown (Tombor_et_al_dataanlysis.Rdm)

# Load libraries

library("stringr")
library("Seurat")
library("dplyr")
library(ggplot2)
library(ggpubr)
library("biomaRt")
library(lubridate)
library(GO.db)
library(monocle)

##### Custom functions ####### 

#' Import and combine several Single cell sequencing experiments into Seurat
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineSeuratObjects <- function(pathways,ids){
  Importer <- function(pathway,id) {
    Ten_X <- Read10X(pathway)
    seuratObject =CreateSeuratObject(raw.data = Ten_X, min.cells= 3, min.genes = 200)
    seuratObject<-RenameCells(object = seuratObject, add.cell.id = id)
    cat("Imported ", length(seuratObject@ident), " cells from ", pathway, "with ID ", id, "\n")
    return(seuratObject)
  }
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  for (i in 1:length(pathways)) {
    if (i<2) { 
      seuratObject1<-Importer(pathways[i],ids[i])
      next
    }
    seuratObject2<-Importer(pathways[i],ids[i])
    seuratObject1 <- MergeSeurat(object1 = seuratObject1, object2 = seuratObject2) 
  }
  cat("Merged Seurat object contains ", length(seuratObject1@ident)," cells\n")
  return(seuratObject1)
}

#' Use barcode annotation in Seurat object to create a metadata entry
#' @author Lukas Tombor
#' @param Seurat A Seurat object
#' @param ids Vector of strings that have been assigned in combineSeuratObjects
#' @param name Name of the meta data entry
#' @return Seurat object with meta data entry

LibrarytoMetaData <- function(Seurat,ids, name){
  barcodes <- rownames(Seurat@meta.data)
  Lab = NULL
  for (i in ids){
    label <- grep(pattern = paste("^",i, sep=""), x = barcodes, value = T)
    if (length(label) == 0) {stop("Could not use Seurat object with given ids")}
    Lab <- c(Lab, rep(i, times = length(label)))
  }
  Lab <- factor(Lab, levels=ids)
  names(Lab) <- barcodes
  Seurat <- AddMetaData(Seurat, metadata = Lab, col.name = name)
  return(Seurat)
}

#' Sets parameter in meta.data to ident
#' @author Lukas Tombor
#' @param Seurat A Seurat object
#' @param meta Entry in the meta.data slot that should be set to ident
#' @return Seurat object with new ident

SetI <- function(Seurat, meta){
  ident <- Seurat@meta.data[, meta]
  Seurat <- SetIdent(Seurat, ident.use = ident)
  return(Seurat)
}
#' Converts a human gene list into mouse
#' @author Lukas Tombor
#' @param x A vector containing a list of human genes
#' @return Converted vector
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
##### Import Data ####### 
cellranger_pipestance_path_homeostasis <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17010/MF17010/outs/filtered_gene_bc_matrices/mm10"
cellranger_pipestance_path_d1 <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17013/MF17013/outs/filtered_gene_bc_matrices/mm10"
cellranger_pipestance_path_d3 <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17014/MF17014/outs/filtered_gene_bc_matrices/mm10"
cellranger_pipestance_path_d5 <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17015/MF17015/outs/filtered_gene_bc_matrices/mm10"
cellranger_pipestance_path_d7 <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17016/MF17016/outs/filtered_gene_bc_matrices/mm10"
cellranger_pipestance_path_d14 <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17017/MF17017/outs/filtered_gene_bc_matrices/mm10"
cellranger_pipestance_path_d28 <- "/media/ATLAS_NGS_storage/Lukas/Publication_Data_Dec_2019/MILENA_BL6_DATA/cellranger/MF17018/MF17018/outs/filtered_gene_bc_matrices/mm10"

pipestances = c(cellranger_pipestance_path_homeostasis, cellranger_pipestance_path_d1, cellranger_pipestance_path_d3, cellranger_pipestance_path_d5, 
                cellranger_pipestance_path_d7, cellranger_pipestance_path_d14, cellranger_pipestance_path_d28)

anno <- c("Homeostasis", "d1_post", "d3_post", "d5_post", "d7_post", "d14_post", "d28_post")

MILENA <- combineSeuratObjects(pipestances, anno)
MILENA <- LibrarytoMetaData(MILENA, anno, "Timepoint")

before<- table(MILENA@meta.data$Timepoint)
# Quality control 

mito.genes <- grep(pattern = "^mt-", x = rownames(x = MILENA@data), value = TRUE)
percent.mito <- Matrix::colSums(MILENA@raw.data[mito.genes, ])/Matrix::colSums(MILENA@raw.data)
MILENA <- AddMetaData(object = MILENA, metadata = percent.mito, col.name = "percent.mito")

MILENA <- SetI(MILENA, "orig.ident")
VlnPlot(MILENA, features.plot = c("nGene", "nUMI", "percent.mito"), nCol=3)
par(mfrow = c(1, 2))
GenePlot(object = MILENA, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.5)
abline(h= 0.1, col= "red")
GenePlot(object = MILENA, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.5)
rect(xleft = 0, xright = 43300, ybottom = 0, ytop = 400, density = 50, col = "white")
rect(xleft = 18000, xright =43300, ybottom = 0, ytop = 6250, density = 50, col = "white")
rect(xleft = 0, xright = 18000, ybottom = 4500, ytop = 6250, density = 50, col = "white")
rect(xleft = 0, xright = 18000, ybottom = 400, ytop = 4500, density =0, col = "red")

hist(MILENA@meta.data$nGene, breaks = 100, main = "nGene", xlab= "nGene")
abline(v = 400, col = "red")
abline(v = 4500, col = "red")
text(x= 5000, y = 1000, labels= "Filtered: \n nGene > 4500 \n nGene < 400")

hist(MILENA@meta.data$nUMI, breaks = 100, main = "nUMI", xlab="nUMI")
abline(v = 18000, col = "red")
text(x= 35000, y = 2500, labels= "Filtered: \n nUMI > 18000 \n")

# Apply filter

MILENA <- FilterCells(MILENA, subset.names = c("nUMI", "nGene", "percent.mito"), low.thresholds = c(0, 400, 0), high.thresholds = c(18000, 4500, 0.1))
after <- table(MILENA@meta.data$Timepoint)
diff <- before-after
p <- round(diff/before, digits = 4)
compare <- data.frame(cbind(before, after, diff, p))


# Filtering for genes 
more.than.2= NULL
variance = NULL
range = NULL
mean = NULL
more.than.2 = NULL
n = length(rownames(MILENA@data))

for(i in rownames(MILENA@data)){
  x = MILENA@raw.data[i, ]
  variance = c(variance, var(x))
  range = c(range, max(x)-min(x))
  mean = c(mean, mean(x))
  n2 = length(x[which(x > 1)])
  more.than.2 = c(more.than.2, n2)
  
  if (length(variance)%%100 == 0)
  {cat(length(variance), "of", n, "genes \n")}
}

Summary <- data.frame(Variance = variance, Range = range, Mean = mean, More.than.1 = more.than.2)
rownames(Summary) <- rownames(MILENA@data)

Nondetectable <- subset(Summary, More.than.1 < 2)
Genes.to.remove <- rownames(Nondetectable)
Genes.to.keep <- rownames(subset(Summary,More.than.1 >= 2))

X <- MILENA@raw.data[Genes.to.keep,rownames(MILENA@meta.data)]

MILENA.filtered <- CreateSeuratObject(X, meta.data = MILENA@meta.data)


# Normalization 

MILENA.filtered <- NormalizeData(object = MILENA.filtered, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

# Find variable genes
MILENA.filtered <- FindVariableGenes(object = MILENA.filtered, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# Scale Data
MILENA.filtered <- ScaleData(object = MILENA.filtered, vars.to.regress = c("nUMI", "percent.mito"))

# Run PCA
MILENA.filtered <- RunPCA(object=MILENA.filtered, pc.genes = MILENA.filtered@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)
PCElbowPlot(MILENA.filtered)
# 10 Dimensions to use for Find Cluster
MILENA.filtered <- FindClusters(object =MILENA.filtered, reduction.type = "pca", dims.use = 1:10, 
                                resolution = 0.6, print.output = 1, save.SNN = F)
# TSNE Clustering
MILENA.filtered <- RunTSNE(MILENA.filtered, do.fast = T)