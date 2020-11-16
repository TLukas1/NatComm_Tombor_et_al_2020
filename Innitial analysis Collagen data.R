# Create Seurat Objects 

#setwd("/media/Storage/Collagen Tracing/Quality Control")

Paths <- c("/media/Helios_scStorage/Lukas/103992/103992-001-001/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-002/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-003/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-004/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-005/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-006/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-007/Solo.out/Gene/filtered",
           "/media/Helios_scStorage/Lukas/103992/103992-001-008/Solo.out/Gene/filtered")

Sample = c("S_1013", "S_1015", "S_1017", "S_1019", "S_1022", "S_1023", "S_1024", "S_1026")

cols <- c("#76BA1B", "#4C9A2A", "#ACDF87", "#522888", "#66369D", "#967BB6", 
          "#CDC7D9", "#FF4E41")

# Additonal Sample information in table
tbl = data.frame(Sample = factor(Sample, levels = Sample),
                 Sex = c("male", "male", "female", "female", "male", "male", "male", "male"),
                 Age = c("Old", "Old", "Old", "Old", "Young", "Young", "Young", "Young"),
                 AMI = c("AMI", "Hom", "AMI", "Hom", "AMI", "AMI", "Hom", "Hom"))

hash <- hash::hash(keys = Paths, values = Sample)

for (i in Paths){
  x <- hash::values(hash[i])
  print(i)
  list.files(i)
  Matrices.10X <- Read10X(i)
  Seurats<- CreateSeuratObject(Matrices.10X, min.cells= 3, min.features = 200)
  assign(x,Seurats)
}
rm(hash, Seurats, Matrices.10X, i, x)

mito.genes <- NULL
for (i in Sample){
  X <- get(i)
  mito.sample <- grep(pattern = "^mt-", x = rownames(x = X), value = TRUE)
  mito.genes <- c(mito.genes, mito.sample)
  rm(mito.sample, X)
}

mito.genes <- unique(mito.genes)
percent.mito = NULL
UMI = NULL
nGene = NULL
mito.p = NULL
label = NULL
for (i in Sample){
  print(i)
  X <- get(i)
  #X <- get("S_1019")
  
  raw.data <- GetAssayData(X, slot = "counts")
  percent.mito <- Matrix::colSums(raw.data[mito.genes, ])/Matrix::colSums(raw.data)
  X$percent.mito <- percent.mito
  file.name <- paste("Quality_plots_",i,".svg", sep ="")
  if(!file.exists(file.name)){
    svg(filename = file.name, width = 10, height = 5.5)
    #par(mfrow = c(1, 2))
    p1 <- FeatureScatter(object = X, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = "black")
    p2 <- FeatureScatter(object = X, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "black")
    
    print(ggarrange(p1,p2, legend = "none"))
    
    dev.off()
  }
  UMI <- c(UMI, X$nCount_RNA)
  nGene <- c(nGene, X$nFeature_RNA)
  mito.p <- c(mito.p, X$percent.mito)
  label <- c(label, rep(i, times = length(rownames(X@meta.data))))
  assign(i, X)
  rm(X, file.name, percent.mito)
}

Quality <- data.frame(UMI, nGene, mito.p, label = factor(label, levels = Sample))
rm(UMI,nGene,mito.p, label)

Quantile.low.UMI <- Quality %>% group_by(label) %>%
  summarise(UMI = list(enframe(quantile(UMI,probs = 0.01)))) %>%
  unnest

Quantile.high.UMI <- Quality %>% group_by(label) %>%
  summarise(UMI = list(enframe(quantile(UMI,probs = 0.95)))) %>%
  unnest

p1 <- ggplot(Quality, aes(x=UMI, fill = label))+
  geom_density()+
  geom_vline(data = Quantile.high.UMI, aes(xintercept = value), color = "red", size = .5)+
  geom_text(data = Quantile.high.UMI, aes(x=65000, y= 1e-04, label = name), color = "red", size = 4)+
  geom_text(data = Quantile.high.UMI, aes(x=55000, y= 1e-04, label = value), color = "red", size = 4)+
  geom_vline(data = Quantile.low.UMI, aes(xintercept = value), color = "blue", size = .5)+
  geom_text(data = Quantile.low.UMI, aes(x=40000, y= 1e-04, label = name), color = "blue", size = 4)+
  geom_text(data = Quantile.low.UMI, aes(x=30000, y= 1e-04, label = value), color = "blue", size = 4)+
  scale_fill_manual(values = cols)+
  coord_flip(xlim = c(0, 75000))+
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
  geom_vline(data = Quantile.high.Gene, aes(xintercept = value), color = "red", size = .5)+
  geom_text(data = Quantile.high.Gene, aes(x=10000, y= 4e-04, label = name), color = "red", size = 4)+
  geom_text(data = Quantile.high.Gene, aes(x=9000, y= 4e-04, label = value), color = "red", size = 4)+
  geom_vline(data = Quantile.low.Gene, aes(xintercept = value), color = "blue", size = .5)+
  geom_text(data = Quantile.low.Gene, aes(x=7500, y= 4e-04, label = name), color = "blue", size = 4)+
  geom_text(data = Quantile.low.Gene, aes(x=6500, y= 4e-04, label = value), color = "blue", size = 4)+
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


for (i in Sample){
  A <- Quality %>% filter(label == i) %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = c(0.01, 0.95)))), 
              UMI = list(enframe(quantile(UMI,probs = c(0.01, 0.95))))) %>%
    unnest
  
  low.Gene <- A$value[1]
  low.UMI <- A$value1[1]
  high.Gene <- A$value[2]
  high.UMI <- A$value1[2]
  test <- 0.1
  
  X <- get(i)
  
  X<- subset(x= X, subset = nCount_RNA < high.UMI & 
               nCount_RNA > low.UMI & 
               nFeature_RNA < high.Gene & 
               nFeature_RNA > low.Gene & 
               percent.mito < test)
  
  assign(i, X)
  
  rm(A,X,low.Gene,low.UMI,high.Gene,high.UMI, test )
}

Number_after_filter <- NULL

for(i in Sample){
  X <- get(i)
  Number_after_filter[i] <- length(X$orig.ident)
}

A <- data.frame(tbl, nCells = Number_after_filter)

df.plot <- ggplot(A, aes(x= Sample, y= nCells, color = Sample))+
  geom_point(size = 3)+
  stat_summary(geom = "text", fun.y = function(x) x+200, size = 5, color = "darkred", aes(label = nCells))+
  theme_bw()+
  scale_y_continuous(position = "right")+
  scale_color_manual(name = "Animal ID", values = cols)+
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks = element_blank())+
  labs(title = "Number of cells after filtering")


df.table1 <- ggplot(A, aes(x= Sample, y= 0, label = Sex))+
  geom_text(size = 4, color = "grey70")+
  theme_minimal() + 
  scale_y_continuous(breaks=NULL, name = "Sex") +
  theme(panel.grid.major = element_blank(), legend.position = "none",
        panel.border = element_blank(), axis.text.x =  element_blank(),
        axis.ticks =  element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0, hjust = 0))

df.table2 <- ggplot(A, aes(x= Sample, y= 0, label = Age))+
  geom_text(size = 4, color = "lightblue")+
  theme_minimal() + 
  scale_y_continuous(breaks=NULL, name = "Age(m)") +
  theme(panel.grid.major = element_blank(), legend.position = "none",
        panel.border = element_blank(), axis.text.x =  element_blank(),
        axis.ticks =  element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0, hjust = 0))

df.table3 <- ggplot(A, aes(x= Sample, y= 0, label = AMI))+
  geom_text(size = 4, color = "darkorange")+
  theme_minimal() + 
  scale_y_continuous(breaks=NULL, name = "AMI") +
  theme(panel.grid.major = element_blank(), legend.position = "none",
        panel.border = element_blank(), axis.text.x =  element_blank(),
        axis.ticks =  element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0, hjust = 0))

gA <- ggplotGrob(df.plot)
gB <- ggplotGrob(df.table1)
gC <- ggplotGrob(df.table2)
gD <- ggplotGrob(df.table3)

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3], gC$widths[2:3], gD$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
gC$widths[2:3] <- as.list(maxWidth)
gD$widths[2:3] <- as.list(maxWidth)



file.name = c("NumberOfCellAfterFiltering.svg")
if(!file.exists(file.name)){
  svg(filename = file.name, width = 10, height = 10)
  print(ggarrange(p1,p2, nrow = 2))
  print(grid.arrange(gA, gB, gC, gD,  ncol=1, heights=c(10, .8, .8, .8)))
  dev.off()
}

'%!in%' <- function(x,y)!('%in%'(x,y))
# filter genes on the number of cells expressing
# modifies the raw.data slot as well now
#' Filters out Genes which are not expressed 
#' @param object A Seurat object
#' @param min.cells Number of cells which at least must hit the UMI threshold 
#' @param UMI.tsh UMI threshold, keeps genes with > UMI.tsh
#' @return A filtered Seurat object with @expressed.genes and @removed.genes
#' @author Lukas Tombor (based on https://github.com/satijalab/seurat/issues/147)

RemoveGenes <- function(object, min.cells = 2, UMI.tsh = 1){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  Data_slot <- as.matrix(GetAssayData(object))
  genes.use <- rownames(Data_slot)
  num.cells <- rowSums(Data_slot > UMI.tsh)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  object@assays$RNA@counts <- object@assays$RNA@counts[genes.use, ]
  object@assays$RNA@data <- object@assays$RNA@data[genes.use, ]
  removed.genes <- names(num.cells[which(names(num.cells)%!in%genes.use)])
  attributes(object)$removed.genes<- removed.genes
  attributes(object)$expressed.genes <- genes.use
  return(object)
}

for(i in Sample){
  X <- get(i)
  X <- RemoveGenes(X)
  assign(i,X)
}

rm(X)

# Standard Seurat pipeline 

for(i in Sample){
  X <- get(i)
  X <- NormalizeData(X, verbose = F)
  X <- FindVariableFeatures(X, selection.method = "vst", verbose = F)
  assign(i, X)
}


### Define Number of GFP positive cells

GFP.threshold <- NULL

for(i in Sample){
  
  X <- get(i)
  X <- get("S_1013")
  
  
  GFP.df <- data.frame(GetAssayData(X)["EGFP", ])
  colnames(GFP.df) <- "GFP.normalized"
  
  test <- density(GFP.df$GFP.normalized, n= length(GFP.df$GFP.normalized))
  
  y.minimum <- min(test$y[test$x > 2 & test$x < 4])
  x.minimum <- which(test$y == y.minimum)
  
  file.name = paste(i,"_GFP.threshold.svg")
  if(!file.exists(file.name)){
    svg(filename = file.name, width = 10, height = 4)
    
    print(ggplot(GFP.df, aes(x= GFP.normalized)) + 
            geom_density() + 
            geom_vline(xintercept = test$x[x.minimum], color = "blue3", size = 1)+
            geom_label(x = test$x[x.minimum]*0.9, y= y.minimum + 0.2, label = round(test$x[x.minimum], digits = 4), color = "blue")+
            ggtitle(i))
    dev.off()
  }
  GFP.threshold[i] <- test$x[x.minimum]
  
}


for (i in Sample){
  distinguisher <- NULL
  X <- get(i)
  for (j in colnames(GetAssayData(X))){
    
    if(GetAssayData(X)["EGFP",j] >= GFP.threshold[i]){
      distinguisher[j] <- "GFP.positive"
    }
    else{distinguisher[j] <- "GFP.negative"}
  }
  
  X$GFP <- distinguisher
  assign(i,X)
  rm(distinguisher)
}

Summary = NULL
for(i in Sample){
  X <- get(i)
  A <- as.data.frame(table(X$GFP))
  A$Sample <- rep(i, 2)
  A$percent <- A$Freq/sum(A$Freq)
  Summary <- data.frame(rbind(Summary, A))
}

colnames(Summary)[1] <- "GFP"

file.name = c("NumberOfGFPPositive.svg")
if(!file.exists(file.name)){
  svg(filename = file.name, width = 10, height = 10)
  p<-ggplot(subset(Summary, subset = GFP == "GFP.positive"), aes(x=Sample, y= percent, color = Sample))+ 
    geom_point(size = 3)+
    scale_color_manual(name = "Animal ID", values = cols)+
    scale_y_continuous(labels = scales::percent_format())+
    labs(title = "Number of GFP+ cells", y = "Percent of GFP+ cells per timepoint")+
    facet_grid(~Sample, scales = "free_x")+  
    theme_bw() + 
    theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold", size = 20))
  
  print(p)
  dev.off()
}


for(i in Sample){
  X <- get(i)
  X@meta.data$orig.ident <- rep(i, nrow(X@meta.data))
  assign(i,X)
}
Col.tracing <- list(S_1013, S_1015, S_1017, S_1019, S_1022, S_1023, S_1024, S_1026)

Col.anchors <- FindIntegrationAnchors(object.list = Col.tracing, dims = 1:40, k.anchor = 7)
Col.integrated <- IntegrateData(anchorset = Col.anchors, dims = 1:40)

DefaultAssay(Col.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Col.integrated <- ScaleData(Col.integrated, verbose = FALSE)
Col.integrated <- RunPCA(Col.integrated, npcs = 40, verbose = FALSE)
Col.integrated <- RunUMAP(Col.integrated, reduction = "pca", dims = 1:40)
Col.integrated <- FindNeighbors(Col.integrated, reduction = "pca", dims = 1:20)
Col.integrated <- FindClusters(Col.integrated, resolution = 0.7)
DimPlot(Col.integrated, group.by = "GFP", cols = c("grey90", "darkgreen"))

DefaultAssay(Col.integrated) <- "RNA"
FeaturePlot(Col.integrated, features = c("Vwf", "Pecam1", "Rgs5", "Lyve1", "Prox1", "Cdh5", "Col1a2", "Acta2", "EGFP"), cols = c("grey90", "blue"))

FeaturePlot(Col.integrated, features = c("EGFP"), cols = c("grey90", "blue"), min.cutoff = 0)
VlnPlot(Col.integrated, features = c("EGFP"), group.by = "seurat_clusters")


DimPlot(Col.integrated, group.by = "seurat_clusters") + theme(legend.position = "bottom") + labs(title = "Collagen Tracing Clusters")

Marker.cluster <- FindAllMarkers(Col.integrated, test.use = "bimod")

Marker.panel <- c("Acta2", "Postn", "Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Col4a1", "Cd74", "Lyz2", "Cd14", "Cxcr2", "Lgals3", "Napsa", "Vwf", "Eng", "Emcn", "Fabp4", "Pecam1", "Cdh5", "Tie1", "Lyve1", "Pdpn", "Prox1", "Ccl5", "Nkg7", "Ptprc", "Klre1", "Ctla4", "Icos", "Cd3e", "Lat", "Cd79a", "Cd79b", "Pax5", "Plp1", "Kcna1", "Cd59a", "Rgs5", "Pdgfrb", "Tagln", "Des", "Krt19", "Krt8", "Nkain4", "EGFP")

DotPlot(Col.integrated, features = Marker.panel) + coord_flip()

Ec.marker <- FindMarkers(Col.integrated, ident.1 = "3", ident.2 = "23")
Ec.marker$Gene <- rownames(Ec.marker)

AMI <- NULL
Age <- NULL
Sex <- NULL
Condition <- NULL

tbl = data.frame(Sample = factor(Sample, levels = Sample),
                 Sex = c("male", "male", "female", "female", "male", "male", "male", "male"),
                 Age = c("Old", "Old", "Old", "Old", "Young", "Young", "Young", "Young"),
                 AMI = c("AMI", "Hom", "AMI", "Hom", "AMI", "AMI", "Hom", "Hom"), 
                 Condition = c("Old AMI", "Old Hom", "Old AMI", "Old Hom", "Young AMI", "Young AMI", "Young Hom", "Young Hom"))

rownames(tbl) <- tbl$Sample

for(i in rownames(Col.integrated@meta.data)){
  S <- as.character(Col.integrated@meta.data[i, "orig.ident"])
  AMI[i] <- as.character(tbl[S, "AMI"])
  Age[i] <- as.character(tbl[S, "Age"])
  Sex[i] <- as.character(tbl[S, "Sex"])
  Condition[i] <- as.character(tbl[S, "Condition"])
}

Col.integrated$AMI <- AMI
Col.integrated$Age <- Age
Col.integrated$Sex <- Sex
Col.integrated$Condition <- Condition

save(Col.integrated, file = "Collagen_July_2020.Rds")
