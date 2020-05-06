setwd("STARmap_AllenVISp/")
library(Seurat)
library(ggplot2)

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")
starmap.imputed <- readRDS("data/seurat_objects/20180505_BY3_1kgenes_imputed.rds")
allen <- readRDS("data/seurat_objects/allen_brain.rds")

DefaultAssay(starmap.imputed) <- "RNA"
genes.leaveout <- intersect(rownames(starmap),rownames(allen))

# Original
for(i in 1:length(genes.leaveout)) {
  if (genes.leaveout[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap.imputed, genes.leaveout[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/Original/', genes.leaveout[[i]],'.pdf'))
}

# Seurat_Predicted
Seurat_Predicted <- read.csv(file = 'Results/Seurat_LeaveOneOut.csv',header = TRUE, check.names = FALSE, row.names = 1)
Seurat_Predicted <- Seurat_Predicted[genes.leaveout,]
colnames(Seurat_Predicted) <- colnames(starmap.imputed@assays$ss2[,])
starmap.imputed[['ss2']] <- CreateAssayObject(data = as.matrix(Seurat_Predicted))
DefaultAssay(starmap.imputed) <- 'ss2'
for(i in 1:length(genes.leaveout)) {
  if (genes.leaveout[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap.imputed, genes.leaveout[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/Seurat_Predicted/', genes.leaveout[[i]],'.pdf'))
}

# SpaGE_Predicted
SpaGE_Predicted <- read.csv(file = 'Results/SpaGE_LeaveOneOut.csv',header = TRUE, check.names = FALSE, row.names = 1)
SpaGE_Predicted <- t(SpaGE_Predicted)
SpaGE_Predicted <- SpaGE_Predicted[genes.leaveout,]
colnames(SpaGE_Predicted) <- colnames(starmap@assays$RNA[,])
starmap[['ss2']] <- CreateAssayObject(data = as.matrix(SpaGE_Predicted))
DefaultAssay(starmap) <- 'ss2'
for(i in 1:length(genes.leaveout)) {
  if (genes.leaveout[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap, genes.leaveout[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/SpaGE_Predicted/', genes.leaveout[[i]],'.pdf'))
}

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")

# Liger_Predicted
Liger_Predicted <- read.csv(file = 'Results/Liger_LeaveOneOut.csv',header = TRUE, check.names = FALSE, row.names = 1)
Liger_Predicted <- Liger_Predicted[genes.leaveout,]
colnames(Liger_Predicted) <- colnames(starmap@assays$RNA[,])
starmap[['ss2']] <- CreateAssayObject(data = as.matrix(Liger_Predicted))
DefaultAssay(starmap) <- 'ss2'
for(i in 1:length(genes.leaveout)) {
  if (genes.leaveout[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap, genes.leaveout[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/Liger_Predicted/', genes.leaveout[[i]],'.pdf'))
}

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")

# gimVI_Predicted
gimVI_Predicted <- read.csv(file = 'Results/gimVI_LeaveOneOut.csv',header = TRUE, check.names = FALSE, row.names = 1)
gimVI_Predicted <- t(gimVI_Predicted)
gimVI_Predicted <- gimVI_Predicted[toupper(genes.leaveout),]
colnames(gimVI_Predicted) <- colnames(starmap@assays$RNA[,])
rownames(gimVI_Predicted) <- genes.leaveout
starmap[['ss2']] <- CreateAssayObject(data = as.matrix(gimVI_Predicted))
DefaultAssay(starmap) <- 'ss2'
for(i in 1:length(genes.leaveout)) {
  if (genes.leaveout[[i]] %in% c('Bsg', 'Rab3c')) {
    qcut = 'q40'
  } else {
    qcut = 0
  }
  p <- PolyFeaturePlot(starmap, genes.leaveout[[i]], flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/gimVI_Predicted/', genes.leaveout[[i]],'.pdf'))
}


### NEw genes 
# show genes not in the starmap dataset
new.genes <- c('Tesc', 'Pvrl3', 'Sox10', 'Grm2', 'Tcrb')
starmap.imputed <- readRDS("data/seurat_objects/20180505_BY3_1kgenes_imputed.rds")
DefaultAssay(starmap.imputed) <- "ss2"

for(i in new.genes) {
  p <- PolyFeaturePlot(starmap.imputed, i, flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/Seurat_Predicted/New_', i,'.pdf'))
}

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")
#SpaGE
SpaGE_New <- read.csv(file = 'Results/SpaGE_New_genes.csv',header = TRUE, check.names = FALSE, row.names = 1)
SpaGE_New <- t(SpaGE_New)
new.genes <- c('Tesc', 'Pvrl3', 'Sox10', 'Grm2', 'Tcrb','Ttyh2','Cldn11','Tmem88b')
SpaGE_New2 <- SpaGE_New[new.genes,]
colnames(SpaGE_New2) <- colnames(starmap@assays$RNA[,])
starmap[['ss2']] <- CreateAssayObject(data = as.matrix(SpaGE_New2))
DefaultAssay(starmap) <- 'ss2'
for(i in new.genes) {
  p <- PolyFeaturePlot(starmap, i, flip.coords = TRUE, min.cutoff = 0, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/SpaGE_Predicted/New_', i,'.pdf'))
}

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")
# Liger
Liger_New <- read.csv(file = 'Results/Liger_New_genes.csv',header = TRUE, check.names = FALSE, row.names = 1)
Liger_New <- Liger_New[new.genes,]
colnames(Liger_New) <- colnames(starmap@assays$RNA[,])
starmap[['ss2']] <- CreateAssayObject(data = as.matrix(Liger_New))
DefaultAssay(starmap) <- 'ss2'
for(i in new.genes) {
  p <- PolyFeaturePlot(starmap, i, flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/Liger_Predicted/New_', i,'.pdf'))
}

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")
# gimVI
gimVI_New <- read.csv(file = 'Results/gimVI_New_genes.csv',header = TRUE, check.names = FALSE, row.names = 1)
gimVI_New <- t(gimVI_New)
gimVI_New <- gimVI_New[toupper(new.genes),]
colnames(gimVI_New) <- colnames(starmap@assays$RNA[,])
rownames(gimVI_New) <- new.genes
starmap[['ss2']] <- CreateAssayObject(data = as.matrix(gimVI_New))
DefaultAssay(starmap) <- 'ss2'
for(i in new.genes) {
  p <- PolyFeaturePlot(starmap, i, flip.coords = TRUE, min.cutoff = qcut, max.cutoff = 'q99') + SpatialTheme()
  ggsave(plot = p, filename = paste0('Figures/gimVI_Predicted/New_', i,'.pdf'))
}