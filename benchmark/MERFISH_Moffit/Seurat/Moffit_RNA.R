setwd("MERFISH_Moffit/")
library(Seurat)

Moffit <- Read10X("data/Moffit_RNA/GSE113576/")

Mo <- CreateSeuratObject(counts = Moffit, project = 'POR', min.cells = 10)
Mo <- NormalizeData(object = Mo)
Mo <- FindVariableFeatures(object = Mo, nfeatures = 2000)
Mo <- ScaleData(object = Mo)
Mo <- RunPCA(object = Mo, npcs = 50, verbose = FALSE)
Mo <- RunUMAP(object = Mo, dims = 1:50, nneighbors = 5)
saveRDS(object = Mo, file = paste0("data/seurat_objects/","Moffit_RNA.rds"))