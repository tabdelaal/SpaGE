setwd("osmFISH_AllenSSp/")
library(Seurat)

allen <- read.table(file = "data/Allen_SSp/SSp_exons_matrix.csv",
                    row.names = 1,sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen <- as.matrix(x = allen)

al <- CreateSeuratObject(counts = allen, project = 'SSp', min.cells = 10)
al <- NormalizeData(object = al)
al <- FindVariableFeatures(object = al, nfeatures = 2000)
al <- ScaleData(object = al)
al <- RunPCA(object = al, npcs = 50, verbose = FALSE)
al <- RunUMAP(object = al, dims = 1:50, nneighbors = 5)
saveRDS(object = al, file = paste0("data/seurat_objects/","allen_brain_SSp.rds"))