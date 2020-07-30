setwd("osmFISH_AllenSSp/")
library(Seurat)

allen <- read.table(file = "data/Allen_SSp/SSp_exons_matrix.csv",
                    row.names = 1,sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen <- as.matrix(x = allen)

meta.data <- read.csv(file = "data/Allen_SSp/AllenSSp_metadata.csv",
                      row.names = 1, stringsAsFactors = FALSE)

al <- CreateSeuratObject(counts = allen, project = 'SSp', min.cells = 10)
ok.cells <- colnames(al[,meta.data$class_label != 'Exclude'])
al <- al[, ok.cells]
al <- NormalizeData(object = al)
al <- FindVariableFeatures(object = al, nfeatures = 2000)
al <- ScaleData(object = al)
al <- RunPCA(object = al, npcs = 50, verbose = FALSE)
al <- RunUMAP(object = al, dims = 1:50, nneighbors = 5)
saveRDS(object = al, file = paste0("data/seurat_objects/","allen_brain_SSp.rds"))
