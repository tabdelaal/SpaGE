setwd("osmFISH_AllenVISp/")
library(Seurat)

allen <- read.table(file = "data/Allen_VISp/mouse_VISp_2018-06-14_exon-matrix.csv",
                    row.names = 1,sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen <- as.matrix(x = allen)
genes <- read.table(file = "data/Allen_VISp/mouse_VISp_2018-06-14_genes-rows.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
rownames(x = allen) <- make.unique(names = genes$gene_symbol)
meta.data <- read.csv(file = "data/Allen_VISp/mouse_VISp_2018-06-14_samples-columns.csv",
                      row.names = 1, stringsAsFactors = FALSE)

al <- CreateSeuratObject(counts = allen, project = 'VISp', meta.data = meta.data, min.cells = 10)
low.q.cells <- rownames(x = meta.data[meta.data$class %in% c('Low Quality', 'No Class'), ])
ok.cells <- rownames(x = meta.data)[!(rownames(x = meta.data) %in% low.q.cells)]
al <- al[, ok.cells]
al <- NormalizeData(object = al)
al <- FindVariableFeatures(object = al, nfeatures = 2000)
al <- ScaleData(object = al)
al <- RunPCA(object = al, npcs = 50, verbose = FALSE)
al <- RunUMAP(object = al, dims = 1:50, nneighbors = 5)
saveRDS(object = al, file = paste0("data/seurat_objects/","allen_brain.rds"))