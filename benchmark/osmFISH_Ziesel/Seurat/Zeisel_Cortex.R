setwd("osmFISH_Ziesel/")
library(Seurat)

Zeisel <- read.delim(file = "data/Zeisel/expression_mRNA_17-Aug-2014.txt",header = FALSE)

meta.data <- Zeisel[1:10,]
rownames(meta.data) <- meta.data[,2]
meta.data <- meta.data[,-c(1,2)]
meta.data <- as.data.frame(t(meta.data))

Zeisel <- Zeisel[12:19983,]
gene_names <- Zeisel[,1]
Zeisel <- Zeisel[,-c(1,2)]
Zeisel <- apply(Zeisel,2,as.numeric)
Zeisel <- as.matrix(Zeisel)
rownames(Zeisel) <- gene_names
Zeisel <- Zeisel[,meta.data$tissue=='sscortex']
meta.data <- meta.data[meta.data$tissue=='sscortex',]

colnames(Zeisel) <- paste0('SMSC_',c(1:dim(Zeisel)[2]))
rownames(meta.data) <- paste0('SMSC_',c(1:dim(meta.data)[1]))

Zeisel <- CreateSeuratObject(counts = Zeisel, project = 'SMSC', meta.data = meta.data, min.cells = 10)
Zeisel <- NormalizeData(object = Zeisel)
Zeisel <- FindVariableFeatures(object = Zeisel, nfeatures = 2000)
Zeisel <- ScaleData(object = Zeisel)
Zeisel <- RunPCA(object = Zeisel, npcs = 50, verbose = FALSE)
Zeisel <- RunUMAP(object = Zeisel, dims = 1:50, nneighbors = 5)
saveRDS(object = Zeisel, file = paste0("data/seurat_objects/","Zeisel_SMSC.rds"))