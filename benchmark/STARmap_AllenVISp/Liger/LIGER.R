setwd("STARmap_AllenVISp/")
library(liger)
library(Seurat)
library(ggplot2)

# allen VISp
allen <- read.table(file = "data/Allen_VISp/mouse_VISp_2018-06-14_exon-matrix.csv",
                    row.names = 1, sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen <- as.matrix(x = allen)
genes <- read.table(file = "data/Allen_VISp/mouse_VISp_2018-06-14_genes-rows.csv",
                    sep = ',', stringsAsFactors = FALSE, header = TRUE)
rownames(x = allen) <- make.unique(names = genes$gene_symbol)
meta.data <- read.csv(file = "data/Allen_VISp/mouse_VISp_2018-06-14_samples-columns.csv",
                      row.names = 1, stringsAsFactors = FALSE)

Genes_count = rowSums(allen > 0)
allen <- allen[Genes_count>=10,]

low.q.cells <- rownames(x = meta.data[meta.data$class %in% c('Low Quality', 'No Class'), ])
ok.cells <- rownames(x = meta.data)[!(rownames(x = meta.data) %in% low.q.cells)]
allen <-allen[,ok.cells]

# STARmap
counts <- read.table(file = "data/Starmap/visual_1020/20180505_BY3_1kgenes/cell_barcode_count.csv",
                     sep = ",", stringsAsFactors = FALSE)
gene.names <- read.table(file = "data/Starmap/visual_1020/20180505_BY3_1kgenes/genes.csv",
                         sep = ",", stringsAsFactors = FALSE)

qhulls <- read.table(file = "data/Starmap/visual_1020/20180505_BY3_1kgenes/qhulls.tsv",
                     sep = '\t', col.names = c('cell', 'y', 'x'), stringsAsFactors = FALSE)
centroids <- read.table(file = "data/Starmap/visual_1020/20180505_BY3_1kgenes/centroids.tsv",
                        sep = "\t", col.names = c("spatial2", "spatial1"), stringsAsFactors = FALSE)
colnames(counts) <- gene.names$V1
rownames(counts) <- paste0('starmap', seq(1:nrow(counts)))
counts <- as.matrix(counts)
rownames(centroids) <- rownames(counts)
centroids <- as.matrix(centroids)
counts <- t(counts)

Gene_set <- intersect(rownames(counts),rownames(allen))

#### New genes prediction
Ligerex <- createLiger(list(STARmap = counts, AllenVISp = allen))
Ligerex <- normalize(Ligerex)
Ligerex@var.genes <- Gene_set
Ligerex <- scaleNotCenter(Ligerex)

# suggestK(Ligerex) # K = 25
# suggestLambda(Ligerex, k = 25) # Lambda = 35

Ligerex <- optimizeALS(Ligerex,k = 25, lambda = 35)
Ligerex <- quantileAlignSNF(Ligerex)

Imputation <- imputeKNN(Ligerex,reference = 'AllenVISp', queries = list('STARmap'), norm = TRUE, scale = FALSE)

new.genes <- c('Tesc', 'Pvrl3', 'Sox10', 'Grm2', 'Tcrb')
Imp_New_genes <- matrix(0,nrow= length(new.genes),ncol = dim(Imputation@norm.data$STARmap)[2])
rownames(Imp_New_genes) <- new.genes

for (i in c(1:length(new.genes))){
  Imp_New_genes[new.genes[[i]],] = as.vector(Imputation@norm.data$STARmap[new.genes[i],])
}
write.csv(Imp_New_genes,file = 'Results/Liger_New_genes.csv')

# leave-one-out validation
genes.leaveout <- Gene_set
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(counts)[2])
rownames(Imp_genes) <- genes.leaveout
colnames(Imp_genes) <- colnames(counts)
NMF_time <- vector(mode= "numeric")
knn_time <- vector(mode= "numeric")

for(i in 1:length(genes.leaveout)) {
  print(i)
  start_time <- Sys.time()
  Ligerex.leaveout <- createLiger(list(STARmap = counts[-which(rownames(counts) %in% genes.leaveout[i]),], AllenVISp = allen))
  Ligerex.leaveout <- normalize(Ligerex.leaveout)
  Ligerex.leaveout@var.genes <- setdiff(Gene_set,genes.leaveout[i])
  Ligerex.leaveout <- scaleNotCenter(Ligerex.leaveout)
  Ligerex.leaveout <- optimizeALS(Ligerex.leaveout,k = 25, lambda = 35)
  Ligerex.leaveout <- quantileAlignSNF(Ligerex.leaveout)
  end_time <- Sys.time()
  NMF_time <- c(NMF_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  
  start_time <- Sys.time()
  Imputation <- imputeKNN(Ligerex.leaveout,reference = 'AllenVISp', queries = list('STARmap'), norm = TRUE, scale = FALSE)
  Imp_genes[genes.leaveout[[i]],] = as.vector(Imputation@norm.data$STARmap[genes.leaveout[i],])
  end_time <- Sys.time()
  knn_time <- c(knn_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
}
write.csv(Imp_genes,file = 'Results/Liger_LeaveOneOut.csv')
write.csv(NMF_time,file = 'Results/Liger_NMF_time.csv',row.names = FALSE)
write.csv(knn_time,file = 'Results/Liger_knn_time.csv',row.names = FALSE)