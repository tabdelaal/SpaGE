setwd("MERFISH_Moffit/")
library(liger)
library(Seurat)
library(ggplot2)

# Moffit RNA
Moffit <- Read10X("data/Moffit_RNA/GSE113576/")

Moffit <- as.matrix(Moffit)
Genes_count = rowSums(Moffit > 0)
Moffit <- Moffit[Genes_count>=10,]

# MERFISH
MERFISH <- read.csv(file = "data/MERFISH/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv", header = TRUE)
MERFISH_1 <- MERFISH[MERFISH['Animal_ID']==1,]
MERFISH_1 <- MERFISH_1[MERFISH_1['Cell_class']!='Ambiguous',]
MERFISH_meta <- MERFISH_1[,c(1:9)]
MERFISH_data <- MERFISH_1[,c(10:170)]
drops <- c('Blank_1','Blank_2','Blank_3','Blank_4','Blank_5','Fos')
MERFISH_data <- MERFISH_data[ , !(colnames(MERFISH_data) %in% drops)]
MERFISH_data <- t(MERFISH_data)

Gene_set <- intersect(rownames(MERFISH_data),rownames(Moffit))

#### New genes prediction
Ligerex <- createLiger(list(MERFISH = MERFISH_data, Moffit_RNA = Moffit))
Ligerex <- normalize(Ligerex)
Ligerex@var.genes <- Gene_set
Ligerex <- scaleNotCenter(Ligerex)

# suggestK(Ligerex) # K = 25
# suggestLambda(Ligerex, k = 25)

Ligerex <- optimizeALS(Ligerex,k = 25)
Ligerex <- quantileAlignSNF(Ligerex)

# leave-one-out validation
genes.leaveout <- Gene_set
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(MERFISH_data)[2])
rownames(Imp_genes) <- genes.leaveout
colnames(Imp_genes) <- colnames(MERFISH_data)
NMF_time <- vector(mode= "numeric")
knn_time <- vector(mode= "numeric")

for(i in 1:length(genes.leaveout)) {
  print(i)
  start_time <- Sys.time()
  Ligerex.leaveout <- createLiger(list(MERFISH = MERFISH_data[-which(rownames(MERFISH_data) %in% genes.leaveout[i]),], Moffit_RNA = Moffit))
  Ligerex.leaveout <- normalize(Ligerex.leaveout)
  Ligerex.leaveout@var.genes <- setdiff(Gene_set,genes.leaveout[i])
  Ligerex.leaveout <- scaleNotCenter(Ligerex.leaveout)
  Ligerex.leaveout <- optimizeALS(Ligerex.leaveout,k = 25)
  Ligerex.leaveout <- quantileAlignSNF(Ligerex.leaveout)
  end_time <- Sys.time()
  NMF_time <- c(NMF_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  
  start_time <- Sys.time()
  Imputation <- imputeKNN(Ligerex.leaveout,reference = 'Moffit_RNA', queries = list('MERFISH'), norm = TRUE, scale = FALSE)
  Imp_genes[genes.leaveout[[i]],] = as.vector(Imputation@norm.data$MERFISH[genes.leaveout[i],])
  end_time <- Sys.time()
  knn_time <- c(knn_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
}
write.csv(Imp_genes,file = 'Results/Liger_LeaveOneOut.csv')
write.csv(NMF_time,file = 'Results/Liger_NMF_time.csv',row.names = FALSE)
write.csv(knn_time,file = 'Results/Liger_knn_time.csv',row.names = FALSE)