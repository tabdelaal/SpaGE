setwd("osmFISH_AllenSSp/")
library(liger)
library(hdf5r)
library(methods)

# allen VISp
allen <- read.table(file = "data/Allen_SSp/SSp_exons_matrix.csv",
                    row.names = 1, sep = ',', stringsAsFactors = FALSE, header = TRUE)
allen <- as.matrix(x = allen)

Genes_count = rowSums(allen > 0)
allen <- allen[Genes_count>=10,]

# osmFISH
osm <- H5File$new("data/osmFISH/osmFISH_SScortex_mouse_all_cells.loom")
mat <- osm[['matrix']][,]
colnames(mat) <- osm[['row_attrs']][['Gene']][]
rownames(mat) <- paste0('osm_', osm[['col_attrs']][['CellID']][])
x_dim <- osm[['col_attrs']][['X']][]
y_dim <- osm[['col_attrs']][['Y']][]
region <- osm[['col_attrs']][['Region']][]
cluster <- osm[['col_attrs']][['ClusterName']][]
osm$close_all()
spatial <- data.frame(spatial1 = x_dim, spatial2 = y_dim)
rownames(spatial) <- rownames(mat)
spatial <- as.matrix(spatial)
osmFISH <- t(mat)
osmFISH <- osmFISH[,!is.element(region,c('Excluded','Internal Capsule Caudoputamen','White matter','Hippocampus','Ventricle'))]

Gene_set <- intersect(rownames(osmFISH),rownames(allen))

#### New genes prediction
Ligerex <- createLiger(list(SMSC_FISH = osmFISH, SMSC_RNA = allen))
Ligerex <- normalize(Ligerex)
Ligerex@var.genes <- Gene_set
Ligerex <- scaleNotCenter(Ligerex)

# suggestK(Ligerex, k.test= seq(5,30,5)) # K = 20
# suggestLambda(Ligerex, k = 20) # Lambda = 20

Ligerex <- optimizeALS(Ligerex,k = 20, lambda = 20)
Ligerex <- quantileAlignSNF(Ligerex)

Imputation <- imputeKNN(Ligerex,reference = 'SMSC_RNA', queries = list('SMSC_FISH'), norm = TRUE, scale = FALSE)
new.genes <- c('Tesc', 'Pvrl3', 'Grm2')
Imp_New_genes <- matrix(0,nrow = length(new.genes),ncol = dim(Imputation@norm.data$SMSC_FISH)[2])
rownames(Imp_New_genes) <- new.genes
for(i in 1:length(new.genes)) {
  Imp_New_genes[new.genes[[i]],] = as.vector(Imputation@norm.data$SMSC_FISH[new.genes[i],])
}
write.csv(Imp_New_genes,file = 'Results/Liger_New_genes.csv')

# leave-one-out validation
genes.leaveout <- intersect(rownames(osmFISH),rownames(allen))
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(osmFISH)[2])
rownames(Imp_genes) <- genes.leaveout
colnames(Imp_genes) <- colnames(osmFISH)
NMF_time <- vector(mode= "numeric")
knn_time <- vector(mode= "numeric")

for(i in 1:length(genes.leaveout)) {
  print(i)
  start_time <- Sys.time()
  Ligerex.leaveout <- createLiger(list(SMSC_FISH = osmFISH[-which(rownames(osmFISH) %in% genes.leaveout[i]),], SMSC_RNA = allen[rownames(osmFISH),]))
  Ligerex.leaveout <- normalize(Ligerex.leaveout)
  Ligerex.leaveout@var.genes <- setdiff(Gene_set,genes.leaveout[i])
  Ligerex.leaveout <- scaleNotCenter(Ligerex.leaveout)
  if(dim(Ligerex.leaveout@norm.data$SMSC_FISH)[2]!=dim(osmFISH)[2]){
    next
  }
  Ligerex.leaveout <- optimizeALS(Ligerex.leaveout,k = 20, lambda = 20)
  Ligerex.leaveout <- quantileAlignSNF(Ligerex.leaveout)
  end_time <- Sys.time()
  NMF_time <- c(NMF_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  
  start_time <- Sys.time()
  Imputation <- imputeKNN(Ligerex.leaveout,reference = 'SMSC_RNA', queries = list('SMSC_FISH'), norm = TRUE, scale = FALSE, knn_k = 30)
  Imp_genes[genes.leaveout[[i]],] = as.vector(Imputation@norm.data$SMSC_FISH[genes.leaveout[i],])
  end_time <- Sys.time()
  knn_time <- c(knn_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
}
write.csv(Imp_genes,file = 'Results/Liger_LeaveOneOut.csv')
write.csv(NMF_time,file = 'Results/Liger_NMF_time.csv',row.names = FALSE)
write.csv(knn_time,file = 'Results/Liger_knn_time.csv',row.names = FALSE)