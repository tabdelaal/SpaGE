setwd("osmFISH_AllenVISp/")
library(Seurat)
library(ggplot2)

osmFISH <- readRDS("data/seurat_objects/osmFISH_Cortex.rds")
osmFISH.imputed <- readRDS("data/seurat_objects/osmFISH_Cortex_imputed.rds")
allen <- readRDS("data/seurat_objects/allen_brain.rds")

genes.leaveout <- intersect(rownames(osmFISH),rownames(allen))
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(osmFISH@assays$RNA)[2])
rownames(Imp_genes) <- genes.leaveout
anchor_time <- vector(mode= "numeric")
Transfer_time <- vector(mode= "numeric")

run_imputation <- function(ref.obj, query.obj, feature.remove) {
  message(paste0('removing ', feature.remove))
  features <- setdiff(rownames(query.obj), feature.remove)
  DefaultAssay(ref.obj) <- 'RNA'
  DefaultAssay(query.obj) <- 'RNA'
  
  start_time <- Sys.time()
  anchors <- FindTransferAnchors(
    reference = ref.obj,
    query = query.obj,
    features = features,
    dims = 1:30,
    reduction = 'cca'
  )
  end_time <- Sys.time()
  anchor_time <<- c(anchor_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  
  refdata <- GetAssayData(
    object = ref.obj,
    assay = 'RNA',
    slot = 'data'
  )
  
  start_time <- Sys.time()
  imputation <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    weight.reduction = 'pca'
  )
  query.obj[['seq']] <- imputation
  end_time <- Sys.time()
  Transfer_time <<- c(Transfer_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  return(query.obj)
}

for(i in 1:length(genes.leaveout)) {
  imputed.ss2 <- run_imputation(ref.obj = allen, query.obj = osmFISH, feature.remove = genes.leaveout[[i]])
  osmFISH[['ss2']] <- imputed.ss2[, colnames(osmFISH)][['seq']]
  Imp_genes[genes.leaveout[[i]],] = as.vector(osmFISH@assays$ss2[genes.leaveout[i],])
}
write.csv(Imp_genes,file = 'Results/Seurat_LeaveOneOut.csv')
write.csv(anchor_time,file = 'Results/Seurat_anchor_time.csv',row.names = FALSE)
write.csv(Transfer_time,file = 'Results/Seurat_transfer_time.csv',row.names = FALSE)

# show genes not in the osmFISH dataset
DefaultAssay(osmFISH.imputed) <- "ss2"
new.genes <- c('Tesc', 'Pvrl3', 'Grm2')
Imp_New_genes <- matrix(0,nrow = length(new.genes),ncol = dim(osmFISH.imputed@assays$ss2)[2])
rownames(Imp_New_genes) <- new.genes
for(i in 1:length(new.genes)) {
  Imp_New_genes[new.genes[[i]],] = as.vector(osmFISH.imputed@assays$ss2[new.genes[i],])
}
write.csv(Imp_New_genes,file = 'Results/Seurat_New_genes.csv')