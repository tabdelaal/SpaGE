setwd("MERFISH_Moffit/")
library(Seurat)
library(ggplot2)

MERFISH <- readRDS("data/seurat_objects/MERFISH.rds")
Moffit <- readRDS("data/seurat_objects/Moffit_RNA.rds")

genes.leaveout <- intersect(rownames(MERFISH),rownames(Moffit))
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(MERFISH@assays$RNA)[2])
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
  anchor_time <- c(anchor_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  
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
  Transfer_time <- c(Transfer_time,as.numeric(difftime(end_time,start_time,units = 'secs')))
  return(query.obj)
}

for(i in 1:length(genes.leaveout)) {
  imputed.ss2 <- run_imputation(ref.obj = Moffit, query.obj = MERFISH, feature.remove = genes.leaveout[[i]])
  MERFISH[['ss2']] <- imputed.ss2[, colnames(MERFISH)][['seq']]
  Imp_genes[genes.leaveout[[i]],] = as.vector(MERFISH@assays$ss2[genes.leaveout[i],])
}
write.csv(Imp_genes,file = 'Results/Seurat_LeaveOneOut.csv')
write.csv(anchor_time,file = 'Results/Seurat_anchor_time.csv',row.names = FALSE)
write.csv(Transfer_time,file = 'Results/Seurat_transfer_time.csv',row.names = FALSE)