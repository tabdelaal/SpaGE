setwd("STARmap_AllenVISp/")
library(Seurat)
library(ggplot2)

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")
starmap.imputed <- readRDS("data/seurat_objects/20180505_BY3_1kgenes_imputed.rds")
allen <- readRDS("data/seurat_objects/allen_brain.rds")

# remove HPC from starmap
class_labels <- read.table(
  file = "data/Starmap/visual_1020/20180505_BY3_1kgenes/class_labels.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)

class_labels$cellname <- paste0('starmap', rownames(class_labels))

class_labels$ClusterName <- ifelse(is.na(class_labels$ClusterName), 'Other', class_labels$ClusterName)

hpc <- class_labels[class_labels$ClusterName == 'HPC', ]$cellname

accept.cells <- setdiff(colnames(starmap), hpc)

starmap <- starmap[, accept.cells]

starmap@misc$spatial <- starmap@misc$spatial[starmap@misc$spatial$cell %in% accept.cells, ]

genes.leaveout <- intersect(rownames(starmap),rownames(allen))
Imp_genes <- matrix(0,nrow = length(genes.leaveout),ncol = dim(starmap@assays$RNA)[2])
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
  imputed.ss2 <- run_imputation(ref.obj = allen, query.obj = starmap, feature.remove = genes.leaveout[[i]])
  starmap[['ss2']] <- imputed.ss2[, colnames(starmap)][['seq']]
  Imp_genes[genes.leaveout[[i]],] = as.vector(starmap@assays$ss2[genes.leaveout[i],])
}
write.csv(Imp_genes,file = 'Results/Seurat_LeaveOneOut.csv')
write.csv(anchor_time,file = 'Results/Seurat_anchor_time.csv',row.names = FALSE)
write.csv(Transfer_time,file = 'Results/Seurat_transfer_time.csv',row.names = FALSE)

# show genes not in the starmap dataset
DefaultAssay(starmap.imputed) <- "ss2"
new.genes <- c('Tesc', 'Pvrl3', 'Sox10', 'Grm2', 'Tcrb')
Imp_New_genes <- matrix(0,nrow = length(new.genes),ncol = dim(starmap.imputed@assays$ss2)[2])
rownames(Imp_New_genes) <- new.genes
for(i in 1:length(new.genes)) {
  Imp_New_genes[new.genes[[i]],] = as.vector(starmap.imputed@assays$ss2[new.genes[i],])
}
write.csv(Imp_New_genes,file = 'Results/Seurat_New_genes.csv')