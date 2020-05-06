setwd("STARmap_AllenVISp/")
library(Seurat)

starmap <- readRDS("data/seurat_objects/20180505_BY3_1kgenes.rds")
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

#Project on allen labels
i2 <- FindTransferAnchors(
  reference = allen,
  query = starmap,
  features = rownames(starmap),
  reduction = 'cca',
  reference.assay = 'RNA',
  query.assay = 'RNA'
)

refdata <- GetAssayData(
  object = allen,
  assay = 'RNA',
  slot = 'data'
)

imputation <- TransferData(
  anchorset = i2,
  refdata = refdata,
  weight.reduction = 'pca'
)

starmap[['ss2']] <- imputation
saveRDS(starmap, 'data/seurat_objects/20180505_BY3_1kgenes_imputed.rds')