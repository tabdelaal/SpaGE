setwd("osmFISH_AllenSSp/")
library(Seurat)

osmFISH <- readRDS("data/seurat_objects/osmFISH_Cortex.rds")
allen <- readRDS("data/seurat_objects/allen_brain_SSp.rds")

#Project on allen labels
i2 <- FindTransferAnchors(
  reference = allen,
  query = osmFISH,
  features = rownames(osmFISH),
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

osmFISH[['ss2']] <- imputation
saveRDS(osmFISH, 'data/seurat_objects/osmFISH_Cortex_imputed.rds')