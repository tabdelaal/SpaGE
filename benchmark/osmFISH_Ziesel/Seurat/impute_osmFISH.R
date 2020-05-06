setwd("osmFISH_Ziesel/")
library(Seurat)

osmFISH <- readRDS("data/seurat_objects/osmFISH_Cortex.rds")
Zeisel <- readRDS("data/seurat_objects/Zeisel_SMSC.rds")

#Project on Zeisel labels
i2 <- FindTransferAnchors(
  reference = Zeisel,
  query = osmFISH,
  features = rownames(osmFISH),
  reduction = 'cca',
  reference.assay = 'RNA',
  query.assay = 'RNA'
)

refdata <- GetAssayData(
  object = Zeisel,
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