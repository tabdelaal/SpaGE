setwd("osmFISH_AllenVISp/")
library(Seurat)
library(hdf5r)
library(methods)

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
mat <- t(mat)

osm_seurat <- CreateSeuratObject(counts = mat, project = 'osmFISH', assay = 'RNA', min.cells = -1, min.features = -1)
names(region) <- colnames(osm_seurat)
names(cluster) <- colnames(osm_seurat)
osm_seurat <- AddMetaData(osm_seurat, region, col.name = 'region')
osm_seurat <- AddMetaData(osm_seurat, cluster, col.name = 'cluster')
osm_seurat[['spatial']] <- CreateDimReducObject(embeddings = spatial, key = 'spatial', assay = 'RNA')
Idents(osm_seurat) <- 'region'
osm_seurat <- SubsetData(osm_seurat, ident.remove = c('Excluded','Internal Capsule Caudoputamen','White matter','Hippocampus','Ventricle'))
total.counts = colSums(x = as.matrix(osm_seurat@assays$RNA@counts))
osm_seurat <- NormalizeData(osm_seurat, scale.factor = median(x = total.counts))
osm_seurat <- ScaleData(osm_seurat)
saveRDS(object = osm_seurat, file = 'data/seurat_objects/osmFISH_Cortex.rds')