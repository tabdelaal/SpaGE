setwd("STARmap_AllenVISp/")
library(Seurat)
library(Matrix)

read_data <- function(base_path, project) {
  counts <- read.table(
    file = paste0(base_path, "cell_barcode_count.csv"),
    sep = ",",
    stringsAsFactors = FALSE
  )
  gene.names <- read.table(
    file = paste0(base_path, "genes.csv"),
    sep = ",",
    stringsAsFactors = FALSE
  )
  qhulls <- read.table(
    file = paste0(base_path, "qhulls.tsv"),
    sep = '\t',
    col.names = c('cell', 'y', 'x'),
    stringsAsFactors = FALSE
  )
  centroids <- read.table(
    file = paste0(base_path, "centroids.tsv"),
    sep = "\t",
    col.names = c("spatial2", "spatial1"),
    stringsAsFactors = FALSE
  )
  colnames(x = counts) <- gene.names$V1
  rownames(x = counts) <- paste0('starmap', seq(1:nrow(x = counts)))
  counts <- as.matrix(x = counts)
  rownames(x = centroids) <- rownames(x = counts)
  centroids <- as.matrix(x = centroids)
  total.counts = rowSums(x = counts)
  
  obj <- CreateSeuratObject(
    counts = t(x = counts),
    project = project,
    min.cells = -1,
    min.features = -1
  )
  obj <- NormalizeData(
    object = obj,
    scale.factor = median(x = total.counts)
  )
  obj <- ScaleData(
    object = obj,
    features = rownames(x = obj)
  )
  obj <- RunPCA(
    object = obj,
    features = rownames(x = obj),
    npcs = 30
  )
  obj <- RunUMAP(
    object = obj,
    dims = 1:30
  )
  obj[['spatial']] <- CreateDimReducObject(
    embeddings <- centroids, 
    assay = "RNA",
    key = 'spatial'
  )
  qhulls$cell <- paste0('starmap', qhulls$cell)
  obj@misc[['spatial']] <- qhulls
  return(obj)
}

# experiments <- c("visual_1020/20180505_BY3_1kgenes/", "visual_1020/20180410-BY3_1kgenes/")
experiments <- c("data/Starmap/visual_1020/20180505_BY3_1kgenes/")

for(i in 1:length(experiments)) {
  project.name <- unlist(x = strsplit(x = experiments[[i]], "/"))[[4]]
  dat <- read_data(base_path =  experiments[[i]],project = project.name)
  saveRDS(object = dat, file = paste0("data/seurat_objects/", project.name, ".rds"))
}