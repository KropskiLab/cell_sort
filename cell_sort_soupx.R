library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(knitr)
gc
set.seed(1234)

parent_directory <- "~/Desktop/8745/"

directories <- c(
  "8745_1", "8745_2"
)


for (dir in directories) {
  filtered_file <- file.path(parent_directory, dir, 'filtered_feature_bc_matrix.h5')
  raw_file <- file.path(parent_directory, dir, 'raw_feature_bc_matrix.h5')
  filt.matrix <- Read10X_h5(filtered_file,use.names = T)
  raw.matrix <- Read10X_h5(raw_file,use.names = T)
  srat  <- CreateSeuratObject(counts = filt.matrix)
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  soup.channel
  srat    <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  head(meta)
  soup.channel  <- autoEstCont(soup.channel)
  head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 20)
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  DropletUtils:::write10xCounts(file.path(parent_directory, dir, "soupx"), adj.matrix)
}