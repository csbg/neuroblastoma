
nb_rna <- 
  readRDS(
    file = "data_generated/rna_decontaminated.rds")
nb_rna_meta <- 
  readRDS(
    file = "data_generated/metadata.rds")

# rna_integrated_monocle.rds
# sce <- SingleCellExperiment(list(counts=counts))
# exprs(nb_rna)
# pData(nb_rna) %>% 
#   data.frame() %>% 
#   View()
# reducedDims(nb_rna)[['PCA']]
# reducedDims(nb_rna)[["Aligned"]]

colsum_nb_rna <- Matrix::colSums(exprs(nb_rna))
sampled_cell_isd <- sample(colsum_nb_rna,73872)
mat <- exprs(nb_rna)
mat1 <- mat[, which(colnames(mat) %in% names(sampled_cell_isd))]
mat2 <- nb_rna_meta[which(colnames(mat) %in% names(sampled_cell_isd)),] 

# sce <- 
#   SingleCellExperiment(
#   list(counts=exprs(nb_rna)),
#   metadata = pData(nb_rna)
#   )


sce <- 
  SingleCellExperiment(
    list(counts=mat1),
    metadata =mat2
  )


sce <-
  SingleCellExperiment(
    list(counts=exprs(nb_rna)),
    metadata = nb_rna_meta
  )

col_order <- colnames(sce)

pca <- reducedDims(nb_rna)[['PCA']]
pca <- pca[which(rownames(pca) %in% names(sampled_cell_isd)),]
pca <- pca[order(factor(rownames(pca), levels=col_order)),]

umap <- reducedDims(nb_rna)[['UMAP']]
umap <- umap[which(rownames(umap) %in% names(sampled_cell_isd)),]
umap <- umap[order(factor(rownames(umap), levels=col_order)),]

tsne <- reducedDims(nb_rna)[['tSNE']]
tsne <- tsne[which(rownames(tsne) %in% names(sampled_cell_isd)),]
tsne <- tsne[order(factor(rownames(tsne), levels=col_order)),]

aligned <- reducedDims(nb_rna)[['Aligned']]
aligned <- aligned[which(rownames(aligned) %in% names(sampled_cell_isd)), ]
aligned <- aligned[order(factor(rownames(aligned), levels=col_order)),]

# reducedDim(sce, "PCA") <- pca         #reducedDims(nb_rna)[['PCA']]
# reducedDim(sce, "UMAP") <- umap       #reducedDims(nb_rna)[['UMAP']]
# reducedDim(sce, "tSNE") <- tsne       #reducedDims(nb_rna)[['tSNE']]
# reducedDim(sce, "Aligned") <- aligned #reducedDims(nb_rna)[['Aligned']]

reducedDim(sce, "PCA") <- reducedDims(nb_rna)[['PCA']]
reducedDim(sce, "UMAP") <- reducedDims(nb_rna)[['UMAP']]
reducedDim(sce, "tSNE") <- reducedDims(nb_rna)[['tSNE']]
reducedDim(sce, "Aligned") <- reducedDims(nb_rna)[['Aligned']]

## Save H5AD file

sceasy::convertFormat(
  Cluster_myeloid, 
  main_layer = 'counts', 
  from="seurat", 
  to="anndata",
  outFile='myeloid.h5ad')



library(zellkonverter)
SCE2AnnData(sce)

sce_rna <- 
  SCE2AnnData(
    sce, 
    X_name = 'counts', 
    #obs = TRUE,
    skip_assays = FALSE)

anndata::write_h5ad(
  sce_rna,
  filename  = "sce_rna_wes_sampled.h5ad",
  compression = "gzip",
  compression_opts=9,
  as_dense = "X"
)

myeloid.only <- sce[, sce@metadata$cellont_abbr == "M"]
ncol(counts(myeloid.only))
colData(myeloid.only)

mat3 <- mat[, which(colnames(mat) %in% colnames(myeloid.only))]
mat4 <- nb_rna_meta[which(colnames(mat) %in% colnames(myeloid.only)),] 

sce_myeloid <- 
  SingleCellExperiment(
    list(counts = mat3),
    metadata = mat4
  )

col_order <- colnames(sce_myeloid)

pca <- reducedDims(nb_rna)[['PCA']]
pca <- pca[which(rownames(pca) %in% colnames(myeloid.only)),]
pca <- pca[order(factor(rownames(pca), levels=col_order)),]

umap <- reducedDims(nb_rna)[['UMAP']]
umap <- umap[which(rownames(umap) %in% colnames(myeloid.only)),]
umap <- umap[order(factor(rownames(umap), levels=col_order)),]

tsne <- reducedDims(nb_rna)[['tSNE']]
tsne <- tsne[which(rownames(tsne) %in% colnames(myeloid.only)),]
tsne <- tsne[order(factor(rownames(tsne), levels=col_order)),]

aligned <- reducedDims(nb_rna)[['Aligned']]
aligned <- aligned[which(rownames(aligned) %in% colnames(myeloid.only)), ]
aligned <- aligned[order(factor(rownames(aligned), levels=col_order)),]

reducedDim(sce_myeloid, "PCA") <- pca          #reducedDims(nb_rna)[['PCA']]
reducedDim(sce_myeloid, "UMAP") <- umap        #reducedDims(nb_rna)[['UMAP']]
reducedDim(sce_myeloid, "tSNE") <- tsne        #reducedDims(nb_rna)[['tSNE']]
reducedDim(sce_myeloid, "Aligned") <- aligned  #reducedDims(nb_rna)[['Aligned']]

library(zellkonverter)

sce_rna <- 
  SCE2AnnData(
    sce_myeloid, 
    X_name = 'counts', 
    #obs = TRUE,
    skip_assays = FALSE)

anndata::write_h5ad(
  sce_rna,
  filename  = "sce_rna_myeloid.h5ad",
  compression = "gzip",
  compression_opts=9,
  as_dense = "X"
)

## Save H5AD file

sceasy::convertFormat(
  sce_myeloid, 
  #main_layer = 'counts', 
  from="sce", 
  to="anndata",
  outFile='myeloid_rna.h5ad')




  
