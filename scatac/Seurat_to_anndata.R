pacman::p_load(Signac, Seurat, 
               zellkonverter, anndata, 
               GenomeInfoDb, EnsDb.Hsapiens.v75, 
               ggplot2, patchwork, hrbrthemes)
# Load Seurat object from disk
nblast <- 
  readRDS(
    file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized_motifs_added.rds")
DimPlot(nblast)
nblast <- RenameIdents(object = nblast, `T_cell` = "T")
nblast <- RenameIdents(object = nblast, `B_cell` = "B")
nblast <- RenameIdents(object = nblast, `Memory_B_cell` = "MB")
nblast <- RenameIdents(object = nblast, `Monocytes` = "M")
nblast <- RenameIdents(object = nblast, `NKT` = "NK")
nblast <- RenameIdents(object = nblast, `Erythroblasts` = "E")
DimPlot(nblast)

Cluster_myeloid <- 
  subset(x = nblast, 
         subset = active.ident == "Monocytes")

colsum_m <- Matrix::colSums(Cluster_myeloid@assays$peaks@counts)

rowsum_m <- Matrix::rowSums(Cluster_myeloid@assays$peaks@counts)
rowsum_zero <- rowsum_m[rowsum_m == 0]
rowsum_nozero <- rowsum_m[rowsum_m != 0]

Cluster_myeloid@assays$peaks@counts <- 
  Cluster_myeloid@assays$peaks@counts[
    rownames(Cluster_myeloid@assays$peaks@counts) %in% names(rowsum_nozero), ]

Cluster_myeloid@assays$peaks@data <- 
  Cluster_myeloid@assays$peaks@data[
    rownames(Cluster_myeloid@assays$peaks@data) %in% names(rowsum_nozero), ]

rowsum_Cluster_myeloid_new <- 
  Matrix::rowSums(Cluster_myeloid@assays$peaks@counts)
rowsum_Cluster_myeloid_new[rowsum_Cluster_myeloid_new != 0]

Cluster_myeloid <- 
  as.loom(Cluster_myeloid, 
                     filename = "~/Cluster_myeloid_atac.loom", 
          verbose = FALSE)
pbmc.loom

## Save H5AD file

Cluster_myeloid[['RNA']] <- NULL

sceasy::convertFormat(
  Cluster_myeloid, 
  main_layer = 'counts', 
  from="seurat", 
  to="anndata",
  outFile='sce_atac_myeloid.h5ad')


Assay = 'peaks'
# Convert Seurat object to SingleCellExperiment object
sce_peaks_m <- 
  as.SingleCellExperiment(
    Cluster_myeloid, 
    assay = 'peaks')

# Convert SingleCellExperiment object to AnnData object
sce_atac_m <- 
  SCE2AnnData(
    sce_peaks_m, 
    X_name = 'counts', 
    skip_assays = FALSE)



# Save AnnData object
write_h5ad(
  sce_atac_m,
  filename  = "sce_atac_m.h5ad",
  compression = "gzip",
  compression_opts=9,
  as_dense = "X"
)




