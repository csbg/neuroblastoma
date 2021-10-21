source("analysis_tumor/packages.R") 

# import data
cds <- readRDS("analysis_tumor/data_generated/cds_raw.rds")

# pre-process (normalization + PCA)
cds_unaligned <- preprocess_cds(
  cds,
  method = "PCA",
  num_dim = 50,
  norm_method = "log"
)

# align (remove patient & cell-cycle phase effect)
cds <- align_cds(
  cds_unaligned,
  preprocess_method = "PCA",
  alignment_group = "phase",
  residual_model_formula_str = 
    "~sample + norm.scores_S + norm.scores_G1 + norm.scores_G2M",
  alignment_k = 20,
  verbose = TRUE
)

# dimensional reduction (aligned data)
cds <- reduce_dimension(
  cds, 
  reduction_method = "UMAP",
  preprocess_method = "Aligned"
)  

# dimensional reduction (unaligned data)
cds_unaligned <- reduce_dimension(
  cds_unaligned, 
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

# cluster unaligned data
cds_unaligned <- cluster_cells(
  cds_unaligned, 
  resolution = 1e-4
)

# cluster aligned data to 4 clusters
set.seed(0)
cds <- cluster_cells(
  cds, 
  reduction_method = "UMAP",
  resolution = 4e-4
)

plot_cells(cds, color_cells_by = "cluster", cell_size = 0.8)

# assign aligned clusters as metadata to both CDS
cds@colData$cluster <- cds@clusters$UMAP$clusters
cds_unaligned@colData$aligned_cluster <- cds@clusters$UMAP$clusters

metadata <- 
  cds@colData %>% 
  as_tibble() %>%
  select(1:6, 17:40)

# export data
saveRDS(cds, "analysis_tumor/data_generated/cds.rds")
saveRDS(metadata, "analysis_tumor/data_generated/metadata_tumor.rds")
saveRDS(cds_unaligned, "analysis_tumor/data_generated/cds_unaligned.rds")
