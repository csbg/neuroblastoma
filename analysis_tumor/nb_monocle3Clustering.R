source("nb_Packages.R") 

################################################################################
# monocle 3: Normalize, pre-process,align & reduce dimensions of the selected CDS
# (UMAP -> core=1 and umap.fast_sgd = FALSE for accurate and repeatable results
################################################################################

#import data
cds <- readRDS(
  paste0(path,"nb_CDS_raw.rds")
  ) #change link

# Pre-Process (normalization + PCA)
cds_unaligned <- preprocess_cds(
  cds,
  method = "PCA",
  num_dim = 50,
  norm_method = "log"
)

# Alignment (remove patient & cell-cycle phase effect)
cds <- align_cds(
  cds_unaligned,
  preprocess_method = "PCA",
  alignment_group = "phase",
  residual_model_formula_str = 
    "~sample + norm.scores_S + norm.scores_G1 + norm.scores_G2M",
  alignment_k = 20,
  verbose = TRUE
)

#Dimensional reduction (preprocessed & corrected data)
cds <- reduce_dimension(
  cds, 
  reduction_method = "UMAP",
  preprocess_method = "Aligned"
)  

#Dimensional reduction (unaligned data)
cds_unaligned <- reduce_dimension(
  cds_unaligned, 
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

#cluster unaligned data
cds_unaligned <- cluster_cells(
  cds_unaligned, 
  resolution = 1e-4
)

#cluster aligned data to 4 clusters
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
    select(1:6,17:40)

# export data
saveRDS(cds,
        file = paste0(path,"/nb_CDS.rds") #change link
        )

saveRDS(metdata,
        file = paste0(path,"/nb_Metadata_tumor.rds") #change link
        )

saveRDS(cds_unaligned,
        file = paste0(path,"/nb_CDS_unaligned.rds") #change link
        )

