# Integrate scRNA-seq datasets Seurat (sctransform and reciprocal PCA).
#
# Generates rna_integrated_seurat.rds, containing the intrgrated dataset.
#
# @DEPI rna_qcpassed.rds
# @DEPO rna_integrated_seurat.rds

library(Seurat)
library(sctransform)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# the QC-filtered datasets
filtered_datasets <- "data_generated/rna_qcpassed.rds"

# folder where results are saved
out_dir <- "data_generated"



# Integrate data ----------------------------------------------------------

nb_list <-
  readRDS(filtered_datasets) %>% 
  map(
    SCTransform,
    vars.to.regress = "percent_mito",
    method = "qpoisson",
    verbose = FALSE
  )

nb_features <- SelectIntegrationFeatures(nb_list)

nb_list <-
  nb_list %>%
  PrepSCTIntegration(anchor.features = nb_features) %>% 
  map(RunPCA, features = nb_features)

nb_anchors <- FindIntegrationAnchors(
  nb_list,
  normalization.method = "SCT",
  anchor.features = nb_features,
  reduction = "rpca"
)

nb_integrated <- IntegrateData(
  nb_anchors,
  normalization.method = "SCT"
)
message("Integration finished")



# Save results ------------------------------------------------------------

saveRDS(nb_integrated, path_join(c(out_dir, "rna_integrated_seurat.rds")))
