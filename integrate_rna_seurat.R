# Integrate scRNA-seq datasets via SCTransform and Seurat.
#
# Generates rna_integrated_seurat.rds, containing the intrgrated dataset.

library(Seurat)
library(sctransform)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

# the merged dataset
merged_dataset <- "data_generated/rna_merged.rds"

# folder where results are saved
out_dir <- "data_generated"



# Integrate data ----------------------------------------------------------

nb_merged <- readRDS(merged_dataset)

nb_list <-
  nb_merged %>% 
  SplitObject("sample") %>% 
  map(
    SCTransform,
    vars.to.regress = "percent.mt",
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
