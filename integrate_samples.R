# Integrate scRNA-seq datasets via one of the methods provided by Seurat
# For large datasets, run this script in batch mode on the HPC cluster
# (via sge_job_integrate.sh)

library(Seurat)
library(sctransform)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

data_dir <- "data_raw"

# definition of parameters for different datasets
# each definition requires the following variables:
#   * workflow, either "sctrans" or "standard"
#   * out_dir, stores the integrated dataset (nb_integrated.rds)
#              and integration anchors (nb_anchors.rds)
#   * samples, data frame with two columns 'sample_name' and 'sample_file'

# (A) integrate all datasets
workflow <- "sctrans"
out_dir <- "data_generated/all_datasets_current"
samples <-
  tibble(
    sample_file = dir_ls(
      data_dir,
      recurse = TRUE,
      regexp = "filtered_feature_bc_matrix.h5"
    ),
    bsf_id = str_match(sample_file, "/(.*)_trans")[, 2]
  ) %>%
  left_join(
    read_csv("data_raw/metadata/sample_groups.csv", comment = "#"),
    by = "bsf_id"
  ) %>% 
    select(sample_file, sample_name = sample) %>% 
    drop_na()


# (B) test workflow on a small dataset
# workflow <- "sctrans"
# out_dir <- "data_generated/3_datasets"
# samples <- c(
#   GNM_2020_1288 = "R1_GNM_2020_1288_transcriptome/filtered_feature_bc_matrix.h5",
#   NB_2018_6056 = "R2_NB_2018_6056_transcriptome/filtered_feature_bc_matrix.h5",
#   NB_BM_Patient1 = "MF220_NB_BM_Patient1_transcriptome/filtered_feature_bc_matrix.h5"
# ) %>%
#   enframe("sample_name", "sample_file") %>% 
#   mutate(sample_file = str_glue("{data_dir}/{sample_file}"))



# Alternative 1: Standard workflow ----------------------------------------

if (workflow == "standard") {
  nb_list <- pmap(
    samples,
    function(sample_file, sample_name) {
      sample_file %>% 
        Read10X_h5() %>% 
        CreateSeuratObject(sample_name) %>% 
        AddMetaData(sample_name, col.name = "sample") %>% 
        NormalizeData() %>% 
        FindVariableFeatures()
    }
  )
  
  nb_anchors <- FindIntegrationAnchors(nb_list, dims = 1:30)
  saveRDS(nb_anchors, path_join(c(out_dir, "nb_anchors.rds")))
  
  nb_integrated <- IntegrateData(nb_anchors, dims = 1:30)
  nb_integrated <- ScaleData(nb_integrated)
  message("Integration finished")
}



# Alternative 2: SCTransform ----------------------------------------------

if (workflow == "sctrans") {
  nb_list <- pmap(
    samples,
    function(sample_file, sample_name) {
      message("SCTransforming ", sample_name)
      sample_file %>% 
        Read10X_h5() %>% 
        CreateSeuratObject(sample_name) %>% 
        AddMetaData(sample_name, col.name = "sample") %>%
        PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
        subset(subset =
                 nFeature_RNA > 200 &
                 nFeature_RNA < 5000 &
                 percent.mt < 10
        ) %>% 
        SCTransform(
          vars.to.regress = "percent.mt",
          verbose = FALSE,
          method = "qpoisson"
        )
    }
  )
  
  nb_features <- SelectIntegrationFeatures(nb_list, nfeatures = 3000)
  nb_list <- PrepSCTIntegration(nb_list, anchor.features = nb_features)
  
  nb_anchors <- FindIntegrationAnchors(
    nb_list,
    normalization.method = "SCT",
    anchor.features = nb_features
  )
  saveRDS(nb_anchors, path_join(c(out_dir, "nb_anchors.rds")))
  
  nb_integrated <- IntegrateData(
    nb_anchors,
    normalization.method = "SCT"
  )
  message("Integration finished")
}



# Save results ------------------------------------------------------------

saveRDS(nb_integrated, path_join(c(out_dir, "nb_integrated.rds")))
