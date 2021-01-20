# Merge scRNA-seq datasets without integration
# Run this script in batch mode on the HPC cluster
# (via sge_job_merge.sh)

library(Seurat)
library(sctransform)
library(tidyverse)
library(fs)



# Parameters --------------------------------------------------------------

data_dir <- "data_raw"
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



# Merge datasets ----------------------------------------------------------

nb_list <- pmap(
  samples,
  function(sample_file, sample_name) {
    message("Loading ", sample_name)
    sample_file %>% 
      Read10X_h5() %>% 
      CreateSeuratObject(sample_name) %>% 
      AddMetaData(sample_name, col.name = "sample") %>%
      PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
      subset(subset =
               nFeature_RNA > 200 &
               nFeature_RNA < 5000 &
               percent.mt < 10
      )
  }
)

nb <-
  merge(nb_list[[1]], nb_list[-1]) %>%  
  SCTransform(vars.to.regress = "percent.mt") %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30) %>% 
  RunTSNE(dims = 1:30)


# Save results ------------------------------------------------------------

saveRDS(nb, path_join(c(out_dir, "nb_merged.rds")))
