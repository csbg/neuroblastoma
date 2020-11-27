# Integrate scRNAseq datasets vie one of the methods provided by Seurat
# For large datasets, run this script in batch mode in the HPC cluster
# (via sge_job_integrate.sh)

library(tidyverse)
library(fs)
library(Seurat)
library(sctransform)

# library(future)
# options(future.globals.maxSize = 4000 * 1024^2)
# plan(multicore)



# Parameters --------------------------------------------------------------

workflow <- "standard"
# workflow <- "sctrans"

data_dir <- "data_raw"


out_dir <- "data_generated/all_samples"
use_samples <- dir_ls(
  data_dir,
  recurse = TRUE,
  regexp = "filtered_feature_bc_matrix.h5"
)

# out_dir <- "data_generated/2_samples"
# use_samples <- c(
#   "R1_GNM_2020_1288_transcriptome/filtered_feature_bc_matrix.h5",
#   "MF220_NB_BM_Patient1_transcriptome/filtered_feature_bc_matrix.h5"
# ) %>%
#  {str_glue("{data_dir}/{.}")}



# Alternative 1: Standard workflow ----------------------------------------

if (workflow == "standard") {
  nb_list <- map(
    use_samples,
    function(sample) {
      sample_name <- str_match(sample, "_(.*)_trans")[,2]
      sample %>% 
        Read10X_h5() %>% 
        CreateSeuratObject() %>% 
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
  nb_list <- map(
    use_samples,
    function(sample) {
      sample_name <- str_match(sample, "_(.*)_trans")[,2]
      sample %>% 
        Read10X_h5() %>% 
        CreateSeuratObject() %>% 
        AddMetaData(sample_name, col.name = "sample") %>% 
        SCTransform()
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



# Integrated analysis -----------------------------------------------------

nb_integrated <- 
  nb_integrated %>% 
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30) %>% 
  RunTSNE(dims = 1:30)

saveRDS(nb_integrated, path_join(c(out_dir, "nb_integrated.rds")))

DimPlot(nb_integrated, group.by = "sample")
ggsave(path_join(c(out_dir, "tsne.pdf")), width = 14, height = 7)

DimPlot(nb_integrated, split.by = "sample") + NoLegend()
ggsave(path_join(c(out_dir, "tsne_split.pdf")), width = 20, height = 7)
