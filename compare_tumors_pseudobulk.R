# @DEPI rna_decontaminated.rds
# @DEPI tumor_data_dong.rds
# @DEPO tumor_data_pseudobulk.rds

library(Seurat)
library(scuttle)
library(scran)
library(batchelor)
library(muscat)
library(scico)
library(ComplexHeatmap)
library(tidyverse)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/rna_decontaminated.rds")
nb@colData <-
  readRDS("data_generated/metadata.rds") %>%
  mutate(dataset = "nb") %>% 
  column_to_rownames("cell") %>%
  as("DataFrame")

nb <- nb[, colData(nb)$cellont_abbr == "NB"]


metadata_dong <- read_csv("metadata/samples_dong.csv", comment = "#")

tumor_dong <-
  readRDS("data_generated/tumor_data_dong.rds") %>% 
  AddMetaData("dong", "dataset")

tumor_dong <- SingleCellExperiment(
  list(counts = tumor_dong$RNA@counts),
  colData = tumor_dong@meta.data %>% as("DataFrame")
)



# Merge data --------------------------------------------------------------

common_genes <- intersect(rownames(nb), rownames(tumor_dong))

out <- multiBatchNorm(
  nb[common_genes, ],
  tumor_dong[common_genes, ],
  assay.type = "counts"
)
gc()

tumor_data_merged <- correctExperiments(
  out[[1]],
  out[[2]],
  PARAM = NoCorrectParam()
)
gc()

colData(tumor_data_merged)$sample <-
  colData(tumor_data_merged)$sample %>% 
  rename_patients()

colData(tumor_data_merged)$group <-
  colData(tumor_data_merged)$sample %>% 
  str_sub(1, 1) %>%
  str_replace("N", "T")

colData(tumor_data_merged)$mycn_status <- if_else(
  colData(tumor_data_merged)$sample %in% c("T162", "T200", "T230") |
    str_starts(colData(tumor_data_merged)$sample, "M"),
  "amplified",
  "normal"
)



# Save data ---------------------------------------------------------------

list(
  pseudobulk_counts =
    tumor_data_merged %>% 
    aggregateData(by = "sample", fun = "mean") %>% 
    assay(1),
  highly_variable_genes =
    tumor_data_merged %>%
    modelGeneVar() %>% 
    getTopHVGs()
) %>%
  saveRDS("data_generated/tumor_data_pseudobulk.rds")

