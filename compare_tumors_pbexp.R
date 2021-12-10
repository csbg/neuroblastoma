# Create pseudobulk correlation heatmap of own and Dong datasets
#
# @DEPI rna_decontaminated.rds
# @DEPI tumor_data_dong.rds

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

# metadata_jansky <-
#   read_csv("data_wip/metadata_samples_jansky.csv", comment = "#")
# 
# used_samples_jansky <- 
#   metadata_jansky %>% 
#   filter(
#     clinical_subtype %in% c("MYCN", "TERT", "ALT"),
#     !mesenchymal_features
#   ) %>% 
#   pull(sample)
# 
# tumor_jansky <-
#   readRDS("data_wip/tumor_data_jansky.rds") %>% 
#   subset(subset = anno_new == "Tumor cells") %>% 
#   AddMetaData("jansky", "dataset")
# 
# tumor_jansky <-
#   SingleCellExperiment(
#     list(counts = tumor_jansky$RNA@counts),
#     colData =
#       tumor_jansky@meta.data %>% 
#       select(sample = patientID, dataset) %>% 
#       as("DataFrame")
#   ) %>% 
#   magrittr::extract(, colData(.)$sample %in% used_samples_jansky)



# Merge data --------------------------------------------------------------

common_genes <-
  rownames(nb) %>% 
  intersect(rownames(tumor_dong))
  # intersect(rownames(tumor_jansky))

out <- multiBatchNorm(
  nb[common_genes, ],
  tumor_dong[common_genes, ],
  # tumor_jansky[common_genes, ],
  assay.type = "counts"
)
gc()

tumor_data_merged <- correctExperiments(
  out[[1]],
  out[[2]],
  # out[[3]],
  PARAM = RescaleParam()
)
gc()

colData(tumor_data_merged)$sample <-
  colData(tumor_data_merged)$sample %>% 
  rename_patients()

colData(tumor_data_merged)$group <-
  colData(tumor_data_merged)$sample %>% 
  str_sub(1, 1) %>%
  str_replace("N", "T")

# use if Jansky data is used
# colData(tumor_data_merged)$mycn_status <- if_else(
#   colData(tumor_data_merged)$sample %in% c("T162", "T200", "T230",
#                 "NB01", "NB02", "NB08", "NB11", "NB14") |
#     str_starts(colData(tumor_data_merged)$sample, "M"),
#   "amplified",
#   "normal"
# )

colData(tumor_data_merged)$mycn_status <- if_else(
  colData(tumor_data_merged)$sample %in% c("T162", "T200", "T230") |
    str_starts(colData(tumor_data_merged)$sample, "M"),
  "amplified",
  "normal"
)

pb_tumor <-
  tumor_data_merged %>% 
  aggregateData(by = "sample", fun = "mean")



# Plot heatmap ------------------------------------------------------------

hvgs <-
  tumor_data_merged %>%
  modelGeneVar() %>% 
  getTopHVGs()

corr_mat <-
  assay(pb_tumor, 1) %>%
  magrittr::extract(hvgs, ) %>%
  cor(use = "pairwise.complete.obs")

distance <- as.dist(1 - corr_mat)

col_metadata <- tibble(
  sample = colnames(corr_mat),
  group = str_sub(sample, 1, 1) %>% str_replace("N", "T"),
  mycn_status = if_else(
    sample %in% c("T162", "T200", "T230",
                  "NB01", "NB02", "NB08", "NB11", "NB14") |
    str_starts(sample, "M"),
    "amplified",
    "normal"
  )
)

p <- Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), max(corr_mat), length.out = 9),
    scico(9, palette = "davos", direction = -1),
  ),
  name = "correlation of\npseudobulk\nexpression",
  heatmap_legend_param = list(
    at = round(c(min(corr_mat), max(corr_mat)), 2)
  ),
  
  clustering_distance_rows = distance,
  clustering_distance_columns = distance,
  
  # width = unit(3.7, "mm") * nrow(corr_mat),
  # height = unit(3.7, "mm") * nrow(corr_mat),
  width = unit(80, "mm"),
  height = unit(80, "mm"),
  
  show_column_dend = FALSE,
  
  left_annotation = rowAnnotation(
    group = col_metadata$group,
    mycn = col_metadata$mycn_status,
    col = list(
      group = c(GROUP_COLORS, "T" = "#433447"),
      mycn = c("normal" = "gray90", "amplified" = "#d35f5f")
    ),
    show_annotation_name = FALSE,
    show_legend = TRUE,
    annotation_legend_param = list(
      group = list(
        title = "group"
      ),
      mycn = list(
        title = "MYCN status"
      )
    )
  ),
)
p
ggsave_default("comparison/correlation_pseudobulk", plot = p)
