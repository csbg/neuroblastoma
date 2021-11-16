library(Seurat)
library(batchelor)
library(muscat)
library(scico)
library(ComplexHeatmap)
library(tidyverse)
source("common_functions.R")
source("styling.R")



# Original code -----------------------------------------------------------

tumor_data_merged <-
  readRDS("data_wip/tumor_data.rds") %>% 
  {merge(.[[1]], .[-1])} %>% 
  as.SingleCellExperiment() %>% 
  logNormCounts()

colData(tumor_data_merged)$sample <- 
  colData(tumor_data_merged)$sample %>%
  rename_patients()

pb_tumor <-
  tumor_data_merged %>% 
  aggregateData(
    assay = "logcounts",
    fun = "mean",
    by = "sample"
  )

hvgs <-
  tumor_data_merged %>%
  scran::modelGeneVar() %>% 
  scran::getTopHVGs()

corr_mat <-
  assay(pb_tumor, 1) %>%
  magrittr::extract(hvgs, ) %>%
  cor(use = "pairwise.complete.obs")

distance <- as.dist(1 - corr_mat)

group_names <-
  colnames(corr_mat) %>%
  map_chr(str_sub, 1, 1)

mycn_status <- 
  files_dong %>%
  select(tumor_id, mycn_amplified) %>%
  deframe() %>% 
  append(
    PATIENT_ORDER %>% 
      str_detect("M") %>% 
      set_names(PATIENT_ORDER)
  ) %>% 
  magrittr::extract(colnames(corr_mat)) %>% 
  if_else("amplified", "normal")

p <- Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), 1, length.out = 9),
    scico(9, palette = "davos", direction = -1),
  ),
  name = "correlation of\npseudobulk\nexpression",
  heatmap_legend_param = list(
    at = c(round(min(corr_mat), 2), 0.9, 1)
  ),
  
  clustering_distance_rows = distance,
  clustering_distance_columns = distance,
  
  width = unit(80, "mm"),
  height = unit(80, "mm"),
  
  show_column_dend = FALSE,
  
  left_annotation = rowAnnotation(
    group = group_names,
    mycn = mycn_status,
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
ggsave_default("comparison/correlation", plot = p)






# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/rna_decontaminated.rds")
nb@colData <-
  readRDS("data_generated/metadata.rds") %>%
  mutate(dataset = "nb") %>% 
  column_to_rownames("cell") %>%
  as("DataFrame")

nb <- nb[, colData(nb)$cellont_abbr == "NB"]

tumor_dong <-
  readRDS("data_wip/tumor_data_dong.rds") %>% 
  # {merge(.[[1]], .[2:3])} %>%  # a subset of tumors for testing
  {merge(.[[1]], .[-1])} %>%  # all tumors
  AddMetaData("dong", "dataset")

tumor_dong <- SingleCellExperiment(
  list(counts = tumor_dong$RNA@counts),
  colData = tumor_dong@meta.data %>% as("DataFrame")
)
  
tumor_jansky <-
  readRDS("data_wip/tumor_data_jansky.rds") %>% 
  subset(subset = anno_new == "Tumor cells") %>% 
  AddMetaData("jansky", "dataset")

tumor_jansky <- SingleCellExperiment(
  list(counts = tumor_jansky$RNA@counts),
  colData =
    tumor_jansky@meta.data %>% 
    select(sample = patientID, dataset) %>% 
    as("DataFrame")
)

tumor_jansky <- tumor_jansky[
  ,
  colData(tumor_jansky)$sample %in% c("NB01", "NB03", "NB05",
                                      "NB08", "NB10", "NB11")
]


# Merge data --------------------------------------------------------------

common_genes <-
  rownames(nb) %>% 
  intersect(rownames(tumor_dong)) %>%
  intersect(rownames(tumor_jansky))

out <- multiBatchNorm(
  nb[common_genes, ],
  tumor_dong[common_genes, ],
  tumor_jansky[common_genes, ]
)

tumor_data_merged <- correctExperiments(
  out[[1]], out[[2]], out[[3]],
  PARAM = RescaleParam()
)

colData(tumor_data_merged)$sample =
  colData(tumor_data_merged)$sample %>% 
  rename_patients()



# Plot heatmap ------------------------------------------------------------

pb_tumor <-
  tumor_data_merged %>% 
  aggregateData(
    assay = "corrected",
    fun = "mean",
    by = "sample"
  )

hvgs <-
  tumor_data_merged %>%
  scran::modelGeneVar() %>% 
  scran::getTopHVGs()

corr_mat <-
  assay(pb_tumor, 1) %>%
  magrittr::extract(hvgs, ) %>%
  cor(use = "pairwise.complete.obs")

distance <- as.dist(1 - corr_mat)

col_metadata <-
  tibble(
    sample = colnames(corr_mat),
    group = case_when(
      sample %>% str_starts("M") ~ "M",
      sample %in% c("T162", "T200", "T230", "NB01", "NB08", "NB11") ~ "M",
      sample %>% str_starts("[TAS]") ~ "A+S",
      sample %in% c("NB03", "NB05", "NB10", "NB12") ~ "A+S",
      TRUE ~ "other"
    ),
    tumor = if_else(str_starts(sample, "[AMS]"), "DTC", "primary")
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
  
  width = unit(150, "mm"),
  height = unit(150, "mm"),
  
  show_column_dend = FALSE,
  
  left_annotation = rowAnnotation(
    group = col_metadata$group,
    tumor = col_metadata$tumor,
    col = list(
      group = c(
        GROUP_COLORS,
        "A+S" = "#dd9a59",
        "T-A" = "#8a504e",
        "T-M" = "#768688",
        "T-S" = "#967F52",
        "T-A+S" = "#906850",
        other = "grey80",
        low_risk = "grey80",
        mesenchymal = "#c61f3c"
      ),
      tumor = c(DTC = "black", primary = "grey80")
    )
  ),
)
p
ggsave_default("comparison/correlation_jansky", plot = p)
