library(Seurat)
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
  column_to_rownames("cell") %>%
  as("DataFrame")

nb <- nb[, colData(nb)$cellont_abbr == "NB"]

tumor_dong <- readRDS("data_wip/tumor_data_dong.rds")

tumor_jansky <-
  readRDS("data_wip/tumor_data_jansky.rds") %>% 
  subset(subset = anno_new == "Tumor cells")

tumor_jansky <- CreateSeuratObject(
  tumor_jansky$RNA@counts,
  meta.data =
    tumor_jansky@meta.data["patientID"] %>%
    magrittr::set_colnames("sample")
)



# Merge data --------------------------------------------------------------

tumor_dong$nb <- CreateSeuratObject(counts(nb))
tumor_dong$nb@meta.data$sample <- rename_patients(colData(nb)$sample)

tumor_data_merged <-
  merge(tumor_jansky, tumor_dong) %>% 
  as.SingleCellExperiment() %>% 
  scuttle::logNormCounts()



# Plot heatmap ------------------------------------------------------------

pb_tumor <-
  tumor_data_merged %>% 
  # tumor_data_merged[common_genes, ] %>% 
  aggregateData(
    assay = "logcounts",
    fun = "mean",
    by = "sample"
  )

hvgs <-
  tumor_data_merged %>%
  # tumor_data_merged[common_genes, ] %>%
  scran::modelGeneVar() %>% 
  scran::getTopHVGs()

corr_mat <-
  assay(pb_tumor, 1) %>%
  magrittr::extract(hvgs, ) %>%
  cor(use = "pairwise.complete.obs")

distance <- as.dist(1 - corr_mat)

# group_names <-
#   colnames(corr_mat) %>%
#   map_chr(str_sub, 1, 1)
# 
# mycn_status <- 
#   files_dong %>%
#   select(tumor_id, mycn_amplified) %>%
#   deframe() %>% 
#   append(
#     PATIENT_ORDER %>% 
#       str_detect("M") %>% 
#       set_names(PATIENT_ORDER)
#   ) %>% 
#   magrittr::extract(colnames(corr_mat)) %>% 
#   if_else("amplified", "normal")

p <- Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), 1, length.out = 9),
    scico(9, palette = "davos", direction = -1),
  ),
  name = "correlation of\npseudobulk\nexpression",
  # heatmap_legend_param = list(
  #   at = c(round(min(corr_mat), 2), 0.9, 1)
  # ),
  
  clustering_distance_rows = distance,
  clustering_distance_columns = distance,
  
  width = unit(150, "mm"),
  height = unit(150, "mm"),
  
  show_column_dend = FALSE,
  
  # left_annotation = rowAnnotation(
  #   group = group_names,
  #   mycn = mycn_status,
  #   col = list(
  #     group = c(GROUP_COLORS, "T" = "#433447"),
  #     mycn = c("normal" = "gray90", "amplified" = "#d35f5f")
  #   ),
  #   show_annotation_name = FALSE,
  #   show_legend = TRUE,
  #   annotation_legend_param = list(
  #     group = list(
  #       title = "group"
  #     ),
  #     mycn = list(
  #       title = "MYCN status"
  #     )
  #   )
  # ),
)
p
ggsave_default("comparison/correlation_jansky", plot = p)




# tumor_data_merged %>% counts() %>% dim()
# 
# 
# jansky_genes <- tibble(
#   gene = rownames(tumor_jansky$RNA@counts),
#   jansky = TRUE
# )
# nb_genes <- tibble(
#   gene = rownames(nb),
#   nb = TRUE
# )
# dong_genes <- tibble(
#   gene = rownames(tumor_dong[[1]]$RNA@counts),
#   dong = TRUE
# )
# 
# 
# intersect(nb_genes$gene, dong_genes$gene) %>% length()
# 
# common_genes <- 
#   jansky_genes %>% 
#   full_join(nb_genes, by = "gene") %>% 
#   full_join(dong_genes, by = "gene") %>% 
#   mutate(gene_in_all = jansky & nb) %>% 
#   # arrange(gene) %>% 
#   # dplyr::count(gene_in_all) %>% 
#   filter(gene_in_all) %>% 
#   pull(gene)
#   View()
# 
#   
# tumor_data_merged[common_genes, ]
