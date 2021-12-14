library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(scico)
source("common_functions.R")
source("styling.R")



# # Load data ---------------------------------------------------------------
# 
# # run compare_tumors_pbexp.R, section "load data
# 
# colData(nb)$sample <- rename_patients(colData(nb)$sample)
# 
# pb_nb <-
#   nb %>% 
#   aggregateData(by = "sample", fun = "mean") %>% 
#   assay(1) %>% 
#   as_tibble(rownames = "gene")
# 
# metadata_nb <- tibble(
#   sample = colnames(pb_nb)[-1],
#   site = "DTC",
#   mycn_status = if_else(str_starts(sample, "M"), "amplified", "normal"),
#   method = "sc"
# )
# 
# 
# pb_dong <-
#   tumor_dong %>% 
#   aggregateData(by = "sample", fun = "mean") %>% 
#   assay(1) %>% 
#   as_tibble(rownames = "gene")
# 
# metadata_dong <- tibble(
#   sample = colnames(pb_dong)[-1],
#   site = "primary",
#   mycn_status = if_else(sample %in% c("T162", "T200", "T230"), "amplified", "normal"),
#   method = "sc"
# )
# 
# 
# load("data_raw/bulk_rnaseq/decon_eda_keyDataFrames.RData")
# 
# metadata_bulk <-
#   read_xlsx(
#     "metadata/original_tables/BulkRNA_seq_HR_NB_DX_TU_DTC_MNC_annotations_STM_090921.xlsx",
#     na = "NA"
#   ) %>% 
#   filter(
#     OMICS_ID %in% metadata_df$omics_id,
#     Material_general %in% c("DTC", "TUM"),
#     !(is.na(MNA) & is.na(ATRX))
#   ) %>%
#   mutate(
#     mycn_status = case_when(
#       MNA == "YES" ~ "amplified",
#       TRUE         ~ "normal"
#     ),
#     site = Material_general %>% fct_recode(primary = "TUM"),
#     method = "bulk"
#   ) %>%
#   select(sample = OMICS_ID, site, mycn_status, method) %>% 
#   arrange(site, mycn_status)
# 
# bulk_data <- 
#   vst_counts_removedBatchEffect_df %>%
#   as_tibble() %>% 
#   select(ensembl_id, all_of(metadata_bulk$sample)) %>%
#   left_join(
#     ensemblAnnot %>% select(ensembl_id, hgnc_symbol),
#     by = "ensembl_id"
#   ) %>%
#   select(!ensembl_id) %>%
#   relocate(gene = hgnc_symbol)
# 
# 
# 
# # Remove batch effect -----------------------------------------------------
# 
# count_mat <-
#   pb_nb %>% 
#   inner_join(pb_dong, by = "gene") %>% 
#   inner_join(bulk_data, by = "gene") %>% 
#   mutate(gene = make.unique(gene)) %>% 
#   column_to_rownames("gene") %>% 
#   as.matrix()
# 
# keep <- filterByExpr(count_mat)
# 
# col_metadata <-
#   tibble(sample = colnames(count_mat)) %>% 
#   left_join(
#     bind_rows(metadata_nb, metadata_dong, metadata_bulk),
#     by = "sample"
#   ) %>% 
#   mutate(across(everything(), as_factor))
# 
# count_mat_nobatch <- ComBat_seq(count_mat[keep, ], batch = col_metadata$method)
# 
# count_mat_nobatch <- removeBatchEffect(
#   count_mat[keep, ],
#   batch = col_metadata$method,
#   design = model.matrix(~site + mycn_status, data = col_metadata)
# )


# Bulk data only ----------------------------------------------------------

load("data_raw/bulk_rnaseq/decon_eda_keyDataFrames.RData")

col_metadata <-
  read_xlsx(
    "metadata/original_tables/BulkRNA_seq_HR_NB_DX_TU_DTC_MNC_annotations_STM_090921.xlsx",
    na = "NA"
  ) %>% 
  filter(
    OMICS_ID %in% metadata_df$omics_id,
    Material_general %in% c("DTC", "TUM"),
    !(is.na(MNA) & is.na(ATRX))
  ) %>%
  mutate(
    group =
      case_when(
        MNA == "YES"       ~ "M",
        ATRX == "Deletion" ~ "A",
        TRUE               ~ "S"
      ) %>% 
      factor(levels = c("M", "A", "S")),
    mycn_status = case_when(
      MNA == "YES" ~ "amplified",
      TRUE         ~ "normal"
    ),
    site = Material_general %>% fct_recode(primary = "TUM")
  ) %>%
  select(sample = OMICS_ID, group, mycn_status, site)

count_mat_nobatch <-
  vst_counts_removedBatchEffect_df %>%
  as_tibble() %>%
  select(ensembl_id, all_of(metadata_bulk$sample)) %>%
  left_join(
    ensemblAnnot %>% select(ensembl_id, hgnc_symbol),
    by = "ensembl_id"
  ) %>%
  mutate(gene = make.unique(hgnc_symbol)) %>%
  select(!c(ensembl_id, hgnc_symbol)) %>%
  column_to_rownames("gene") %>% 
  as.matrix()

keep <- edgeR::filterByExpr(count_mat_nobatch, group = col_metadata$group)
count_mat_nobatch <- count_mat_nobatch[keep, ]



# Analyse data ------------------------------------------------------------

corr_mat <-
  count_mat_nobatch %>%
  cor(use = "pairwise.complete.obs")

distance <- as.dist(1 - corr_mat)

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

  width = unit(3.7, "mm") * nrow(corr_mat),
  height = unit(3.7, "mm") * nrow(corr_mat),
  # width = unit(80, "mm"),
  # height = unit(80, "mm"),

  show_column_dend = FALSE,
  show_column_names = FALSE,

  left_annotation = rowAnnotation(
    group = col_metadata$group,
    site = col_metadata$site,
    mycn = col_metadata$mycn_status,
    col = list(
      group = c(GROUP_COLORS, "T" = "#433447"),
      site = c(primary = "black", DTC = "gray75"),
      mycn = c(normal = "gray90", amplified = "#d35f5f")
    ),
    show_annotation_name = FALSE,
    show_legend = TRUE
  ),
)
p
ggsave_default("comparison/correlation_pseudobulk_bulk", plot = p, height = 420, width = 420)







