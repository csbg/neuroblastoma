library(readxl)
library(tidyverse)
library(ComplexHeatmap)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

load("data_raw/bulk_rnaseq/decon_eda_keyDataFrames.RData")

sample_metadata <-
  read_xlsx(
    "metadata/tables_ira/BulkRNA_seq_HR_NB_DX_TU_DTC_MNC_annotations_STM_090921.xlsx",
    na = "NA"
  ) %>% 
  filter(
    OMICS_ID %in% metadata_df$omics_id,
    Material_general %in% c("DTC", "TUM"),
    !(is.na(MNA) & is.na(ATRX))
  ) %>%
  mutate(
    group = case_when(
      MNA == "YES"       ~ "M",
      ATRX == "Deletion" ~ "A",
      TRUE               ~ "S"
    ) %>% 
      factor(levels = c("M", "A", "S")),
    site = factor(Material_general, levels = c("DTC", "TUM"))
  ) %>%
  select(sample = OMICS_ID, site, group) %>% 
  arrange(site, group)

markers <- read_csv("metadata/markers_adrenal_medulla.csv", comment = "#")



# Heatmap -----------------------------------------------------------------

plot_heatmap <- function() {
  neuroblast_markers <- c(
    "MKI67", "TOP2A",                         # cycling
    "NEFM", "GAP43", "STMN2", "ISL1", "ALK",  # all
    "SYN3", "IL7"                             # late, ALK decreases
  )
  
  mat <-
    vst_counts_removedBatchEffect_df %>%
    as_tibble() %>% 
    select(ensembl_id, all_of(sample_metadata$sample)) %>% 
    left_join(
      ensemblAnnot %>% select(ensembl_id, hgnc_symbol),
      by = "ensembl_id"
    ) %>%
    select(!ensembl_id) %>%
    filter(hgnc_symbol %in% neuroblast_markers) %>%
    column_to_rownames("hgnc_symbol") %>%
    magrittr::extract(neuroblast_markers, sample_metadata$sample) %>%
    as.matrix()
  
  p <- Heatmap(
    mat,
    name = "expression",
    col = circlize::colorRamp2(
      seq(min(mat), max(mat), length.out = 9),
      RColorBrewer::brewer.pal(9, "Reds"), # Reds or RdBu
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    
    bottom_annotation = HeatmapAnnotation(
      group = sample_metadata$group,
      site = sample_metadata$site,
      col = list(
        group = GROUP_COLORS[c("M", "A", "S")],
        site = c(TUM = "black", DTC = "gray75")
      ),
      show_annotation_name = TRUE,
      show_legend = TRUE
    )
  )
  p
}

(p <- plot_heatmap())
ggsave_default("comparison/bulk_heatmap", plot = p, height = 90, width = 350)



# Heatmap of all markers --------------------------------------------------

ADRENAL_MEDULLA_COLORS <- c(
  "cycling Neuroblasts" = "#aacedc",
  "Neuroblasts" = "#244a94",
  "late Neuroblasts" = "#303d63"
)

dtc_counts <-
  vst_counts_removedBatchEffect_df %>%
  as_tibble() %>% 
  select(ensembl_id, all_of(sample_metadata$sample)) %>% 
  left_join(
    ensemblAnnot %>% select(ensembl_id, hgnc_symbol),
    by = "ensembl_id"  
  ) %>% 
  select(!ensembl_id)

all_neuroblast_markers <- 
  markers %>% 
  filter(
    gene %in% dtc_counts$hgnc_symbol,
    cell_type %>% str_detect("Neuroblast")
  ) %>% 
  distinct(gene, .keep_all = TRUE)

mat <-
  dtc_counts %>% 
  filter(hgnc_symbol %in% all_neuroblast_markers$gene) %>%
  column_to_rownames("hgnc_symbol") %>%
  magrittr::extract(all_neuroblast_markers$gene, sample_metadata$sample) %>%
  as.matrix() %>%
  t() %>% scale() %>% t()

p <- Heatmap(
  mat,
  name = "scaled\nexpression",
  col = circlize::colorRamp2(
    seq(min(mat), max(mat), length.out = 9),
    RColorBrewer::brewer.pal(9, "RdBu"), # Reds or RdBu
  ),
  
  show_row_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_split = all_neuroblast_markers$cell_type,
  show_row_dend = TRUE,
  row_title_rot = 0,
  
  bottom_annotation = HeatmapAnnotation(
    group = sample_metadata$group,
    site = sample_metadata$site,
    col = list(
      group = GROUP_COLORS[c("M", "A", "S")],
      site = c(TUM = "black", DTC = "gray75")
    ),
    show_annotation_name = TRUE,
    show_legend = TRUE
  )
)
p
ggsave_default(
  "comparison/bulk_heatmap_all_markers",
  plot = p, width = 350
)




# Single-cell expression of Jansky markers --------------------------------

# run section 'Load data' of `plot_dge_mm.R`

plotExpression(
  dge$cds[, dge$cds$cellont_abbr == "NB"],
  c("MKI67", "TOP2A", "NEFM", "GAP43", "STMN2", "ISL1", "ALK", "SYN3", "IL7"),
  x = "sample",
  colour_by = "group"
)
ggsave_default("comparison/sc_violins", height = 420)
