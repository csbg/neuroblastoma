library(scico)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)
library(fs)
source("common_functions.R")

ht_opt(message = FALSE, show_parent_dend_line = FALSE)


ifob <- readRDS("~/Desktop/run.final.infercnv_obj")
nb_data <- readRDS("data_generated/metadata.rds")


plot_cnv <- function(infercnv_obj, metadata, cells = c("tumor", "ref"),
                     cells_per_sample = 50L, filename = NULL) {
  cells <- match.arg(cells)
  
  # data frame with tumor or reference cell indices
  if (cells == "tumor") {
    # order of tumor samples
    tumor_samples <-
      metadata %>% 
      filter(group != "I") %>% 
      pull(sample) %>% 
      fct_drop() %>% 
      levels()
    
    cell_metadata <- 
      infercnv_obj@tumor_subclusters$subclusters %>%
      enframe() %>% 
      filter(str_starts(name, "malignant")) %>% 
      extract(
        name,
        into = c("sample", "group"),
        regex = "malignant_(.*)_([IV]+)"
      ) %>% 
      mutate(sample = as_factor(sample) %>% fct_relevel(tumor_samples)) %>% 
      arrange(sample)
  } else {
    cell_metadata <- 
      infercnv_obj@tumor_subclusters$subclusters %>%
      enframe(name = "sample") %>% 
      filter(!str_starts(sample, "malignant"))
  }
  
  cell_metadata <- 
    cell_metadata %>% 
    unnest_longer(value) %>%
    select(!value_id) %>%
    unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>% 
    left_join(metadata %>% select(cell, patient = sample), by = "cell")
  
  # optional: subset n cells per sample
  if (!is.null(cells_per_sample)) {
    cell_metadata <- 
      cell_metadata %>% 
      group_by(sample) %>% 
      slice_sample(n = cells_per_sample) %>% 
      ungroup()
  }
  
  # matrix with residual expression, rows ordered
  cnv_mat <-
    infercnv_obj@expr.data %>% 
    magrittr::extract(, cell_metadata$cell_index) %>% 
    t()
  
  # annotation (group for tumor cells, sample for reference cells)
  if (cells == "tumor") {
    left_annotation <- rowAnnotation(
      group = cell_metadata$group,
      col = list(
        # colours.cafe 502
        group = c("II" = "#154d42", "III" = "#dcbc65", "IV" = "#b7336c")
      )
    )
  } else {
    patient_ids <- unique(cell_metadata$patient)
    left_annotation <- rowAnnotation(
      sample = cell_metadata$patient,
      col = list(
        sample = rainbow(length(patient_ids)) %>% set_names(patient_ids)
      )
    )
  }
  
  # draw heatmap
  p <- 
    Heatmap(
      cnv_mat,
      name = "residual\nexpression",
      col = circlize::colorRamp2(
        breaks = seq(0.85, 1.15, length.out = 7),
        colors = rev(brewer.pal(7, "RdBu"))
      ),
      show_row_names = FALSE,
      show_column_names = FALSE,
      
      row_split = cell_metadata$sample,
      row_gap = unit(0, "mm"),
      row_title_rot = 0,
      
      cluster_columns = FALSE,
      column_split =
        infercnv_obj@gene_order$chr %>%
        str_sub(4) %>%
        as_factor(),
      cluster_column_slices = FALSE,
      column_gap = unit(0, "mm"),
      
      border = TRUE,
      
      left_annotation = left_annotation
    ) %>%
    draw(
      row_title = if (cells == "tumor") "sample",
      column_title = "chromosome"
    ) 
  
  ggsave_default(filename, plot = p)
}

plot_cnv(ifob, nb_data, filename = "cnv/cnv")
plot_cnv(ifob, nb_data, cells = "ref", filename = "cnv/ref")





