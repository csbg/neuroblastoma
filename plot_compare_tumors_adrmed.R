# @DEPI adrmed_class_*.csv

library(ComplexHeatmap)
library(scico)
library(tidyverse)
library(fs)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

## Dong ----

metadata_dong <-
  read_csv("metadata/samples_dong.csv", comment = "#") %>% 
  filter(high_risk) %>% 
  transmute(sample = tumor_id, group = if_else(mycn_amplified, "T-M", "T-S"))

singler_results_dong <- 
  dir_ls("data_generated/adrmed", glob = "*class_dong*") %>% 
  map_dfr(read_csv, .id = "file") %>% 
  extract(file, into = "sample", regex = "dong_(.+)\\.") %>% 
  left_join(metadata_dong, by = "sample") %>% 
  relocate(group, .after = sample) %>% 
  mutate(dataset = "dong", .before = 1)



## Own data ----

metadata_nb <-
  readRDS("data_generated/metadata.rds") %>% 
  transmute(cell, sample = rename_patients(sample), group = rename_groups(group))

singler_results_nb <- 
  read_csv("data_generated/adrmed/adrmed_class_nb.csv") %>% 
  left_join(metadata_nb, by = "cell") %>% 
  relocate(sample, group) %>% 
  mutate(dataset = "nb", .before = 1)


## Combine data ----

singler_results <-
  bind_rows(
    singler_results_dong,
    singler_results_nb
  ) %>% 
  mutate(
    mycn_status = if_else(group %in% c("M", "T-M"), "amplified", "normal"),
    tumor = if_else(group %in% c("A", "M", "S"), "DTC", "primary"),
    pruned_labels =
      pruned_labels %>% 
      factor(levels = names(ADRMED_CELLS_COLORS)) %>% 
      fct_rev() %>% 
      fct_explicit_na()
  )



# Analyze -----------------------------------------------------------------

## Cell type profile ----

plot_profile <- function(group_col = group) {
  sample_order <- 
    singler_results %>% 
    distinct(group = {{group_col}}, sample) %>% 
    arrange(group, sample) %>% 
    pull(sample)
  
  singler_results %>% 
    mutate(sample = factor(sample, levels = sample_order)) %>% 
    group_by(group = {{group_col}}, sample) %>%
    count(cell_type = pruned_labels) %>% 
    mutate(n_rel = n / sum(n) * 100) %>% 
    ggplot(aes(sample, n_rel)) +
    geom_col(aes(fill = group)) +
    geom_hline(yintercept = 0) +
    ylab("relative abundance (%)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_grid(
      vars(cell_type),
      vars(group),
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank()
    )  
}

plot_profile()
ggsave_default("comparison/cell_types_group", height = 420)


plot_agg_profile <- function(data, group_col = group) {
  data %>%
    filter(dataset != "jansky") %>% 
    group_by(group = {{group_col}}) %>% 
    count(cell_type = pruned_labels) %>% 
    mutate(n_rel = n / sum(n) * 100) %>% 
    ungroup() %>% 
    ggplot(aes(cell_type, n_rel)) +
    geom_col(aes(fill = cell_type), show.legend = FALSE) +
    scale_x_discrete(NULL, drop = FALSE) +
    ylab("Percentage of cells") +
    scale_fill_manual(values = c(ADRMED_CELLS_COLORS, "(Missing)" = "black")) +
    coord_flip() +
    facet_wrap(vars(group), nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}

singler_results %>% 
  mutate(
    group =
      fct_cross(tumor, mycn_status) %>%
      fct_relevel("DTC:amplified", "DTC:normal")
  ) %>% 
  plot_agg_profile()
ggsave_default("comparison/profile_group", height = 100)



## Patient heatmap ----

plot_patient_heatmap <- function() {
  mat <- 
    singler_results %>% 
    filter(dataset != "jansky") %>%
    group_by(sample) %>%
    count(cell_type = pruned_labels) %>% 
    mutate(n = n / sum(n) * 100) %>% 
    ungroup() %>% 
    pivot_wider(names_from = sample, values_from = n) %>% 
    column_to_rownames("cell_type") %>%
    as.matrix() %>% 
    replace_na(0)
  
  col_metadata <- 
    tibble(sample = colnames(mat)) %>% 
    left_join(
      singler_results %>% distinct(sample, mycn_status, tumor),
      by = "sample"
    )
  
  Heatmap(
    mat,
    col = RColorBrewer::brewer.pal(9, "YlOrBr"),
    name = "relative\nabundance",
    heatmap_legend_param = list(
      at = round(c(min(mat), max(mat)), 2)
    ),
    
    column_split =
      fct_cross(col_metadata$tumor, col_metadata$mycn_status) %>% 
      fct_relevel("DTC:amplified", "DTC:normal"),
    cluster_column_slices = FALSE,
    column_title = NULL,
    
    top_annotation = HeatmapAnnotation(
      tumor = col_metadata$tumor,
      mycn = col_metadata$mycn_status,
      col = list(
        mycn = c("normal" = "gray90", "amplified" = "#d35f5f"),
        tumor = c(DTC = "black", primary = "grey80")
      ),
      annotation_legend_param = list(
        mycn = list(title = "MYCN status"),
        tumor = list(title = "tumor")
      )
    ),
  )
}
(p <- plot_patient_heatmap())

ggsave_default("comparison/heatmap_adrmed_celltypes",
               plot = p, height = 100, width = 230)
