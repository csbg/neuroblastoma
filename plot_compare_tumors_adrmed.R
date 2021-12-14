# Evaluate results of adrenal medullar classification.
#
# @DEPI adrmed_class_*.csv

library(ComplexHeatmap)
library(scico)
library(tidyverse)
library(fs)
source("common_functions.R")
source("styling.R")


ADRMED_CELLS_COLORS <- c(
  "late SCPs" = "#921813",
  "SCPs" = "#be202e",
  "cycling SCPs" = "#be6867",
  "Bridge" = "#ef845f",
  "connecting Chromaffin cells" = "#33ac75",
  "Chromaffin cells" = "#006e3b",
  "late Chromaffin cells" = "#686c58",
  "cycling Neuroblasts" = "#aacedc",
  "Neuroblasts" = "#244a94",
  "late Neuroblasts" = "#303d63"
)



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


# ## Jansky ----
# 
# metadata_jansky <-
#   readRDS("data_wip/tumor_data_jansky.rds")@meta.data %>%
#   as_tibble(rownames = "cell") %>% 
#   filter(anno_new == "Tumor cells") %>% 
#   select(cell, sample = patientID)
# 
# samples_jansky <-
#   read_csv("data_wip/metadata_samples_jansky.csv", comment = "#") %>% 
#   transmute(
#     sample,
#     time_point,
#     group = recode(clinical_subtype,
#                    MYCN = "T-M", LR = "low_risk", ALT = "T-A", TERT = "T-S"),
#     group = if_else(mesenchymal_features, "mesenchymal", group)
#   )
# 
# singler_results_jansky <-
#   read_csv("data_wip/adrmed_class_jansky.csv") %>% 
#   inner_join(metadata_jansky, by = "cell") %>% 
#   left_join(samples_jansky, by = "sample") %>% 
#   relocate(sample, group) %>% 
#   mutate(dataset = "jansky", .before = 1) %>% 
#   filter(time_point == "primary") %>% 
#   select(!time_point)


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
    # singler_results_jansky,
    singler_results_nb
  ) %>% 
  mutate(
    # group = fct_relevel(group, "low_risk", "mesenchymal"),
    mycn_status = if_else(group %in% c("M", "T-M"), "amplified", "normal"),
    tumor = if_else(group %in% c("A", "M", "S"), "DTC", "primary"),
    pruned_labels =
      pruned_labels %>% 
      factor(levels = names(ADRMED_CELLS_COLORS)) %>% 
      fct_rev() %>% 
      fct_explicit_na()
  )



# Analyze -----------------------------------------------------------------

# # this should resemble Jansky's figure 3d
# plot_jansky_like <- function() {
#   singler_results_jansky %>%
#     mutate(
#       cell_type =
#         pruned_labels %>% 
#         factor(levels = names(ADRMED_CELLS_COLORS)) %>% 
#         fct_rev(),
#       group = if_else(group %in% c("T-A", "T-S"), "TERT/ALT", group)
#     ) %>% 
#     group_by(group) %>% 
#     count(cell_type) %>% 
#     mutate(n_rel = n / sum(n) * 100) %>% 
#     ungroup() %>% 
#     mutate(group = fct_relevel(group, "low_risk", "TERT/ALT", "T-M")) %>% 
#     ggplot(aes(cell_type, n_rel)) +
#     geom_col(aes(fill = cell_type), show.legend = FALSE) +
#     scale_x_discrete(NULL, drop = FALSE) +
#     ylab("Percentage of cells") +
#     scale_fill_manual(values = ADRMED_CELLS_COLORS) +
#     coord_flip() +
#     facet_wrap(vars(group), nrow = 1, scales = "free_x") +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       strip.background = element_blank()
#     )
# }
# 
# plot_jansky_like()
# ggsave_default("comparison/jansky_3d", height = 80, width = 200)



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



# Publication figures -----------------------------------------------------

## Figure 2c ----

plot_adrmed_profile <- function() {
  singler_results %>% 
    group_by(tumor, mycn_status) %>% 
    count(cell_type = pruned_labels) %>% 
    mutate(n_rel = n / sum(n) * 100) %>% 
    ungroup() %>% 
    filter(cell_type != "(Missing)") %>% 
    mutate(mycn_status = fct_relevel(mycn_status, "normal")) %>% 
    ggplot(aes(cell_type, n_rel)) +
    geom_col(aes(fill = cell_type), show.legend = FALSE) +
    scale_x_discrete(NULL) +
    ylab("Percentage of cells") +
    scale_fill_manual(values = ADRMED_CELLS_COLORS) +
    coord_flip() +
    facet_wrap(vars(tumor, mycn_status), nrow = 1, scales = "free_x") +
    theme_nb(grid = FALSE) +
    theme()
}

plot_adrmed_profile()
ggsave_publication("2c_adrmed_profile", width = 8, height = 3)



## Figure S2b ----

plot_adrmed_heatmap <- function() {
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "pt"),
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(1, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(1, "mm")
  )
  
  mat <- 
    singler_results %>% 
    filter(dataset != "jansky") %>%
    group_by(sample) %>%
    count(cell_type = pruned_labels) %>% 
    mutate(n = n / sum(n) * 100) %>% 
    ungroup() %>% 
    filter(cell_type != "(Missing)") %>% 
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
      at = round(c(min(mat), max(mat)), 2),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    column_split =
      fct_cross(col_metadata$tumor, col_metadata$mycn_status) %>% 
      fct_relevel("DTC:amplified", "DTC:normal"),
    cluster_column_slices = FALSE,
    column_title = NULL,
    
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    row_dend_width = unit(3, "mm"),
    
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_dend_gp = gpar(lwd = 0.5),
    column_dend_height = unit(3, "mm"),
    
    width = unit(60, "mm"),
    height = unit(27, "mm"),
    column_gap = unit(0.5, "mm"),
    border = FALSE,
    
    top_annotation = HeatmapAnnotation(
      tumor = col_metadata$tumor,
      mycn = col_metadata$mycn_status,
      col = list(
        mycn = c("normal" = "gray90", "amplified" = "#d35f5f"),
        tumor = c(DTC = "black", primary = "grey80")
      ),
      annotation_label = list(
        mycn = "MYCN status",
        title = "tumor"
      ),
      annotation_name_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      annotation_legend_param = list(
        mycn = list(
          title = "MYCN status",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        ),
        tumor = list(
          title = "tumor",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  )
}

(p <- plot_adrmed_heatmap())
ggsave_publication("S2b_adrmed_heatmap", plot = p, height = 5, width = 12)
