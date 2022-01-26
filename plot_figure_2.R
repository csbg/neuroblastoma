library(scuttle)
library(monocle3)
library(ComplexHeatmap)
library(RColorBrewer)
library(scico)
library(latex2exp)
library(tidyverse)
source("common_functions.R")
source("styling.R")

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(1, "mm")
)



# Load data ---------------------------------------------------------------

## General ----

nb_metadata <- readRDS("data_generated/metadata.rds")

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")
colData(nb) <-
  nb_metadata %>%
  mutate(Size_Factor = colData(nb)$Size_Factor) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

selected_markers <- read_csv("metadata/cell_markers.csv", comment = "#")

tumor_data <- readRDS("data_generated/tumor_data_pseudobulk.rds")


## Adrenal medullary cell types ----

singler_results_dong <- 
  dir_ls("data_generated/adrmed", glob = "*class_dong*") %>% 
  map_dfr(read_csv, .id = "file") %>% 
  extract(file, into = "sample", regex = "dong_(.+)\\.") %>% 
  left_join(
    read_csv("metadata/samples_dong.csv", comment = "#") %>% 
      filter(high_risk) %>% 
      transmute(
        sample = tumor_id,
        group = if_else(mycn_amplified, "T-M", "T-S")
      ),
    by = "sample"
  ) %>% 
  relocate(group, .after = sample) %>% 
  mutate(dataset = "dong", .before = 1)

singler_results_nb <- 
  read_csv("data_generated/adrmed/adrmed_class_nb.csv") %>% 
  left_join(
    nb_metadata %>% 
      transmute(
        cell,
        sample = rename_patients(sample),
        group = rename_groups(group)
      ),
    by = "cell"
  ) %>% 
  relocate(sample, group) %>% 
  mutate(dataset = "nb", .before = 1)

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



# Figures -----------------------------------------------------------------

## 2a ----

plot_patient_dots <- function(min_exp = -2.5, max_exp = 2.5) {
  scale_and_limit <- function(x) {
    scale(x)[,1] %>% 
      pmax(min_exp) %>% 
      pmin(max_exp)
  }
  
  nb <- nb[, colData(nb)$cellont_abbr != "other"]
  colData(nb)$cellont_abbr <- fct_drop(colData(nb)$cellont_abbr)
  colData(nb)$cellont_cluster <- fct_drop(colData(nb)$cellont_cluster)
  
  features <- 
    selected_markers %>%
    filter(cell_type == "NB") %>%
    pull(gene)
  
  groups <-
    if_else(
      colData(nb)$cellont_abbr == "NB",
      as.character(colData(nb)$sample),
      as.character(colData(nb)$cellont_abbr)
    ) %>% 
    rename_patients() %>% 
    as_factor() %>% 
    fct_relevel(PATIENT_ORDER, names(CELL_TYPE_ABBREVIATIONS))
  
  vis_data <- 
    logcounts(nb)[features, ] %>% 
    Matrix::t() %>% 
    as.matrix() %>% 
    as_tibble(rownames = "cell") %>% 
    group_by(id = groups) %>% 
    summarise(
      across(
        where(is.numeric),
        list(
          avg_exp = ~mean(expm1(.)),
          pct_exp = ~length(.[. > 0]) / length(.) * 100
        ),
        .names = "{.col}__{.fn}"
      )
    ) %>% 
    mutate(across(ends_with("avg_exp"), scale_and_limit)) %>%
    pivot_longer(
      !id,
      names_to = c("feature", ".value"),
      names_pattern = "(.+)__(.+)"
    ) %>% 
    mutate(
      feature = factor(feature, levels = features),
      type = if_else(
        id %in% PATIENT_ORDER,
        "cancer cells",
        "microenvironment cells"
      )
    )
  
  ggplot(vis_data, aes(id, feature)) +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_x_discrete("patient or cell type") +
    scale_y_discrete(NULL) +
    scale_color_dotplot(
      "scaled average expression",
      limits = c(-0.5, 2.5),
      breaks = c(-0.5, 2.5),
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(10, "mm"),
        label.position = "top",
        ticks = FALSE,
        title.vjust = 0.2
      )
    ) +
    scale_radius(
      "% expressed",
      range = c(0, 2), 
      limits = c(0, 100),
      breaks = c(0, 50, 100)
    ) +
    # coord_fixed() +
    facet_grid(
      cols = vars(type),
      scales = "free_x",
      space = "free_x"
    ) +
    theme_nb(grid = FALSE) +
    theme(
      legend.box.just = "bottom",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = "bottom",
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(-5, 1, 0, 1, "mm"),
      panel.spacing = unit(1, "mm")
      # plot.margin = margin(0, 1, 0, 1, "mm"),
    )
}

plot_patient_dots()
ggsave_publication("2a_nb_markers_patients", width = 7.5, height = 5)



## 2b ----

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
ggsave_publication("2b_adrmed_profile", width = 8, height = 3)



## 2c ----

plot_corr_mat <- function(size, samples = NULL) {
  if (is.null(samples)) {
    corr_mat <-
      tumor_data$pseudobulk_counts %>%
      magrittr::extract(tumor_data$highly_variable_genes, ) %>%
      cor(use = "pairwise.complete.obs")
  } else {
    corr_mat <-
      tumor_data$pseudobulk_counts %>%
      magrittr::extract(tumor_data$highly_variable_genes, samples) %>%
      cor(use = "pairwise.complete.obs")
  }
  
  distance <- as.dist(1 - corr_mat)
  
  col_metadata <- tibble(
    sample = colnames(corr_mat),
    group =
      str_sub(sample, 1, 1) %>% 
      fct_relevel("M", "A", "S", "T"),
    mycn_status = if_else(
      sample %in% c("T162", "T200", "T230") | str_starts(sample, "M"),
      "amplified",
      "normal"
    )
  )
  
  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(min(corr_mat), max(corr_mat), length.out = 9),
      scico(9, palette = "davos", direction = -1),
    ),
    name = "correlation of\npseudobulk\nexpression",
    heatmap_legend_param = list(
      at = round(c(min(corr_mat), max(corr_mat)), 2),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    clustering_distance_rows = distance,
    clustering_distance_columns = distance,
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    row_dend_width = unit(3, "mm"),
    
    width = unit(size, "mm"),
    height = unit(size, "mm"),
    border = FALSE,
    
    show_column_dend = FALSE,
    show_column_names = FALSE,
    
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
          title = "group",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        ),
        mycn = list(
          title = "MYCN status",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    ),
  )
}

(p <- plot_corr_mat(25, PATIENT_ORDER[-1:-5]))
ggsave_publication("2c_pseudobulk_cor", plot = p, height = 3, width = 8)



## S2a ----

plot_gene_signature <- function() {
  vis_data <- 
    nb_metadata %>% 
    mutate(
      cells =
        case_when(
          cellont_cluster == "NB (8)" ~ "tumor",
          TRUE ~ "other"
        ) %>% 
        fct_rev()
    ) %>%
    pivot_longer(
      starts_with("signature"),
      names_to = "signature",
      names_prefix = "signature_"
    ) %>% 
    mutate(
      signature =
        factor(signature) %>% 
        fct_recode("NCC-like" = "ncc_like") %>% 
        fct_relevel("adrenergic", "noradrenergic")
    )
  
  label_data <-
    vis_data %>% 
    distinct(signature) %>% 
    mutate(label = signature, cells = "other", value = 0.95)
  
  ggplot(vis_data, aes(cells, value)) +
    geom_violin(
      aes(fill = cells),
      scale = "width",
      size = BASE_LINE_SIZE
    ) +
    stat_summary(geom = "point", fun = mean, size = .1) +
    geom_text(
      data = label_data,
      aes(label = label),
      size = BASE_TEXT_SIZE_MM,
      hjust = 1
    ) +
    xlab(NULL) +
    scale_y_continuous(
      name = "gene signature score",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      expand = expansion(add = .025)
    ) +
    scale_fill_manual(
      values = c(CELL_TYPE_COLORS["NB"], "gray80")
    ) +
    coord_flip() +
    facet_grid(vars(signature)) +
    theme_nb(grid = FALSE) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey92",
        size = BASE_LINE_SIZE
      ),
    )
}

plot_gene_signature()
ggsave_publication("S2a_gene_signatures", width = 4, height = 4)



## S2b ----

plot_mesenchymal <- function(top_prop = 0.05) {
  data_highlight <-
    bind_rows(
      mesenchymal =
        nb_metadata %>%
        slice_max(signature_mesenchymal, prop = top_prop) %>% 
        select(
          cell, umap_1_monocle, umap_2_monocle,
          value = signature_mesenchymal
        ),
      `NCC-like` =
        nb_metadata %>%
        slice_max(signature_ncc_like, prop = top_prop) %>% 
        select(
          cell, umap_1_monocle, umap_2_monocle,
          value = signature_ncc_like
        ),
      .id = "signature"
    ) %>% 
    arrange(value)
  
  color_breaks <- round(
    c(min(data_highlight$value), max(data_highlight$value)),
    2
  )
  
  nb_metadata %>%
    ggplot(aes(umap_1_monocle, umap_2_monocle)) +
    geom_point(color = "gray90", size = 0.001, shape = 16) +
    geom_point(
      data = data_highlight,
      aes(color = value),
      size = 0.001,
      shape = 16
    ) +
    scale_x_continuous("UMAP1", breaks = c(-10, 0, 10)) +
    scale_y_continuous("UMAP2", breaks = c(-10, 0, 10)) +
    scale_color_distiller(
      name = "gene\nsignature\nscore",
      palette = "RdPu",
      direction = 1,
      breaks = c(min(data_highlight$value), max(data_highlight$value)),
      labels = function(x) round(x, 2),
      guide = guide_colorbar(
        barwidth = unit(2, "mm"),
        barheight = unit(15, "mm")
      )
    ) +
    coord_fixed() +
    facet_wrap(vars(signature)) +
    theme_nb(grid = FALSE) +
    theme(
      legend.position = c(.95, .32)
    )
}

plot_mesenchymal()
ggsave_publication("S2b_signature_umap", type = "png", width = 9, height = 5)



## S2c ----

plot_nb_dots <- function(top_prop = 0.05, min_exp = -2.5, max_exp = 2.5) {
  scale_and_limit <- function(x) {
    scale(x)[,1] %>% 
      pmax(min_exp) %>% 
      pmin(max_exp)
  }
  
  selected_cells <-
    bind_rows(
      NB =
        colData(nb) %>% 
        as_tibble(rownames = "cell") %>% 
        filter(cellont_abbr == "NB"),
      mesenchymal =
        colData(nb) %>% 
        as_tibble(rownames = "cell") %>% 
        slice_max(prop = top_prop, order_by = signature_mesenchymal),
      ncc_like =
        colData(nb) %>% 
        as_tibble(rownames = "cell") %>% 
        slice_max(prop = top_prop, order_by = signature_ncc_like),
      .id = "type"
    ) %>%  
    mutate(
      type = case_when(
        type == "NB" | cluster_50 == "8" ~ "NB",
        TRUE ~ str_c(type, cluster_50, sep = "__")
      )
    )
  
  nb_markers <-
    selected_markers %>% 
    filter(cell_type == "NB") %>% 
    pull(gene)
  
  vis_data <- 
    logcounts(nb)[nb_markers, selected_cells$cell] %>% 
    Matrix::t() %>% 
    as.matrix() %>%
    as_tibble(rownames = "cell") %>%
    group_by(type = selected_cells$type) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          avg_exp = ~mean(expm1(.)),
          pct_exp = ~length(.[. > 0]) / length(.) * 100
        ),
        .names = "{.col}__{.fn}"
      )
    ) %>% 
    mutate(across(ends_with("avg_exp"), scale_and_limit)) %>%
    pivot_longer(
      !type,
      names_to = c("feature", ".value"),
      names_pattern = "(.+)__(.+)"
    ) %>%
    separate(type, into = c("type", "cluster"), sep = "__", fill = "right") %>% 
    replace_na(list(cluster = "8")) %>% 
    mutate(
      feature = factor(feature, levels = nb_markers),
      type = fct_relevel(type, "NB") %>% fct_recode("NCC-like" = "ncc_like"),
      cluster = fct_inseq(cluster)
    )
  
  ggplot(vis_data, aes(cluster, feature)) +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_x_discrete(
      "cluster",
      expand = expansion(add = 0.5)
    ) +
    scale_y_discrete(
      NULL,
      expand = expansion(add = 0.5)
    ) +
    scale_color_dotplot(
      "scaled average expression",
      limits = c(-0.5, 2.5),
      breaks = c(-0.5, 2.5),
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        label.position = "top",
        ticks = FALSE,
        title.vjust = 0.2
      )
    ) +
    scale_radius(
      "% expressed",
      range = c(0, 2), 
      limits = c(0, 100)
    ) +
    facet_grid(
      cols = vars(type),
      scales = "free_x",
      space = "free_x"
    ) +
    theme_nb(grid = FALSE) +
    theme(
      legend.box.just = "bottom",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = "bottom",
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(-5, 1, 0, 1, "mm"),
      panel.spacing = unit(1, "mm")
      # plot.margin = margin(0, 1, 0, 1, "mm"),
    )
}

plot_nb_dots()
ggsave_publication("S2c_signature_dots", width = 9, height = 4)



## S2d ----

plot_adrmed_heatmap <- function() {
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
ggsave_publication("S2d_adrmed_heatmap", plot = p, height = 5, width = 12)



## S2e ----

(p <- plot_corr_mat(40))
ggsave_publication("S2e_pseudobulk_cor_all", plot = p, height = 4.5, width = 8)

