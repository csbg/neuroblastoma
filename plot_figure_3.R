# @DEPI metadata.rds
# @DEPI dge_mm_results.rds

library(scuttle)
library(monocle3)
library(ComplexHeatmap)
library(RColorBrewer)
library(scico)
library(latex2exp)
library(tidyverse)
library(patchwork)
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

nb_metadata <- readRDS("data_generated/metadata.rds")

dge <- readRDS("data_generated/dge_mm_results.rds")



# Figures -----------------------------------------------------------------

## 3a ----

nb_metadata %>%
  group_by(group, sample) %>% 
  count(cell_type = cellont_abbr) %>% 
  mutate(
    n = n / sum(n) * 100,
    sample = rename_patients(sample),
    group = rename_groups(group)
  ) %>%
  complete(nesting(group, sample), cell_type, fill = list(n = 0)) %>% 
  ggplot(aes(cell_type, n, color = group)) +
  geom_boxplot(outlier.shape = NA, size = .25, show.legend = FALSE) +
  geom_point(position = position_jitterdodge(seed = 1), size = .5, alpha = .5) +
  xlab("cell type") +
  ylab("relative abundance (%)") +
  scale_color_manual(
    values = GROUP_COLORS,
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  theme_nb(grid = FALSE) +
  theme(
    legend.key.height = unit(1, "mm"),
    legend.key.width = unit(1, "mm"),
    legend.position = c(.95, .8)
  )
ggsave_publication("3a_cell_type_abundances", width = 13, height = 4)



## 3b ----

plot_consistent_genes <- function() {
  consistent_genes <- 
    dge$results_wide_filtered %>%
    mutate(comparison = rename_contrast(comparison)) %>% 
    filter(
      abs(logFC) > log(2),
      comparison %>% str_ends("c")
    ) %>%
    group_by(cell_type, direction, comparison) %>%
    summarise(genes = list(unique(gene))) %>%
    ungroup() %>% 
    pivot_wider(names_from = comparison, values_from = genes) %>% 
    rowwise() %>%
    mutate(
      shared_AM = Ac %>% intersect(Mc) %>% length(),
      shared_AS = Ac %>% intersect(Sc) %>% length(),
      shared_MS = Mc %>% intersect(Sc) %>% length(),
      shared_AMS = Ac %>% intersect(Mc) %>% intersect(Sc) %>% length()
    ) %>%
    ungroup() %>%
    mutate(
      across(
        starts_with("shared"),
        ~if_else(direction == "up", ., -.)
      )
    ) %>% 
    pivot_longer(
      starts_with("shared"),
      names_to = "shared",
      names_prefix = "shared_",
      values_to = "n"
    ) %>% 
    mutate(
      shared =
        shared %>% 
        fct_relevel("AM", "AS", "MS") %>% 
        fct_recode("A and M" = "AM", "A and S" = "AS",
                   "M and S" = "MS", "A, M, and S" = "AMS"),
      cell_type = factor(cell_type, levels = names(CELL_TYPE_ABBREVIATIONS))
    )
  
  consistent_genes %>% 
    ggplot(aes(cell_type, n, fill = shared)) +
    geom_col() +
    geom_hline(yintercept = 0, size = BASE_LINE_SIZE) +
    scale_fill_manual(
      "shared\nbetween",
      values = CONSISTENT_GENES_COLORS,
      guide = guide_legend(reverse = TRUE)
    ) +
    xlab("cell type") +
    ylab("number of consistently\ndown- and upregulated genes") +
    theme_nb(grid = FALSE) +
    theme(
      axis.ticks.length.x = unit(0, "mm"),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      legend.margin = margin(0, 0, 0, -2, "mm"),
      panel.border = element_blank(),
      panel.grid.major.y = element_line(
        color = "grey92",
        size = BASE_LINE_SIZE
      )
    )
}

plot_consistent_genes()
ggsave_publication("3b_number_of_consistent_genes", width = 5, height = 4)



## 3c ----

plot_logfc_correlation_heatmap <- function() {
  corr_mat <- 
    dge$results_wide %>% 
    filter(abs(logFC) < 8, comparison  %>% str_ends("vs_I")) %>%
    select(gene, cell_type, comparison, logFC) %>%
    mutate(comparison = rename_contrast(comparison)) %>%
    unite(cell_type, comparison, col = "cluster_group", sep = "_") %>%
    pivot_wider(names_from = "cluster_group", values_from = "logFC") %>%
    select(!gene) %>%
    cor(use = "pairwise.complete.obs")
  
  distance <- as.dist(1 - corr_mat)
  
  metadata <- 
    colnames(corr_mat) %>% 
    str_match("(.+)_(.+)")
  
  group_names <- 
    metadata[, 3] %>% 
    factor(levels = names(CONTRAST_COLORS)) %>% 
    fct_drop()
  
  colnames(corr_mat) <- metadata[, 2]
  rownames(corr_mat) <- metadata[, 2]
  
  color_max <- round(max(corr_mat[lower.tri(corr_mat)]), 2)
  
  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(0, color_max, length.out = 9),
      scico(9, palette = "davos", direction = -1),
    ),
    heatmap_legend_param = list(
      at = c(0, color_max),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    name = "log fold change\ncorrelation",
    
    clustering_distance_rows = distance,
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    row_dend_width = unit(3, "mm"),
    
    clustering_distance_columns = distance,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    
    height = unit(4.5, "cm"),
    width = unit(4.5, "cm"),
    border = FALSE,
    
    right_annotation = rowAnnotation(
      group = group_names,
      col = list(group = CONTRAST_COLORS),
      show_annotation_name = FALSE,
      show_legend = TRUE,
      annotation_legend_param = list(
        group = list(
          title = "comparison",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  )
}

(p <- plot_logfc_correlation_heatmap())
ggsave_publication("3c_logfc_correlation_heatmap",
                   plot = p, width = 7, height = 5)



## 3d ----

plot_gsea <- function(db = "MSigDB_Hallmark_2020",
                      comparisons = NULL,
                      cell_types = NULL,
                      top_n_positive = 5L,
                      top_n_negative = 5L,
                      max_p_adj = 0.05,
                      min_abs_NES = 1) {
  data <- 
    dge$gsea %>%
    mutate(comparison = factor(comparison) %>% rename_contrast()) %>% 
    filter(db == {{db}})
  
  if (!is.null(comparisons))
    data <- 
      data %>% 
      filter(comparison %in% {{comparisons}})
  
  if (!is.null(cell_types))
    data <- 
      data %>% 
      filter(cell_type %in% {{cell_types}})
  
  data_top_terms <-
    data %>% 
    filter(
      padj <= max_p_adj,
      abs(NES) >= min_abs_NES
    ) %>% 
    group_by(comparison, cell_type)
  
  top_terms_pos <- 
    data_top_terms %>% 
    slice_max(n = top_n_positive, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  top_terms_neg <- 
    data_top_terms %>% 
    slice_min(n = top_n_negative, order_by = NES, with_ties = FALSE) %>%
    pull(pathway) %>%
    unique()
  
  data_vis <- 
    data %>% 
    filter(pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE),
      cell_type = factor(cell_type, levels = names(CELL_TYPE_ABBREVIATIONS))
    )
  
  if (nlevels(data_vis$pathway) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$pathway), 5),
        size = BASE_LINE_SIZE,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  color_limit <- max(abs(dge$gsea$NES), na.rm = TRUE)
  
  ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    xlab(NULL) +
    ylab(NULL) +
    scale_color_gsea(
      "normalized\nenrichment\nscore",
      limits = c(-color_limit, color_limit),
      breaks = c(-color_limit, 0, color_limit),
      labels = function(x) round(x, 2),
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_size_area(
      name = TeX("-log_{10} p_{adj}"),
      max_size = 2.5,
      limits = c(0, 40),
      breaks = c(0, 20, 40),
      labels = c("0", "20", "40 or higher"),
      oob = scales::oob_squish
    )  +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    theme_nb(grid = FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.margin = margin(),
      panel.spacing = unit(-.5, "pt"),
    )
}

wrap_plots(
  plot_gsea(comparisons = c("M vs C", "A vs C", "S vs C"), cell_types = "B"),
  plot_gsea(comparisons = c("M vs C", "A vs C", "S vs C"), cell_types = "M"),
  plot_gsea(comparisons = "M vs A+S", cell_types = "B"),
  plot_gsea(comparisons = "M vs A+S", cell_types = "M"),
  wrap_elements(
    panel = textGrob(
      "comparison",
      just = "center",
      gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    clip = FALSE
  )
) +
  plot_layout(
    design = "ABCD\nEEEE",
    widths = c(5, 5, 3, 3),
    heights = c(30, 1),
    guides = "collect"
  )
ggsave_publication("3d_gsea", width = 18, height = 6)



## 3e ----

selected_genes <- list(
  "TNF-alpha Signaling via NF-kB" = c("IL1B", "FOS", "NFKB1", "MYC", "IFNGR2",
                                      "CD44", "RELB", "SOD2", "MAP3K8"),
  "Interferon Gamma Response" = c("IFI44", "NFKB1", "MYD88", "HLA-A", "CD86",
                                  "HLA-DQA1", "HLA-DQB1", "TNFSF10"),
  "Myc Targets V1" = c("TRIM28", "PCNA", "CDK4", "MCM5", "HDAC2"),
  "E2F Targets" = c("DNMT1", "EZH2", "MKI67", "PCNA")
)

plot_pathway_genes <- function(db = "MSigDB_Hallmark_2020",
                               pathways = selected_genes,
                               level = c("group", "sample"),
                               cell_types = c("B", "M")) {
  level <- match.arg(level)
  
  row_metadata <-
    dge$gene_sets[[db]][names(pathways)] %>%
    enframe("pathway", "gene") %>%
    unnest_longer(gene) %>%
    mutate(pathway = as_factor(pathway)) %>%
    filter(gene %in% rownames(dge$cds)) %>% 
    semi_join(
      selected_genes %>% 
        enframe("pathway", "gene") %>% 
        unnest_longer(gene)
    )
  
  # collect barcodes per column
  col_metadata <-
    dge$metadata %>% 
    filter(
      cellont_abbr %in% {{cell_types}},
      cell %in% colnames(dge$cds)
    ) %>% 
    mutate(
      cell_type = factor(cellont_abbr, names(CELL_TYPE_ABBREVIATIONS)),
      group = rename_groups(group),
      sample = rename_patients(sample) %>% factor(levels = PATIENT_ORDER)
    ) %>% 
    arrange(group) %>% 
    {
      if (level == "sample")
        group_by(., cell_type, group, sample)
      else
        group_by(., cell_type, group)
    } %>% 
    summarise(cells = list(cell))
  
  # make expression matrix
  mat <-
    dge$cds %>% 
    logcounts() %>% 
    magrittr::extract(row_metadata$gene, ) %>% 
    as.matrix()
  
  mat <- 
    col_metadata %>% 
    pull(cells) %>% 
    map(~mat[, .] %>% rowMeans()) %>%
    purrr::reduce(cbind)
  
  # scale within cell types
  if (level == "sample") {
    mat <- cbind(
      mat[, 1:16] %>% t() %>% scale() %>% t(),
      mat[, 17:32] %>% t() %>% scale() %>% t()
    )
  } else {
    mat <- cbind(
      mat[, 1:4] %>% t() %>% scale() %>% t(),
      mat[, 5:8] %>% t() %>% scale() %>% t()
    )
  }
  
  # plot heatmap
  Heatmap(
    mat,
    name = "expression",
    col = circlize::colorRamp2(
      seq(
        min(mat),
        max(mat),
        length.out = 9
      ),
      viridisLite::cividis(9)
    ),
    border = FALSE,
    heatmap_legend_param = list(
      at = c(min(mat), max(mat)),
      labels = c("low", "high"),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    show_row_names = TRUE,
    row_split = row_metadata$pathway,
    row_title_rot = 0,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_dend = TRUE,
    row_gap = unit(.5, "mm"),
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    cluster_columns = FALSE,
    column_split = col_metadata$cell_type,
    show_column_names = FALSE,
    column_gap = unit(.5, "mm"),
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    width = unit(.5 * (length(cell_types) - 1) + ncol(mat) * 2, "mm"),
    height = unit(.5 * (length(pathways) - 1) + nrow(mat) * 2, "mm"),
    
    top_annotation = HeatmapAnnotation(
      group = col_metadata$group,
      col = list(
        group = GROUP_COLORS
      ),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        group = list(
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  )
}

(p <- plot_pathway_genes())
ggsave_publication("3e_pathway_genes",
                   plot = p, width = 9, height = 6)



## S3a ----

plot_celltype_enrichment <- function(p_lim = 20, or_lim = 3) {
  do_fisher_test <- function(sample, cellont_abbr) {
    res <- 
      nb_metadata %>%
      filter(group == "I" | sample == .env$sample) %>%
      mutate(
        sample = fct_collapse(
          sample,
          yes = .env$sample,
          other_level = "no"
        ),
        cell_type = fct_collapse(
          cellont_abbr,
          yes = .env$cellont_abbr,
          other_level = "no"
        )
      ) %>%
      count(sample, cell_type) %>%
      pivot_wider(names_from = cell_type, values_from = n) %>%
      column_to_rownames("sample") %>%
      as.matrix() %>%
      fisher.test()
    
    list(
      odds_ratio = res$estimate,
      p_val = res$p.value
    )
  }
  
  nb_metadata %>%
    distinct(
      group,
      sample = as.character(sample),
      cell_type = as.character(cellont_abbr)
    ) %>%
    filter(group != "I", cell_type != "NB") %>%
    rowwise() %>%
    mutate(fisher = list(do_fisher_test(sample, cell_type))) %>%
    unnest_wider(fisher) %>% 
    mutate(p_adj = p.adjust(p_val, method = "BH")) %>% 
    filter(cell_type != "other") %>% 
    mutate(
      cell_type =
        as_factor(cell_type) %>%
        fct_relevel(
          levels(nb_metadata$cellont_cluster) %>%
            str_extract("\\w+") %>% 
            unique()
        ),
      log_p_adj = pmin(-log10(p_adj), p_lim),
      log_odds_ratio = log2(odds_ratio),
      log_odds_ratio = case_when(
        log_odds_ratio > or_lim  ~ or_lim,
        log_odds_ratio < -or_lim ~ -or_lim,
        TRUE                     ~ log_odds_ratio
      ),
      sample =
        sample %>% 
        factor(levels = levels(nb_metadata$sample)) %>% 
        rename_patients() %>% 
        fct_rev(),
      group = rename_groups(group)
    ) %>% 
    ggplot(aes(cell_type, sample)) + 
    geom_point(aes(color = log_odds_ratio, size = log_p_adj), shape = 16) + 
    xlab("cell type") +
    ylab("patient") +
    scale_color_gsea(
      name = "enrichment relative\nto control\n(log2 odds ratio)",
      limits = c(-or_lim, or_lim),
      breaks = c(-or_lim, 0, or_lim),
      labels = c(
        str_glue("-{or_lim} or lower"),
        "0",
        str_glue("{or_lim} or higher")
      ),
      guide = guide_colorbar(ticks = FALSE)
    ) +
    scale_radius(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0.25, 3),
      breaks = c(0, 10, 20),
      labels = function(x) c(x[-length(x)], paste(x[length(x)], "or higher"))
    ) +
    facet_grid(vars(group), scales = "free_y", space = "free_y") +
    theme_nb(grid = FALSE) +
    theme(
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      legend.spacing = unit(5, "mm"),
      legend.margin = margin(0, 0, 0, 0, "mm"),
      panel.spacing = unit(-.5, "pt"),
      strip.text.y.right = element_text(angle = 0)
    )
}

plot_celltype_enrichment()
ggsave_publication("S3a_cell_type_enrichment", width = 5, height = 4)



## S3b ----

plot_gsea(comparisons = c("M vs C", "A vs C", "S vs C")) +
  xlab("comparison")
ggsave_publication("S3b_gsea_all_vs_C", width = 10, height = 10)



## S3c ----

plot_gsea(comparisons = "M vs A+S") +
  xlab("M vs A+S\ncomparison") +
  theme(
    axis.ticks.length.x = unit(0, "mm"),
    axis.text.x = element_blank(),
    strip.text = element_text(angle = 90, hjust = 0)
  )
ggsave_publication("S3c_gsea_all_vs_AS", width = 7, height = 7)



## S3d ----

plot_pathway_genes_sc <- function(db = "MSigDB_Hallmark_2020",
                                  pathways = c(
                                    "TNF-alpha Signaling via NF-kB" = "up",
                                    "Interferon Gamma Response" = "up",
                                    "Myc Targets V1" = "down",
                                    "E2F Targets" = "down"
                                  )) {
  log_fold_changes <- 
    dge$results_wide_filtered %>% 
    filter(
      p_adj <= 0.05,
      comparison %in% c("II_vs_I", "III_vs_I", "IV_vs_I")
    ) %>%
    group_by(gene) %>% 
    summarise(min_logFC = min(logFC), max_logFC = max(logFC))
  
  row_metadata <-
    left_join(
      pathways %>% 
        enframe("pathway", "direction"),
      dge$gene_sets[[db]][names(pathways)] %>%
        enframe("pathway", "gene"),
      by = "pathway"
    ) %>% 
    unnest_longer(gene) %>%
    mutate(pathway = as_factor(pathway)) %>%
    filter(gene %in% rownames(dge$cds)) %>%
    left_join(log_fold_changes, by = "gene") %>% 
    filter(
      direction == "up" & max_logFC > log(2) |
      direction == "down" & min_logFC < -log(2)
    )
  
  set.seed(1)
  col_metadata <-
    dge$metadata %>% 
    filter(
      cellont_abbr %in% c("B", "M"),
      cell %in% colnames(dge$cds)
    ) %>% 
    group_by(cellont_abbr, group) %>%
    slice_sample(n = 300) %>%
    # group_by(cellont_abbr, group, sample) %>% 
    # slice_sample(n = 150) %>% 
    # group_by(cellont_abbr, group) %>%
    # slice_sample(prop = 1) %>%
    mutate(
      cell_type = factor(cellont_abbr, names(CELL_TYPE_ABBREVIATIONS)),
      group = rename_groups(group),
      sample = rename_patients(sample) %>% factor(levels = PATIENT_ORDER)
    ) %>% 
    arrange(cell_type, group)
  
  mat <-
    dge$cds %>% 
    logcounts() %>% 
    magrittr::extract(row_metadata$gene, col_metadata$cell) %>% 
    as.matrix() %>% 
    t() %>% scale() %>% t() %>%
    replace_na(min(., na.rm = TRUE))
  
  Heatmap(
    mat,
    name = "expression",
    col = circlize::colorRamp2(
      seq(0, 1.5, length.out = 9),
      # scico(9, palette = "oslo", direction = -1)
      viridisLite::cividis(9)
    ),
    use_raster = FALSE,
    # raster_quality = 3,
    
    heatmap_legend_param = list(
      at = c(0, 1.5),
      labels = c("low", "high"),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    border = FALSE,
    show_row_names = FALSE,
    row_split = row_metadata$pathway,
    row_title_rot = 0,
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    row_gap = unit(.5, "mm"),
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    cluster_columns = FALSE,
    column_split = col_metadata$cell_type,
    show_column_names = FALSE,
    column_gap = unit(.5, "mm"),
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    top_annotation = HeatmapAnnotation(
      group = col_metadata$group,
      col = list(
        group = GROUP_COLORS
      ),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        group = list(
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  )
}

(p <- plot_pathway_genes_sc())
ggsave_publication("S3d_pathway_genes_sc", type = "png",
                   plot = p, width = 9, height = 5)



# Tables ------------------------------------------------------------------

## S4 ----

gene_pathways <- 
  CellChatDB.human$interaction %>% 
  as_tibble() %>% 
  select(interaction_name, pathway_name) %>% 
  mutate(genes = str_split(interaction_name, "_")) %>% 
  unnest_longer(genes) %>% 
  distinct(pathway_name, genes) %>% 
  group_by(Gene = genes) %>% 
  summarise(
    Pathways =
      unique(pathway_name) %>%
      str_sort() %>%
      str_c(collapse = ", ")
  )

used_comparisons <- c(
  "II_vs_I",
  "III_vs_I",
  "IV_vs_I",
  "MNA_vs_other"
)

rename_columns <- function(s) {
  comparison <- str_sub(s, 7)
  case_when(
    str_starts(s, "log") ~ paste0("Log fold change (", comparison, ")"),
    str_starts(s, "p_adj") ~ paste0("Adjusted p-value (", comparison, ")"),
    TRUE ~ s
  )
}

dge$results_wide_filtered %>% 
  arrange(comparison, cell_type, desc(logFC)) %>%
  filter(comparison %in% used_comparisons) %>% 
  mutate(
    logFC = logFC / log(2),  # nebula returns natural log fold changes
    comparison = rename_contrast_long(comparison),
    cell_type = CELL_TYPE_ABBREVIATIONS[cell_type]
  ) %>%
  select(!c(p, frq, frq_ref, direction)) %>%
  pivot_wider(names_from = comparison, values_from = c(logFC, p_adj)) %>%
  relocate(1, 2, 3, 7, 4, 8, 5, 9, 6, 10) %>% 
  rename_with(rename_columns) %>% 
  left_join(gene_pathways, by = c(gene = "Gene")) %>%
  split(.$cell_type) %>%
  map(select, !cell_type) %>%
  save_table("S4_dge")



## S5 ----

dge$gsea %>% 
  arrange(db, comparison, cell_type, desc(NES)) %>%
  filter(comparison %in% used_comparisons) %>% 
  mutate(
    comparison = rename_contrast_long(comparison),
    cell_type = CELL_TYPE_ABBREVIATIONS[cell_type],
    leadingEdge = map_chr(leadingEdge, str_c, collapse = ", ")
  ) %>%
  select(
    "Database" = db,
    "Comparison" = comparison,
    "Cell type" = cell_type,
    "Pathway" = pathway,
    "Normalized Enrichment Score" = NES,
    "Adjusted p-value" = padj,
    "Leading edge genes" = leadingEdge
  ) %>%
  split(.$`Cell type`) %>%
  map(select, !`Cell type`) %>%
  save_table("S5_gsea")
