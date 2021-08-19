# Plot mixed model DGE results.
#
# @DEPI dge_mm_results.rds
# @DEPI dge_pb_results.rds

library(scater)
library(tidyverse)
library(latex2exp)
library(ComplexHeatmap)
library(ggrepel)
library(scico)
library(ggpmisc)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

dge <- readRDS("data_generated/dge_mm_results.rds")
dge_pb <- readRDS("data_generated/dge_pb_results.rds")

# are there problematic genes with convergence <= -20 ?
any(dge$results$convergence <= -20)
any(dge$results_tumor$convergence <= -20)



# Plot results ------------------------------------------------------------

## Volcano plots ----

ggplot(dge$results_wide, aes(logFC, -log10(p))) +
  geom_point(size = 0.1) +
  facet_grid(vars(group), vars(cell_type)) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("dge_mm/volcano")


dge$results_wide %>% 
  slice_sample(prop = 1) %>% 
  ggplot(aes(logFC, -log10(p))) +
  geom_point(aes(color = cell_type), size = 0.1) +
  scale_color_manual(values = CELL_TYPE_COLORS) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("dge_mm/volcano_all")



## Violin plots ----

plot_violin <- function(gene, cluster, left_group) {
  plotExpression(
    dge$cds[, dge$cds$cellont_abbr == cluster &
              dge$cds$group %in% c(left_group, "I")],
    gene,
    x = "sample",
    colour_by = "group"
  )  
}

dge$results_wide %>% 
  arrange(desc(logFC)) %>% 
  filter(abs(logFC) > 20, cell_type == "T")

# high absolute logFC
# ... unfiltered
plot_violin("CKAP4", "pDC", "II")
plot_violin("KIF6", "pDC", "II")
plot_violin("PRF1", "E", "II")
plot_violin("MYCN", "T", "II")
plot_violin("STMN2", "T", "II")
plot_violin("XIST", "T", "III")

# ... after filtering by group frequency
plot_violin("RPS4Y1", "T", "II")
plot_violin("IGLC3", "T", "II")

# ... after filtering by minimum sample frequency
plot_violin("IFI44L", "SC", "II")
plot_violin("MTRNR2L8", "T", "III")
plot_violin("MAP3K8", "M", "IV")

# ... downregulated genes
plot_violin("JUND", "E", "II")
plot_violin("TCL1A", "B", "III")
plot_violin("HIST1H1E", "B", "IV")



## Comparison to pseudobulk ----

plot_comparison <- function(lim = NULL, filename = NULL) {
  p <-
    dge_pb$results %>% 
    filter(contrast != "tif") %>% 
    select(gene, cell_type = cluster_id,
           contrast, logFC_pb = logFC, p_pb = p_val) %>% 
    extract(contrast, into = "group", regex = "(.+)_vs_I") %>% 
    left_join(
      dge$results_wide %>% rename(logFC_mm = logFC, p_mm = p),
      by = c("gene", "cell_type", "group")
    ) %>% 
    ggplot(aes(logFC_pb, logFC_mm)) +
    geom_point(size = 0.1) +
    facet_grid(vars(group), vars(cell_type)) +
    coord_fixed(xlim = lim, ylim = lim) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  ggsave_default(filename)
  p
}

plot_comparison(lim = NULL, filename = "dge_mm/comparison_pb_mm_full")
plot_comparison(lim = c(-10, 10), filename = "dge_mm/comparison_pb_mm")



## Enrichr results ----

#' Generate a dotplot for enrichment terms.
#' 
#' Dots that meet given limits of adjusted p value, overlap size, and odds ratio
#' are bordered.
#'
#' @param data Results from `enrich_all_genes()`.
#' @param db Enrichr database for which results should be plotted.
#' @param direction Plot results for `"up"`- or `"down"`-regulated genes.
#' @param top_n Only plot the terms with the highest odds ratio per cluster.
#' @param max_p_adj Maximum adjusted p value, …
#' @param min_odds_ratio … minimum odds ratio, and …
#' @param min_overlap_size … minimum overlap size required for bordered dots.
#' @param log_odds_cap Upper boundary of the color scale.
#' @param filename Name of output file. If `"auto"`, derive from database name.
#' @param ... Additional parameters passed to `ggsave_default()`.
#'
#' @return A ggplot object.
plot_enrichr_dots <- function(data,
                              db,
                              direction = "up",
                              top_n = 5,
                              max_p_adj = 0.05,
                              min_odds_ratio = 5,
                              min_overlap_size = 1,
                              log_odds_cap = 2,
                              filename = "auto",
                              ...) {
  data_selected <-
    data %>% 
    filter(
      db == {{db}},
      direction == {{direction}},
      Odds.Ratio >= 1
    )
  
  top_terms <- 
    data_selected %>% 
    group_by(group, cell_type) %>%
    slice_max(n = top_n, order_by = Odds.Ratio, with_ties = FALSE) %>%
    pull(Term) %>% 
    unique()
  
  data_vis <- 
    data_selected %>%
    filter(Term %in% top_terms) %>% 
    mutate(
      is_significant =
        Adjusted.P.value <= max_p_adj &
        Odds.Ratio >= min_odds_ratio &
        overlap_size >= min_overlap_size,
      Odds.Ratio = pmin(Odds.Ratio, 10^log_odds_cap),
      Term =
        as_factor(Term) %>%
        fct_reorder(str_c(cell_type, group), n_distinct)
    )
  
  p <- 
    ggplot(data_vis, aes(group, Term, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "Reds",
      direction = 1,
      limits = c(0, log_odds_cap)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cell_type), drop = FALSE, nrow = 1) +
    labs(
      y = "",
      color = TeX("log_{10} (odds ratio)"),
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "Enrichr results ({db})"),
      caption = str_glue(
        "top {top_n} terms per cell type; ",
        "bordered circles: adjusted p value <= {max_p_adj}, ",
        "odds ratio >= {min_odds_ratio}, ",
        "overlap size >= {min_overlap_size}; ",
        "color scale capped at {log_odds_cap}"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  if (filename == "auto") {
    filename <- str_glue("dge_mm/enrichr_{db}")
  }
  ggsave_default(filename, ...)
  p
}

plot_enrichr_dots(dge$enrichr,
                  db = "MSigDB_Hallmark_2020",
                  width = 420)


## GSEA results ----

dge$gsea %>%
  ggplot(aes(NES, -log10(padj))) +
  geom_point(alpha = .25)


#' Generate a dotplot for terms enriched by GSEA.
#' 
#' Select the top "important" terms (i.e., terms that meet given limits of
#' adjusted p-value and normalized enrichment score, NES) per contrast and
#' cluster. Plot dots for these terms in all conditions, border the important
#' ones.
#'
#' @param data Results from `perform_gsea()`.
#' @param db Enrichr database for which results should be plotted.
#' @param top_n_positive Only plot n terms with the highest NES per cluster.
#' @param top_n_negative Only plot n terms with the lowest NES per cluster.
#' @param max_p_adj Maximum adjusted p value and …
#' @param min_NES minimum absolute NES required for bordered dots.
#' @param filename Name of output file. If `"auto"`, derive from database name.
#' @param ... Additional parameters passed to `ggsave_default()`.
#'
#' @return A ggplot object.
plot_gsea_dots <- function(data,
                           db,
                           top_n_positive = 5L,
                           top_n_negative = 5L,
                           max_p_adj = 0.05,
                           min_abs_NES = 1,
                           filename = "auto",
                           ...) {
  data_top_terms <-
    data %>% 
    filter(db == {{db}}, padj <= max_p_adj, abs(NES) >= min_abs_NES) %>% 
    group_by(group, cell_type)
  
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
    filter(db == {{db}}, pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      is_significant =
        padj <= max_p_adj &
        abs(NES) >= min_abs_NES,
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE)
    )
  
  if (nlevels(data_vis$pathway) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(5, nlevels(data_vis$pathway), 5),
        size = 0.5,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  color_limit <- max(abs(data_vis$NES))
  
  p <- 
    ggplot(data_vis, aes(group, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "PiYG",
      direction = 1,
      limits = c(-color_limit, color_limit)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    labs(
      y = "",
      color = "normalized\nenrichment\nscore",
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue("GSEA results ({db})"),
      caption = str_glue(
        "top {top_n_positive} positively and top {top_n_negative} negatively ",
        "enriched terms per cell_type; ",
        "bordered circles: adjusted p value <= {max_p_adj}, ",
        "absolute NES >= {min_abs_NES}"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  if (filename == "auto") {
    filename <- str_glue("dge_mm/gsea_{db}")
  }
  ggsave_default(filename, width = 400, ...)
  p
}

plot_gsea_dots(dge$gsea,
               db = "MSigDB_Hallmark_2020")



## DGE in tumor cluster ----

dge$results_tumor %>% 
  ggplot(aes(logFC_II_vs_IV, -log10(p_II_vs_IV))) +
  geom_point() +
  facet_wrap(vars(subclusters)) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("dge_mm/volcano_tumor")


dge$results_tumor_wide %>% 
  filter(subclusters == "all", logFC < 10) %>%
  inner_join(
    read_csv("metadata/rifatbegovic2018_table_s5.csv", comment = "#"),
    by = "gene"
  ) %>%
  rename(logfc_sc = logFC, q_sc = p_adj, logfc_bulk = logfc, q_bulk = q) %>%
  ggplot(aes(logfc_sc, logfc_bulk)) +
  geom_point(aes(size = -log10(q_sc)), alpha = .5) +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = gene), seed = 42) +
  scale_radius(range = c(0.5, 6)) +
  coord_fixed() +
  labs(
    x = "log FC (scRNA-seq)",
    y = "log FC (bulk)",
    size = TeX("-log_{10} (p_{adj})"),
    caption = "bulk values: Table S5 from Rifatbegovic et al 2018"
  ) +
  theme_bw()

ggsave_default("dge_mm/mycn_bulk_vs_sc")


plot_lfcc_heatmap_subclusters <- function() {
  corr_mat <- 
    dge$results_tumor_wide %>% 
    select(gene, logFC, subclusters) %>%
    filter(abs(logFC) <= 8) %>% 
    pivot_wider(
      names_from = subclusters,
      names_prefix = "NB_",
      values_from = c(logFC)
    ) %>% 
    select(NB_1:NB_4) %>% 
    cor(use = "pairwise.complete.obs")
  
  distance <- as.dist(1 - corr_mat)
  
  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(0, max(corr_mat[lower.tri(corr_mat)]), length.out = 9),
      scico(9, palette = "davos", direction = -1),
    ),
    name = "log fold change\ncorrelation",
    
    clustering_distance_rows = distance,
    clustering_distance_columns = distance,
  )
}

plot_lfcc_heatmap_subclusters()



# Publication figures -----------------------------------------------------

## Figure 3c ----

make_matrix <- function(gene,
                        cell_type,
                        groups,
                        direction = c("up", "down"),
                        row = 1,
                        col = 1) {
  direction <- match.arg(direction)
  
  barcodes <- 
    dge$metadata %>% 
    filter(
      cellont_cluster %in% dge$used_clusters,
      cellont_abbr == {{cell_type}}
    ) %>% 
    pull(cell)
  
  logcounts(dge$cds)[gene, barcodes, drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    magrittr::set_colnames("logexp") %>%
    as_tibble(rownames = "cell") %>%
    left_join(dge$metadata, by = "cell") %>%
    filter(group %in% {{groups}}) %>%
    transmute(
      row = row,
      col = col,
      label = str_glue("{cell_type}, {gene} ({direction})"),
      logexp = logexp / max(logexp),
      sample = rename_patients(sample),
      group = rename_groups(group)
    )
}


plot_violin <- function(genes) {
  plot_data <- pmap_dfr(genes, make_matrix)
  
  ggplot(plot_data, aes(sample, logexp)) +
    geom_violin(
      aes(fill = group),
      size = BASE_LINE_SIZE,
      # color = "gray60",
      scale = "width",
      width = 0.8,
      show.legend = FALSE
    ) +
    stat_summary(geom = "point", fun = mean, size = .2) +
    geom_text_npc(
      data = distinct(plot_data, row, col, label),
      aes(label = label),
      npcx = 0.05,
      npcy = 0.95,
      size = BASE_TEXT_SIZE_MM,
      hjust = 0
    ) +
    xlab("patient") +
    scale_y_continuous(
      "log-normalized expression",
      limits = c(0, 1.1)
    ) +
    scale_fill_manual(values = GROUP_COLORS) +
    facet_grid(vars(row), vars(col), scales = "free_x", space = "free_x") +
    theme_nb(grid = FALSE) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank()
    )
}

tribble(
  ~row, ~col, ~gene, ~cell_type, ~groups, ~direction,
  1,    1,    "IRF9",     "B",  c("I", "II"),  "up",
  1,    2,    "WDR74",    "B",  c("I", "III"), "up",
  1,    3,    "IFI44L",   "SC", c("I", "IV"),  "up",
  2,    1,    "SAP30",    "SC", c("I", "II"),  "down",
  2,    2,    "IRF9",     "B",  c("I", "III"), "up",
  2,    3,    "HIST1H1E", "SC", c("I", "IV"),  "down"
) %>% 
  plot_violin()
ggsave_publication("3c_exp_violin", width = 10, height = 6)



## Figure 3d ----

plot_gsea <- function(data,
                      db,
                      circle_significant = FALSE,
                      top_n_positive = 5L,
                      top_n_negative = 5L,
                      max_p_adj = 0.05,
                      min_abs_NES = 1) {
  data_top_terms <-
    data %>% 
    filter(db == {{db}}, padj <= max_p_adj, abs(NES) >= min_abs_NES) %>% 
    group_by(group, cell_type)
  
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
    filter(
      group != "tif",
      db == {{db}},
      pathway %in% c(top_terms_pos, top_terms_neg)
    ) %>% 
    mutate(
      is_significant =
        padj <= max_p_adj &
        abs(NES) >= min_abs_NES,
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE),
      group = factor(group) %>% rename_groups(),
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
  
  color_limit <- max(abs(data_vis$NES))
  
  if (circle_significant) {
    circle_significant <- geom_point(
      data = data_vis %>% filter(is_significant),
      shape = 1
    )
  } else {
    circle_significant <- NULL
  }
  
  p <- 
    ggplot(data_vis, aes(group, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    circle_significant +
    xlab("vs C (contrast)") +
    ylab(NULL) +
    scale_color_gsea(
      "normalized enrichment score",
      limits = c(-color_limit, color_limit),
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        label.position = "top",
        title.vjust = 0.1
      )
    ) +
    scale_size_area(
      name = TeX("-log_{10} p_{adj}"),
      max_size = 2.5
    )  +
    coord_fixed() +
    facet_wrap(vars(cell_type), nrow = 1) +
    theme_nb(grid = FALSE) +
    theme(
      legend.box.just = "bottom",
      legend.box.margin = margin(0, 0, 0, -25, "mm"),
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = "bottom",
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(-3, 1, 0, 1, "mm"),
      panel.spacing = unit(-.5, "pt"),
    )
  
  p
}
dge$gsea
plot_gsea(dge$gsea, "MSigDB_Hallmark_2020")
ggsave_publication("3d_gsea", width = 8, height = 10)



## Figure 3e ----

plot_logfc_correlation_heatmap <- function() {
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "pt"),
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(1, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(1, "mm")
  )
  
  corr_mat <- 
    dge$results_wide %>% 
    filter(abs(logFC) < 8) %>%
    select(gene, cell_type, group, logFC) %>%
    mutate(group = rename_groups(group)) %>%
    unite(cell_type, group, col = "cluster_group", sep = "_") %>%
    pivot_wider(names_from = "cluster_group", values_from = "logFC") %>%
    select(!gene) %>%
    cor(use = "pairwise.complete.obs")
  
  distance <- as.dist(1 - corr_mat)
  
  metadata <- 
    colnames(corr_mat) %>% 
    str_match("(.+)_(.+)")
  
  group_names <- 
    metadata[, 3] %>% 
    factor(levels = names(GROUP_COLORS)) %>% 
    fct_drop()
  
  colnames(corr_mat) <- metadata[, 2]
  rownames(corr_mat) <- metadata[, 2]
  
  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(0, max(corr_mat[lower.tri(corr_mat)]), length.out = 9),
      scico(9, palette = "davos", direction = -1),
    ),
    heatmap_legend_param = list(
      at = c(0, .4, .8),
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
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_dend_gp = gpar(lwd = 0.5),
    column_dend_height = unit(3, "mm"),
    
    height = unit(3.5, "cm"),
    width = unit(3.5, "cm"),
    border = F,
    
    right_annotation = rowAnnotation(
      group = group_names,
      col = list(group = GROUP_COLORS),
      show_annotation_name = FALSE,
      show_legend = TRUE,
      annotation_legend_param = list(
        group = list(
          title = "vs C\n(contrast)",
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    ),
    
    bottom_annotation = HeatmapAnnotation(
      group = group_names,
      col = list(group = GROUP_COLORS),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
  )
}

p <- plot_logfc_correlation_heatmap()
p
ggsave_publication("3e_logfc_correlation_heatmap",
                   plot = p, width = 6, height = 5)



## Figure S2x ----

plot_bulk_sc_comparison <- function() {
  labeled_genes <- tribble(
    ~gene, ~nudge_x, ~nudge_y,
    "MYCN",   -1,  1,
    "MYCNOS",  0,  2,
    "VIP",     0,  1,
    "NPW",     2,  2,
    "NKAIN4",  0, -3,
    
    "CXorf57", 0,  3,
    "IL7",     1,  1, 
    "MUC15",   0, -3,
    "KRT19",   0, -1
  )
  
  data <- 
    dge$results_tumor_wide %>% 
    filter(subclusters == "all", logFC < 10) %>%
    inner_join(
      read_csv("metadata/rifatbegovic2018_table_s5.csv", comment = "#"),
      by = "gene"
    ) %>%
    rename(
      logfc_sc = logFC,
      q_sc = p_adj,
      logfc_bulk = logfc,
      q_bulk = q
    )
  
  label_data <- 
    data %>%
    filter(gene %in% labeled_genes$gene) %>% 
    mutate(gene = factor(gene, levels = labeled_genes$gene)) %>% 
    arrange(gene)
  
  ggplot(data, aes(logfc_sc, logfc_bulk)) +
    geom_point(
      aes(size = -log10(q_sc)),
      color = "#e7298a",
      shape = 16, 
      alpha = .5
    ) +
    geom_smooth(
      method = "lm",
      fullrange = TRUE,
      color = "#7570b3",
      fill = "#7570b3",
      size = BASE_LINE_SIZE
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = gene),
      nudge_x = labeled_genes$nudge_x,
      nudge_y = labeled_genes$nudge_y,
      seed = 42, 
      size = BASE_TEXT_SIZE_MM,
      segment.size = BASE_LINE_SIZE,
      min.segment.length = unit(0, "mm")
    ) +
    scale_x_continuous(
      "log fold change (scRNA-seq)",
      limits = c(-10, 10),
      expand = expansion(add = 0.1)
    ) +
    scale_y_continuous(
      "log fold change (bulk RNA-seq)",
      limits = c(-10, 10),
      expand = expansion(add = 0.1)
    ) +
    scale_radius(
      TeX("-log_{10} p_{adj}"),
      range = c(0.1, 2)
    ) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = c(.8, .2)
    )
}

plot_bulk_sc_comparison()
ggsave_publication("2x_comparison_bulk_sc", width = 4, height = 4)



# Publication tables ------------------------------------------------------

## Table S3 ----

dge$results_wide_filtered %>% 
  arrange(group, cell_type, desc(logFC)) %>%
  mutate(
    group = rename_contrast_long(group),
    cell_type = CELL_TYPE_ABBREVIATIONS[cell_type]
  ) %>%
  select(
    "Contrast" = group,
    "Cell type" = cell_type,
    "Gene" = gene,
    "Log fold change" = logFC,
    "Adjusted p-value" = p_adj
  ) %>%
  save_table("S3_dge", "Microenvironment")


## Table S4 ----

dge$gsea %>% 
  arrange(db, group, cell_type, desc(NES)) %>%
  mutate(
    group = rename_contrast_long(group),
    cell_type = CELL_TYPE_ABBREVIATIONS[cell_type],
    leadingEdge = map_chr(leadingEdge, str_c, collapse = ", ")
  ) %>%
  select(
    "Database" = db,
    "Contrast" = group,
    "Cell type" = cell_type,
    "Pathway" = pathway,
    "Normalized Enrichment Score" = NES,
    "Adjusted p-value" = padj,
    "Leading edge genes" = leadingEdge
  ) %>%
  save_table("S4_gsea", "Microenvironment")
