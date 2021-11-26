# Plot mixed model DGE results.
#
# @DEPI dge_mm_results.rds
# @DEPI dge_pb_results.rds

library(scater)
library(monocle3)
library(CellChat)
library(tidyverse)
library(latex2exp)
library(ComplexHeatmap)
library(muscat)
library(ggrepel)
library(scico)
library(ggpmisc)
library(patchwork)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

dge <- readRDS("data_generated/dge_mm_results.rds")
dge_pb <- readRDS("data_generated/dge_pb_results.rds")

# are there problematic genes with convergence <= -20 ?
any(dge$results_vs_C$convergence <= -20)
any(dge$results_vs_S$convergence <= -20)
any(dge$results_vs_A$convergence <= -20)
any(dge$results_MNA_vs_other$convergence <= -20)



# Number of DE genes ------------------------------------------------------

bind_rows(
  "p_adj <= 0.05, log2 FC > 1" = 
    dge$results_wide_filtered %>% 
    filter(logFC > log(2), p_adj <= 0.05) %>%
    count(cell_type, comparison),
  "p_adj <= 0.05" = 
    dge$results_wide_filtered %>% 
    filter(p_adj <= 0.05) %>%
    count(cell_type, comparison),
  .id = "filter"
) %>% 
  complete(filter, cell_type, comparison, fill = list(n = 0L)) %>% 
  mutate(
    comparison =
      rename_contrast(comparison) %>%
      factor(levels = names(CONTRAST_COLORS))
  ) %>%
  ggplot(aes(cell_type, n)) +
  geom_col(aes(fill = comparison), position = "dodge") +
  xlab(NULL) +
  ylab("number of DE genes") +
  ylim(0, 2300) +
  scale_fill_manual(
    NULL,
    values = CONTRAST_COLORS,
    guide = guide_legend(nrow = 1)
  ) +
  facet_wrap(vars(filter)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("dge_mm/number_of_de_genes", height = 80, width = 150)


dge$results_wide_filtered %>%
  mutate(comparison = rename_contrast(comparison)) %>% 
  filter(abs(logFC) > log(2), comparison %>% str_ends("c")) %>%
  group_by(cell_type, direction, comparison) %>% 
  summarise(genes = list(unique(gene))) %>% 
  summarise(shared = list(genes)) %>%
  ungroup() %>% 
  rowwise() %>% 
  mutate(shared = shared %>% reduce(intersect) %>% length()) %>%
  mutate(shared = if_else(direction == "up", shared, -shared)) %>% 
  ggplot(aes(cell_type, shared, fill = direction)) +
  geom_col() +
  theme_bw()
ggsave_default("dge_mm/nu ber_of_consistent_genes")


  
# Volcano plots -----------------------------------------------------------

ggplot(dge$results_wide, aes(logFC, -log10(p))) +
  geom_point(size = 0.1) +
  facet_grid(vars(comparison), vars(cell_type)) +
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



# Violin plots ------------------------------------------------------------

plot_violin <- function(gene, cell_type, group, ref_group = "I") {
  plotExpression(
    dge$cds[, dge$cds$cellont_abbr == cell_type &
              dge$cds$group %in% c(group, ref_group)],
    gene,
    x = "sample",
    colour_by = "group"
  )  
}

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

# ... cases that still have extreme log fold changes
plot_violin("HIST1H1B", "SC", "III")
plot_violin("PTPN6", "NK", "III")



# currently, the commented code below does not work!
# plot_top_violins <- function(cell_type, n = 10, direction = c("up", "down")) {
#   direction <- match.arg(direction)
#   
#   dge$results_wide_filtered %>% 
#     filter(
#       p_adj <= 0.05,
#       cell_type == "B",
#       direction == "up"
#     ) %>% 
#     arrange(desc(abs(logFC))) %>% 
#     mutate(col = factor(comparison) %>% as.integer()) %>%
#     extract(
#       comparison,
#       into = c("group", "ref_group"),
#       regex = "(.+)_vs_(.+)"
#     ) %>%
#     group_by(col) %>%
#     mutate(row = row_number()) %>%
#     ungroup() %>%
#     arrange(col, row) %>% 
#     select(row, col, gene, cell_type, ref_group, group, direction) %>%
#     filter(row <= n) %>% 
#     # plot_violin() %>% 
#     {.}
# }
# 
# plot_top_violins("B")
# ggsave_default("dge_mm/violins_B", height = 120, width = 150)
# plot_top_violins("M")
# ggsave_default("dge_mm/violins_M", height = 120, width = 150)



# Comparison to pseudobulk ------------------------------------------------

plot_comparison <- function(data, lim = NULL, filename = NULL) {
  # note that nebula returns natural log fold changes
  # -> add reference line with slope log(2)
  p <-
    dge_pb$results %>% 
    filter(contrast != "tif") %>% 
    select(gene, cell_type = cluster_id,
           contrast, logFC_pb = logFC, p_pb = p_val) %>% 
    extract(contrast, into = "group", regex = "(.+)_vs_I") %>% 
    left_join(
      data %>%
        rename(logFC_mm = logFC, p_mm = p) %>% 
        filter(comparison %in% c("II_vs_I", "III_vs_I", "IV_vs_I")) %>% 
        extract(comparison, "group", "(.+)_vs"),
      by = c("gene", "cell_type", "group")
    ) %>% 
    ggplot(aes(logFC_pb, logFC_mm)) +
    geom_point(size = 0.1) +
    geom_abline(intercept = 0, slope = log(2), color = "blue", alpha = .25) +
    facet_grid(vars(group), vars(cell_type)) +
    coord_fixed(xlim = lim, ylim = lim) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  ggsave_default(filename)
  p
}

plot_comparison(
  dge$results_wide,
  lim = NULL,
  filename = "dge_mm/comparison_pb_mm_full"
)
plot_comparison(
  dge$results_wide_filtered,
  lim = NULL,
  filename = "dge_mm/comparison_pb_mm_full_filtered"
)
plot_comparison(
  dge$results_wide,
  lim = c(-10, 10),
  filename = "dge_mm/comparison_pb_mm"
)



# Enrichr results ---------------------------------------------------------

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
    group_by(comparison, cell_type) %>%
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
        fct_reorder(str_c(cell_type, comparison), n_distinct)
    )
  
  p <- 
    ggplot(data_vis, aes(comparison, Term, size = -log10(Adjusted.P.value))) +
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



# GSEA results ------------------------------------------------------------

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
    filter(db == {{db}}, pathway %in% c(top_terms_pos, top_terms_neg)) %>% 
    mutate(
      is_significant =
        padj <= max_p_adj &
        abs(NES) >= min_abs_NES,
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE),
      cell_type = factor(cell_type, names(CELL_TYPE_ABBREVIATIONS))
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
    ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
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



# Pathway genes -----------------------------------------------------------

pathways <- c(
  "TNF-alpha Signaling via NF-kB",
  "Interferon Gamma Response",
  "Interferon Alpha Response",
  "Inflammatory Response",
  "IL-2/STAT5 Signaling",
  "p53 Pathway"
)


## Single-cell heatmap ----

plot_pathway_heatmap <- function(db,
                                 pathways,
                                 scale_by_cell_type = TRUE) {
  sig_genes <- 
    dge$results_wide_filtered %>% 
    filter(logFC > 1, p_adj <= 0.05) %>% 
    pull(gene)
  
  row_metadata <-
    dge$gene_sets[[db]][pathways] %>%
    enframe("pathway", "gene") %>%
    unnest_longer(gene) %>%
    mutate(pathway = as_factor(pathway)) %>%
    filter(
      gene %in% rownames(dge$cds),
      gene %in% sig_genes
    )
  
  set.seed(1)
  col_metadata <-
    dge$metadata %>% 
    filter(
      cellont_abbr %in% c("T", "NK", "B", "M"),
      cell %in% colnames(dge$cds)
    ) %>% 
    group_by(cellont_abbr, group, sample) %>% 
    slice_sample(n = 150) %>% 
    group_by(cellont_abbr, group) %>%
    slice_sample(prop = 1) %>%
    mutate(
      cell_type = factor(cellont_abbr, names(CELL_TYPE_ABBREVIATIONS)),
      group = rename_groups(group),
      sample = rename_patients(sample) %>% factor(levels = PATIENT_ORDER)
    ) %>% 
    arrange(group)
  
  mat <-
    dge$cds %>% 
    logcounts() %>% 
    magrittr::extract(row_metadata$gene, col_metadata$cell) %>% 
    as.matrix()
  
  if (scale_by_cell_type) {
    mat <-
      map(
        unique(col_metadata$cell_type),
        ~mat[, col_metadata$cell_type == .]
        %>% t() %>% scale() %>% t() %>%
          replace_na(min(., na.rm = TRUE))
      ) %>%
      reduce(cbind)
  } else {
    mat <- 
      mat %>% 
      t() %>% scale() %>% t() %>%
      replace_na(min(., na.rm = TRUE))
  }
  
  
  Heatmap(
    mat,
    name = "scaled\nlogcounts",
    col = circlize::colorRamp2(
      seq(
        # quantile(mat, 0.1, na.rm = TRUE),
        # quantile(mat, 0.9, na.rm = TRUE),
        0,
        1.5,
        length.out = 9
      ),
      scico(9, palette = "oslo", direction = -1),
    ),
    border = TRUE,
    
    show_row_names = FALSE,
    row_split = row_metadata$pathway,
    row_title_rot = 0,
    cluster_row_slices = FALSE,
    show_row_dend = TRUE,
    row_gap = unit(0, "mm"),
    
    cluster_columns = FALSE,
    column_split = col_metadata$cell_type,
    show_column_names = FALSE,
    column_gap = unit(0, "mm"),
    
    top_annotation = HeatmapAnnotation(
      group = col_metadata$group,
      col = list(
        group = GROUP_COLORS
      )
    )
  )
}

(p <- plot_pathway_heatmap("MSigDB_Hallmark_2020", pathways, FALSE))
ggsave_default("dge_mm/pathway_heatmap", plot = p)



## Dotplot ----

plot_pathway_dots <- function(db, pathway,
                              min_exp = -2.5, max_exp = 2.5) {
  scale_and_limit <- function(x) {
    scale(x)[,1] %>% 
      pmax(min_exp) %>% 
      pmin(max_exp)
  }
  
  # assemble row and col metadata
  col_metadata <-
    dge$metadata %>% 
    filter(
      cellont_abbr %in% c("T", "NK", "B", "M"),
      cell %in% colnames(dge$cds)
    ) %>%
    mutate(group_type = fct_cross(group, cellont_abbr))
    # mutate(group_type = fct_cross(cellont_abbr, group))
  
  panel_annotation <- tribble(
    ~xmin, ~xmax, ~fill,
    4.5, 8.5, "grey90",
    12.5, 16.5, "grey90"
  )
  counts <- logcounts(dge$cds[, col_metadata$cell])
  features <- intersect(dge$gene_sets[[db]][[pathway]], rownames(counts))
  groups <- col_metadata$group_type
  
  
  # basically code of plot_dots, with counts, features, groups,
  # and panel_annotation prefilled and features sorted by percent expressed
  known_features <- intersect(features, rownames(counts))
  missing_features <- setdiff(features, rownames(counts))
  if (length(missing_features) > 0)
    warn(
      "The following requested features are missing: ",
      "{str_c(missing_features, collapse = ', ')}"
    )
  
  vis_data <- 
    counts[known_features, , drop = FALSE] %>% 
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
    # new from here
    mutate(
      feature = factor(feature, levels = features) %>% fct_reorder(pct_exp)
    )
  
  if (!is.null(panel_annotation)) {
    panel_bg <- list(
      geom_point(color = "white"),  # initialize discrete coordinate system
      geom_rect(
        data = panel_annotation,
        aes(xmin = xmin, xmax = xmax,
            ymin = 0.5, ymax = nlevels(vis_data$feature) + 0.5,
            fill = fill),
        show.legend = FALSE,
        inherit.aes = FALSE,
      ),
      scale_fill_identity()
    )
  } else {
    panel_bg <- NULL
  }
  
  ggplot(vis_data, aes(id, feature)) +
    panel_bg +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_x_discrete("cluster", expand = expansion(add = 0.5)) +
    scale_y_discrete("feature", expand = expansion(add = 0.5)) +
    scale_color_scico(
      "scaled\naverage\nexpression",
      palette = "oslo",
      direction = -1,
      aesthetics = "color"
    ) +
    scale_radius("% expressed", range = c(0, 6)) +
    coord_fixed(
      # xlim = c(0.5, nlevels(vis_data$id) + 0.5),
      # ylim = c(0.5, nlevels(vis_data$feature) + 0.5),
      clip = "off"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      panel.grid = element_blank()
    )
}

plot_pathway_dots("MSigDB_Hallmark_2020", pathways[1])
ggsave_default("dge_mm/pathway_dots_tnf", height = 900)

plot_pathway_dots("MSigDB_Hallmark_2020", pathways[2])
ggsave_default("dge_mm/pathway_dots_ifng", height = 900)

plot_pathway_dots("MSigDB_Hallmark_2020", pathways[3])
ggsave_default("dge_mm/pathway_dots_ifna", height = 600)

plot_pathway_dots("MSigDB_Hallmark_2020", pathways[4])
ggsave_default("dge_mm/pathway_dots_inflammatory", height = 900)

plot_pathway_dots("MSigDB_Hallmark_2020", pathways[5])
ggsave_default("dge_mm/pathway_dots_il2stat5", height = 900)

plot_pathway_dots("MSigDB_Hallmark_2020", pathways[6])
ggsave_default("dge_mm/pathway_dots_p53", height = 900)



## logFC ----

plot_pathway_genes <- function(db, pathway, filename = "auto", ...) {
  vis_data <- 
    dge$results_wide %>% 
    filter(
      gene %in% dge$gene_sets[[db]][[pathway]],
      abs(logFC) < 8,
      comparison %>% str_ends("_I")
    ) %>% 
    mutate(comparison = rename_contrast(comparison))
  
  p <-
    vis_data %>% 
    mutate(gene = factor(gene) %>% fct_reorder(logFC)) %>% 
    ggplot(aes(comparison, gene, size = -log10(p_adj), color = logFC)) +
    geom_point() +
    scale_color_gsea(limits = c(-3, 3)) +
    scale_size_area() +
    facet_wrap(vars(cell_type), nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = .5),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  if (filename == "auto") {
    pathway_name <- str_replace_all(pathway, "/", "_")
    filename <- str_glue("dge_mm/pathway_genes_{db}_{pathway_name}")
  }
  ggsave_default(filename, width = 120, ...)
  p
}

plot_pathway_genes("MSigDB_Hallmark_2020", pathways[1], height = 700)
plot_pathway_genes("MSigDB_Hallmark_2020", pathways[2], height = 600)
plot_pathway_genes("MSigDB_Hallmark_2020", pathways[3], height = 300)
plot_pathway_genes("MSigDB_Hallmark_2020", pathways[4], height = 600)
plot_pathway_genes("MSigDB_Hallmark_2020", pathways[5], height = 500)
plot_pathway_genes("MSigDB_Hallmark_2020", pathways[6], height = 600)




# Tumor cells -------------------------------------------------------------

## Comparison to bulk data ----

dge$results_wide %>% 
  filter(cell_type == "NB", comparison == "II_vs_IV", logFC < 10) %>% 
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



## Pseudobulk correlation ----

plot_expc_heatmap_samples <- function() {
  cds_tumor <- dge$cds[, colData(dge$cds)$cellont_abbr == "NB"]
  
  pb_tumor <- aggregateData(
    cds_tumor,
    assay = "logcounts",
    fun = "mean",
    by = "sample"
  )
  colnames(pb_tumor) <-
    colnames(pb_tumor) %>%
    rename_patients()
  
  hvgs <-
    cds_tumor %>%
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
  
  Heatmap(
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
    
    show_column_dend = FALSE,
    
    width = unit(60, "mm"),
    height = unit(60, "mm"),
    
    left_annotation = rowAnnotation(
      group = group_names,
      col = list(group = GROUP_COLORS),
      show_annotation_name = FALSE,
      show_legend = TRUE,
      annotation_legend_param = list(
        group = list(
          title = "group"
        )
      )
    )
  )
}

(p <- plot_expc_heatmap_samples())
ggsave_default(
  "dge_mm/pseudobulk_correlation_tumor",
  plot = p
)



# Pseudobulk correlation of all cells -------------------------------------

plot_pbc_heatmap <- function(level = "group") {
  pb <- aggregateData(
    dge$cds,
    assay = "logcounts",
    fun = "mean",
    by = c("cellont_abbr", level)
  )
  
  if (level == "group") {
    colnames(pb) <-
      colnames(pb) %>%
      rename_groups()
  } else {
    colnames(pb) <-
      colnames(pb) %>%
      rename_patients()
  }
  
  hvgs <-
    dge$cds %>%
    scran::modelGeneVar() %>% 
    scran::getTopHVGs()
  
  mat <-
    assays(pb) %>% 
    as.list() %>% 
    imap(~magrittr::set_colnames(., str_c(colnames(.x), .y, sep="_"))) %>% 
    reduce(cbind) %>% 
    magrittr::extract(hvgs, )
  
  missing_cells <- colSums(mat) %>% near(0)
  
  corr_mat <-
    mat[, !missing_cells] %>% 
    cor(use = "pairwise.complete.obs")
  
  distance <- as.dist(1 - corr_mat)
  
  col_metadata <-
    colnames(corr_mat) %>%
    str_match("(.+)_(.+)") %>% 
    magrittr::set_colnames(c("label", "sample", "cell_type")) %>% 
    as_tibble() %>% 
    mutate(group = str_sub(sample, 1, 1))
  
  if (level == "group") {
    left_annotation <- rowAnnotation(
      group = col_metadata$group,
      cell_type = col_metadata$cell_type,
      col = list(
        group = GROUP_COLORS,
        cell_type = CELL_TYPE_COLORS
      )
    )
  } else {
    left_annotation <- rowAnnotation(
      sample = col_metadata$sample,
      group = col_metadata$group,
      cell_type = col_metadata$cell_type,
      col = list(
        sample = PATIENT_COLORS,
        group = GROUP_COLORS,
        cell_type = CELL_TYPE_COLORS
      )
    )
  }
  
  Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(min(corr_mat), 1, length.out = 9),
      scico(9, palette = "davos", direction = -1),
    ),
    name = "correlation of\npseudobulk\nexpression",
    heatmap_legend_param = list(
      at = round(c(min(corr_mat), max(corr_mat)), 2)
    ),
    
    clustering_distance_rows = distance,
    clustering_distance_columns = distance,
    
    show_column_dend = FALSE,
    
    width = unit(200, "mm"),
    height = unit(200, "mm"),
    
    left_annotation = left_annotation
  )
}

(p <- plot_pbc_heatmap())

ggsave_default(
  "dge_mm/pseudobulk_correlation_group",
  plot = p,
  height = 700,
  width = 700
)

(p <- plot_pbc_heatmap("sample"))

ggsave_default(
  "dge_mm/pseudobulk_correlation_sample",
  plot = p,
  height = 700,
  width = 700
)



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
    filter(
      comparison != "tif",
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
      comparison = factor(comparison) %>% rename_contrast(),
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
    ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    circle_significant +
    # xlab("vs C (contrast)") +
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
      axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
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

plot_gsea(dge$gsea, "MSigDB_Hallmark_2020")
ggsave_publication("3d_gsea", width = 13, height = 10)  # was 8 x 10



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
      col = list(group = CONTRAST_COLORS),
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
      col = list(group = CONTRAST_COLORS),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
  )
}

(p <- plot_logfc_correlation_heatmap())
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
    dge$results_wide %>% 
    filter(cell_type == "NB", comparison == "II_vs_IV", logFC < 10) %>% 
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

# dge$results_wide_filtered %>% 
#   arrange(comparison, cell_type, desc(logFC)) %>%
#   mutate(
#     comparison = rename_contrast_long(comparison),
#     cell_type = CELL_TYPE_ABBREVIATIONS[cell_type]
#   ) %>%
#   select(
#     "Comparison" = comparison,
#     "Cell type" = cell_type,
#     "Gene" = gene,
#     "Log fold change" = logFC,
#     "Adjusted p-value" = p_adj
#   ) %>%
#   left_join(gene_pathways, by = "Gene") %>% 
#   save_table("S3_dge", "Microenvironment")

dge$results_wide_filtered %>% 
  arrange(comparison, cell_type, desc(logFC)) %>%
  mutate(
    logFC = logFC / log(2),  # nebula returns natural log fold changes
    comparison = rename_contrast(comparison),
    cell_type = CELL_TYPE_ABBREVIATIONS[cell_type]
  ) %>%
  select(!c(p, frq, frq_ref, direction)) %>% 
  pivot_wider(names_from = comparison, values_from = c(logFC, p_adj)) %>% 
  left_join(gene_pathways, by = c(gene = "Gene")) %>%
  split(.$cell_type) %>% 
  map(select, !cell_type) %>% 
  save_table("S3_dge")



## Table S4 ----

# dge$gsea %>% 
#   arrange(db, comparison, cell_type, desc(NES)) %>%
#   mutate(
#     comparison = rename_contrast_long(comparison),
#     cell_type = CELL_TYPE_ABBREVIATIONS[cell_type],
#     leadingEdge = map_chr(leadingEdge, str_c, collapse = ", ")
#   ) %>%
#   select(
#     "Database" = db,
#     "Comparison" = comparison,
#     "Cell type" = cell_type,
#     "Pathway" = pathway,
#     "Normalized Enrichment Score" = NES,
#     "Adjusted p-value" = padj,
#     "Leading edge genes" = leadingEdge
#   ) %>%
#   save_table("S4_gsea", "GSEA")

dge$gsea %>% 
  arrange(db, comparison, cell_type, desc(NES)) %>%
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
  save_table("S4_gsea")
