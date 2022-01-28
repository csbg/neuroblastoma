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
stopifnot(
  all(dge$results_vs_C$convergence > -20),
  all(dge$results_MNA_vs_other$convergence > -20)
)



# Number of DE genes ------------------------------------------------------

## Counts ----

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


## Consistent genes ----

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
  )


ggplot(consistent_genes, aes(cell_type, shared_AMS)) +
  geom_col(fill = hcl(285, 100, 65)) +
  geom_hline(yintercept = 0, size = BASE_LINE_SIZE) +
  xlab("cell type") +
  ylab("number of consistently down- and upregulated genes") +
  theme_nb()
ggsave_default("dge_mm/number_of_consistent_genes", width = 40, height = 60)

consistent_genes %>% 
  select(!c(Ac, Mc, Sc)) %>% 
  pivot_longer(
    starts_with("shared"),
    names_to = "shared",
    names_prefix = "shared_",
    values_to = "n"
  ) %>% 
  mutate(shared = fct_relevel(shared, "AM", "AS", "MS")) %>% 
  ggplot(aes(cell_type, n, fill = shared)) +
  geom_col() +
  geom_hline(yintercept = 0, size = BASE_LINE_SIZE) +  
  xlab("cell type") +
  ylab("number of consistently down- and upregulated genes") +
  theme_nb()
ggsave_default("dge_mm/number_of_consistent_genes_groups", width = 60, height = 60)


  
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

selected_genes <- list(
  "TNF-alpha Signaling via NF-kB" = c("IL1B", "FOS", "NFKB1", "MYC", "IFNGR2",
                                      "CD44", "RELB", "SOD2", "MAP3K8"),
  "Interferon Gamma Response" = c("IFI44", "NFKB1", "MYD88", "HLA-A", "CD86",
                                  "HLA-DQA1", "HLA-DQB1", "TNFSF10"),
  "Myc Targets V1" = c("TRIM28", "PCNA", "CDK4", "MCM5", "HDAC2"),
  "E2F Targets" = c("DNMT1", "EZH2", "MKI67", "PCNA")
)


## Single-cell heatmap ----

plot_pathway_heatmap_sc <- function(db = "MSigDB_Hallmark_2020",
                                    pathways = c(
                                      "TNF-alpha Signaling via NF-kB",
                                      "Interferon Gamma Response",
                                      "Myc Targets V1",
                                      "E2F Targets"
                                    ),
                                    cell_types = c("B", "M"),
                                    genes = c("all", "selected"),
                                    norm_method = c("quantile", "scale")) {
  norm_method <- match.arg(norm_method)
  genes <- match.arg(genes)
  
  if (genes == "all") {
    sig_genes <- 
      dge$results_wide_filtered %>% 
      filter(
        # cell_type %in% {{cell_types}},
        comparison %in% c("II_vs_I", "III_vs_I", "IV_vs_I"),
        abs(logFC) > 1,
        p_adj <= 0.05
      ) %>% 
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
  } else {
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
  }
    
  set.seed(1)
  col_metadata <-
    dge$metadata %>% 
    filter(
      cellont_abbr %in% {{cell_types}},
      cell %in% colnames(dge$cds)
    ) %>% 
    group_by(cellont_abbr, group, sample) %>% 
    slice_sample(n = 100) %>% 
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
  
  if (norm_method == "quantile") {
    # cutoff at 95th percentile of logcounts
    mat <-
      map(
        unique(col_metadata$cell_type),
        function(cell_type) {
          m <- mat[, col_metadata$cell_type == cell_type]
          q <- apply(m, 1, quantile, .95)
          m %>%
            sweep(1, q, function(x, y) pmin(x / y, 1)) %>% 
            replace_na(min(., na.rm = TRUE))
        }
      ) %>%
      reduce(cbind)
    min_limit <- min(mat)
    max_limit <- max(mat)
  } else {
    # scaled logcounts
    mat <-
      mat %>%
      t() %>% scale() %>% t() %>%
      replace_na(min(., na.rm = TRUE))
    # mat <-
    #   map(
    #     unique(col_metadata$cell_type),
    #     ~mat[, col_metadata$cell_type == .]
    #     %>% t() %>% scale() %>% t() %>%
    #       replace_na(min(., na.rm = TRUE))
    #   ) %>%
    #   reduce(cbind)
    min_limit <- 0
    max_limit <- 1.5
  }
  
  Heatmap(
    mat,
    name = "expression",
    col = circlize::colorRamp2(
      seq(
        min_limit,
        max_limit,
        length.out = 9
      ),
      viridisLite::cividis(9),
      # scico(9, palette = "oslo", direction = -1)
    ),
    border = FALSE,
    heatmap_legend_param = list(
      at = c(min_limit, max_limit),
      labels = c("low", "high"),
      border = FALSE
    ),
    # use_raster = FALSE,
    
    show_row_names = FALSE,
    row_split = row_metadata$pathway,
    row_title_rot = 0,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    row_gap = unit(.5, "mm"),
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
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
      )
    )
  )
}

(p <- plot_pathway_heatmap_sc(norm_method = "scale"))
ggsave_default("dge_mm/pathway_heatmap_sc_scale_allX",
               plot = p, width = 180, height = 100)


(p <- plot_pathway_heatmap_sc(norm_method = "quantile"))
ggsave_default("dge_mm/pathway_heatmap_sc_quantile_allX", plot = p)






## Bulk heatmap ----

plot_pathway_heatmap <- function(db = "MSigDB_Hallmark_2020",
                                 pathways = c(
                                   "TNF-alpha Signaling via NF-kB" = "up",
                                   "Interferon Gamma Response" = "up",
                                   "Myc Targets V1" = "down",
                                   "E2F Targets" = "down"
                                 ),
                                 genes = c("top", "all", "selected"),
                                 level = c("sample", "group"),
                                 cell_types = c("B", "M")) {
  level <- match.arg(level)
  genes <- match.arg(genes)
  
  if (genes == "top") {
    sig_genes <-
      dge$results_wide_filtered %>% 
      filter(
        cell_type %in% {{cell_types}},
        comparison %in% c("II_vs_I", "III_vs_I", "IV_vs_I"),
        abs(logFC) > log(4),
        p_adj <= 0.05
      ) %>% 
      group_by(gene) %>% 
      summarise(mean_logFC = mean(logFC))
    
    row_metadata <- 
      dge$gene_sets[[db]][names(pathways)] %>%
      enframe("pathway", "gene") %>%
      unnest_longer(gene) %>%
      mutate(pathway = as_factor(pathway)) %>%
      filter(
        gene %in% rownames(dge$cds),
      ) %>% 
      left_join(sig_genes, by = "gene") %>%
      group_split(pathway) %>% 
      map2(
        pathways,
        ~if (.y == "up") {
          slice_max(.x, mean_logFC, n = 10)
        } else {
          slice_min(.x, mean_logFC, n = 10)
        }
      ) %>% 
      bind_rows()
  } else if (genes == "all") {
    sig_genes <-
      dge$results_wide_filtered %>%
      filter(
        cell_type %in% {{cell_types}},
        comparison %in% c("II_vs_I", "III_vs_I", "IV_vs_I"),
        abs(logFC) > log(4),
        p_adj <= 0.05
      ) %>%
      pull(gene)
    
    row_metadata <-
      dge$gene_sets[[db]][names(pathways)] %>%
      enframe("pathway", "gene") %>%
      unnest_longer(gene) %>%
      mutate(pathway = as_factor(pathway)) %>%
      filter(
        gene %in% rownames(dge$cds),
        gene %in% sig_genes
      )
  } else {
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
  }
  
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
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "pt"),
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(1, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(1, "mm")
  )
  
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

(p <- plot_pathway_heatmap(level = "sample", pathways = selected_genes,
                           genes = "selected"))
ggsave_default("dge_mm/pathway_heatmap_sample_selected",
               plot = p, width = 130, height = 60)

(p <- plot_pathway_heatmap(level = "group", pathways = selected_genes,
                           genes = "selected"))
ggsave_default("dge_mm/pathway_heatmap_group_selected",
               plot = p, width = 86, height = 60)

# (p <- plot_pathway_heatmap(level = "sample", genes = "all"))
# ggsave_default("dge_mm/pathway_heatmap_sample_all",
#                plot = p, width = 200, height = 200)
# 
# (p <- plot_pathway_heatmap(level = "group", genes = "all"))
# ggsave_default("dge_mm/pathway_heatmap_group_all",
#                plot = p, width = 150, height = 200)
# 
# 
# (p <- plot_pathway_heatmap(level = "sample", genes = "top"))
# ggsave_default("dge_mm/pathway_heatmap_sample_top",
#                plot = p, width = 200, height = 100)
# 
# (p <- plot_pathway_heatmap(level = "group", genes = "top"))
# ggsave_default("dge_mm/pathway_heatmap_group_top",
#                plot = p, width = 150, height = 100)




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
