# Plot mixed model DGE results.
#
# @DEPI dge_mm_results.rds
# @DEPI dge_pb_results.rds

library(scater)
library(tidyverse)
library(latex2exp)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

dge <- readRDS("data_generated/dge_mm_results.rds")
dge_pb <- readRDS("data_generated/dge_pb_results.rds")

# are there problematic genes with convergence <= -20 ?
any(dge$results$convergence <= -20)



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
