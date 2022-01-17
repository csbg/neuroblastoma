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

dge <- readRDS("data_wip/dge_mm_myeloid_results.rds")

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
ggsave_default("dge_myeloid/number_of_de_genes", height = 80, width = 150)


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
ggsave_default("dge_myeloid/number_of_consistent_genes", width = 40, height = 60)

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
ggsave_default("dge_myeloid/number_of_consistent_genes_groups", width = 60, height = 60)



# Volcano plots -----------------------------------------------------------

ggplot(dge$results_wide, aes(logFC, -log10(p))) +
  geom_point(size = 0.1) +
  facet_grid(vars(comparison), vars(cell_type)) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 25)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("dge_myeloid/volcano")



# Violin plots ------------------------------------------------------------

# plot_violin <- function(gene, cell_type, group, ref_group = "I") {
#   plotExpression(
#     dge$cds[, dge$cds$cellont_abbr == cell_type &
#               dge$cds$group %in% c(group, ref_group)],
#     gene,
#     x = "sample",
#     colour_by = "group"
#   )  
# }



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
      cell_type = factor(cell_type)
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
      title = db
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
    filename <- str_glue("dge_myeloid/gsea_{db}")
  }
  ggsave_default(filename, width = 400, ...)
  p
}

dge$gsea %>% 
  filter(comparison != "Mas") %>% 
  mutate(
    cell_type = fct_recode(
      cell_type,
      "class" = "classical mono",
      "nonclass" = "nonclassical mono"
    )
  ) %>% 
  plot_gsea_dots(db = "MSigDB_Hallmark_2020", height = 150)

dge$gsea %>% 
  filter(comparison == "Mas") %>% 
  mutate(
    cell_type = fct_recode(
      cell_type,
      "class" = "classical mono",
      "nonclass" = "nonclassical mono"
    )
  ) %>% 
  plot_gsea_dots(db = "MSigDB_Hallmark_2020",
                 height = 150, filename = "dge_myeloid/gsea_Mas")




# LogFC correlation -------------------------------------------------------

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
    filter(abs(logFC) < 8, comparison %>% str_ends("c"), cell_type != "other") %>%
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
    
    height = unit(3, "cm"),
    width = unit(3, "cm"),
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
ggsave_default("dge_myeloid/logfc_correlation",
                   plot = p, width = 70, height = 70)
