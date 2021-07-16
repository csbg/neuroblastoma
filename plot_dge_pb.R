# Plot pseudo-bulk DGE results.
#
# @DEPI dge_pb_results.RData

library(monocle3)
library(scater)
library(muscat)
library(tidyverse)
library(viridis)
library(patchwork)
library(ggpmisc)
library(ComplexHeatmap)
library(latex2exp)
library(scico)
library(openxlsx)
library(ggrepel)
library(ggbeeswarm)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

load("data_generated/dge_pb_results.RData")



# Plot results ------------------------------------------------------------

## Volcano plot ----

#' Make one volcano plot per cluster for a given contrast.
#'
#' @param data A DGE dataset.
#' @param contrast A contrast present in the dataset.
#' @param max_p_adj see `filter_dge_results()`.
#' @param min_abs_log_fc see `filter_dge_results()`.
#' @param min_freq see `filter_dge_results()`.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_volcano <- function(data,
                         contrast,
                         max_p_adj = 0.05,
                         min_abs_log_fc = 1,
                         min_freq = 0.1,
                         filename = NULL) {
  contrast_groups <- str_match(contrast, "(.+)_vs_(.+)")[1,]
  group_left <- contrast_groups[2]
  group_right <- contrast_groups[3]
  
  data_filtered <- 
    data %>%
    filter(contrast == {{contrast}}) %>% 
    mutate(
      is_significant =
        p_adj.loc < max_p_adj &
        abs(logFC) > min_abs_log_fc,
      is_frequent =
        !!sym(str_glue("{group_left}.frq")) > min_freq |
        !!sym(str_glue("{group_right}.frq")) > min_freq,
      diff_expressed =
        case_when(
          is_significant & is_frequent & logFC > 0 ~ "up",  
          is_significant & is_frequent & logFC < 0 ~ "down",
          TRUE ~ "no change"
        ) %>% 
        factor(levels = c("up", "no change", "down"))
    )

  genes_count <- 
    data_filtered %>% 
    count(cluster_id, diff_expressed, .drop = FALSE) %>% 
    pivot_wider(names_from = "diff_expressed", values_from = "n")

  dot_colors <- c(
    `no change` = "gray80",
    up = "red",
    down = "blue"
  )

  p <-
    ggplot(data_filtered, aes(logFC, -log10(p_val))) +
    geom_point(
      aes(color = diff_expressed),
      alpha = .5,
      size = 0.5
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = down),
      npcx = 0.2,
      npcy = 0.9,
      size = 3,
      color = dot_colors["down"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = `no change`),
      npcx = 0.5,
      npcy = 0.9,
      size = 3,
      color = dot_colors["no change"],
      na.rm = TRUE
    ) +
    geom_text_npc(
      data = genes_count,
      aes(label = up),
      npcx = 0.8,
      npcy = 0.9,
      size = 3,
      color = dot_colors["up"],
      na.rm = TRUE
    ) +
    scale_color_manual(
      "differentially\nexpressed",
      values = dot_colors,
      guide = guide_legend(override.aes = list(alpha = 1, size = 3)),
      drop = FALSE
    ) +
    facet_wrap(vars(cluster_id), ncol = 6, drop = FALSE) +
    labs(
      title = str_glue("DGE using contrast {contrast}"),
      caption = str_glue(
        "gene expression frequency > {min_freq}, ",
        "adjusted p value < {max_p_adj}, ",
        "absolute log fold change > {min_abs_log_fc}"
      )
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  ggsave_default(filename)
  p
}

plot_volcano(dge_results, "II_vs_I", filename = "dge/volcano_II")
plot_volcano(dge_results, "III_vs_I", filename = "dge/volcano_III")
plot_volcano(dge_results, "IV_vs_I", filename = "dge/volcano_IV")



## Violin plot ----

#' Title
#'
#' @param data A DGE dataset.
#' @param contrast A contrast and ...
#' @param cluster ... cluster for which all samples should be plotted.
#' @param top_n Plot the top differentially expressed genes.
#' @param direction Plot up- or downregulated genes.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_violin <- function(data,
                        contrast,
                        cluster,
                        top_n = 10,
                        direction = c("up", "down"),
                        filename = NULL) {
  contrast_groups <- str_match(contrast, "(.+)_vs_(.+)")[1,]
  group_left <- contrast_groups[2]
  group_right <- contrast_groups[3]
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>%
    filter(
      contrast == {{contrast}},
      cluster_id == {{cluster}},
      direction == {{direction}}
    ) %>% 
    slice_max(n = top_n, logFC, with_ties = FALSE) %>% 
    pull(gene)
  
  p <-
    plotExpression(
      nb[, nb$cluster_id == cluster &
           nb$group_id %in% c(group_left, group_right)],
      top_genes,
      x = "sample_id",
      colour_by = "group_id",
      point_size = 1,
      ncol = NULL
    ) +
    labs(
      caption = str_glue(
        "{direction}regulated genes, cluster {cluster}, contrast {contrast}"
      )
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # calculate facet grid dimensions to ensure that all plots have the same size
  n_facet <- wrap_dims(length(top_genes))
  ggsave_default(
    filename,
    height = 30 + 60 * n_facet[1],
    width = 20 + 70 * n_facet[2]
  )
  p
}

# plot_violin(dge_results_filtered, "II_vs_I", cluster = "1",
#             filename = "dge/violin_II_vs_I_c1_up")

dge_results_filtered %>% 
  distinct(contrast, cluster = as.character(cluster_id), direction) %>% 
  pwalk(
    function(contrast, cluster, direction) {
      info("Plotting {direction}regulated genes ",
           "in cluster {cluster} in contrast {contrast}")
      plot_violin(
        dge_results_filtered,
        top_n = Inf,
        contrast = contrast,
        cluster = cluster,
        direction = direction,
        filename = str_glue("dge/violin_{contrast}_c{cluster}_{direction}")
      )
    }
  )



## Pseudo-bulk heatmap ----

#' Plot a pseudobulk heatmap.
#'
#' @param data A DGE dataset.
#' @param pb_data Pseudobulk data as generated by `muscat::aggregateData()`.
#' @param contrast The selected contrast.
#' @param top_n Plot the top differentially expressed genes.
#' @param clusters A character vector denoting clusters that should be included.
#' @param filename Name of output file.
#'
#' @return An object of class Heatmap.
plot_pbheatmap <- function(data,
                           pb_data,
                           contrast,
                           top_n = Inf,
                           clusters = NULL,
                           filename = NULL) {
  contrast_groups <- str_match(contrast, "(.+)_vs_(.+)")[1,]
  group_left <- contrast_groups[2]
  group_right <- contrast_groups[3]
  
  top_genes <- 
    data %>% 
    filter(contrast == {{contrast}}) %>% 
    {
      if (!is.null(clusters))
        filter(., cluster_id %in% {{clusters}})
      else
        .
    } %>% 
    group_by(cluster_id) %>% 
    slice_max(logFC, n = top_n, with_ties = FALSE)
  
  sample_data <- 
    colData(pb_data) %>%
    as_tibble(rownames = "sample") %>%
    filter(group_id %in% c(group_left, group_right))
  
  mat <- 
    top_genes %>% 
    group_map(~assays(pb_data)[[.y$cluster_id]][.x$gene, sample_data$sample]) %>%
    reduce(rbind) %>% 
    t() %>% 
    scale() %>% 
    t()
  
  p <-
    Heatmap(
      mat,
      col = scico(9, palette = "roma", direction = -1),
      
      cluster_rows = FALSE,
      show_row_names = FALSE,
      row_split = top_genes$cluster_id,
      row_title_rot = 0,
      
      cluster_columns = FALSE,
      column_split = sample_data$group_id,
      column_title = "group",
      
      left_annotation = rowAnnotation(
        logFC = top_genes$logFC,
        col = list(
          logFC = circlize::colorRamp2(
            breaks = c(min(data$logFC), 0, max(data$logFC)),
            colors = c("blue", "white", "red")
          )
        )
      ),
      
      top_annotation = HeatmapAnnotation(
        foo = anno_block(
          gp = gpar(fill = 3:4),
          labels = unique(sample_data$group_id),
        )
      ),

      heatmap_legend_param = list(title = "scaled\nmean\nlogcounts")
    ) %>% 
    draw(
      row_title = "cluster",
      column_title = "Pseudobulk heatmap",
      column_title_gp = gpar(fontface = "bold", fontsize = 14)
    )
  
  ggsave_default(filename, plot = p, width = 150, height = 297)
  p
}

pb_log <- aggregateData(
  nb,
  assay = "logcounts",
  fun = "mean",
  by = c("cluster_id", "sample_id")
)

plot_pbheatmap(dge_results_filtered, pb_log,
               contrast = "II_vs_I", filename = "dge/heatmap_II")
plot_pbheatmap(dge_results_filtered, pb_log,
               contrast = "III_vs_I", filename = "dge/heatmap_III")
plot_pbheatmap(dge_results_filtered, pb_log,
               contrast = "IV_vs_I", filename = "dge/heatmap_IV")



## Enrichment dot plot ----

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
    group_by(contrast, cluster) %>%
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
        fct_reorder(str_c(cluster, contrast), n_distinct)
    )
  
  p <- 
    ggplot(data_vis, aes(contrast, Term, size = -log10(Adjusted.P.value))) +
    geom_point(aes(color = log10(Odds.Ratio))) +
    geom_point(data = data_vis %>% filter(is_significant), shape = 1) +
    scale_color_distiller(
      palette = "Reds",
      direction = 1,
      limits = c(0, log_odds_cap)
    ) +
    scale_size_area() +
    coord_fixed() +
    facet_wrap(vars(cluster), drop = FALSE, nrow = 1) +
    labs(
      y = "",
      color = TeX("log_{10} (odds ratio)"),
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue( "Enrichr results ({db})"),
      caption = str_glue(
        "top {top_n} terms per cluster; ",
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
    filename <- str_glue("dge/enrichr_{db}")
  }
  ggsave_default(filename, ...)
  p
}


plot_enrichr_dots(enrichr_results,
                  db = "GO_Biological_Process_2018",
                  height = 400,
                  width = 400)

plot_enrichr_dots(enrichr_results,
                  db = "GO_Cellular_Component_2018")

plot_enrichr_dots(enrichr_results,
                  db = "GO_Molecular_Function_2018",
                  width = 420)

plot_enrichr_dots(enrichr_results,
                  db = "KEGG_2019_Human")

plot_enrichr_dots(enrichr_results,
                  db = "WikiPathways_2019_Human",
                  width = 420,
                  height = 250)

plot_enrichr_dots(enrichr_results,
                  db = "MSigDB_Hallmark_2020",
                  width = 420)



## GSEA results ----

gsea_results %>%
  ggplot(aes(NES, -log10(padj))) +
  geom_point(alpha = .1)


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
    group_by(contrast, cluster)
  
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
        # fct_reorder(NES * is_significant, sum, na.rm = TRUE)
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
    ggplot(data_vis, aes(contrast, pathway, size = -log10(padj))) +
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
    facet_wrap(vars(cluster), nrow = 1) +
    labs(
      y = "",
      color = "normalized\nenrichment\nscore",
      size = TeX("-log_{10} (p_{adj})"),
      title = str_glue("GSEA results ({db})"),
      caption = str_glue(
        "top {top_n_positive} positively and top {top_n_negative} negatively ",
        "enriched terms per cluster; ",
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
    filename <- str_glue("dge/gsea_{db}")
  }
  ggsave_default(filename, width = 400, ...)
  p
}


plot_gsea_dots(gsea_results,
               db = "MSigDB_Hallmark_2020")

plot_gsea_dots(gsea_results,
               db = "GO_Biological_Process_2018",
               height = 500)

plot_gsea_dots(gsea_results,
               db = "KEGG_2019_Human",
               height = 300)

plot_gsea_dots(gsea_results,
               db = "WikiPathways_2019_Human",
               height = 300)

plot_gsea_dots(gsea_results,
               db = "TRRUST_Transcription_Factors_2019",
               height = 300)



## DGE in tumor cluster ----

plot_volcano(dge_results_tumor, "II_vs_IV", "dge/volcano_tumor_II_vs_IV")

plot_violin(
  dge_results_filtered_tumor,
  contrast = "II_vs_IV",
  cluster = "NB",
  direction = "up",
  top_n = Inf,
  filename = "dge/violin_tumor_II_vs_IV_up"
)
plot_violin(
  dge_results_filtered_tumor,
  contrast = "II_vs_IV",
  cluster = "NB",
  direction = "down",
  top_n = Inf,
  filename = "dge/violin_tumor_II_vs_IV_down"
)

pb_log_tumor <- aggregateData(
  nb_tumor,
  assay = "logcounts",
  fun = "mean",
  by = c("cluster_id", "sample_id")
)

plot_pbheatmap(dge_results_filtered_tumor, pb_log_tumor,
               contrast = "II_vs_IV", filename = "dge/heatmap_tumor_II_vs_IV")



## Comparson of bulk and sc data ---- 

dge_results_tumor %>% 
  filter(contrast == "II_vs_IV") %>%
  inner_join(
    read_csv("metadata/rifatbegovic2018_table_s5.csv", comment = "#")
  ) %>%
  rename(logfc_sc = logFC, q_sc = p_adj.loc, logfc_bulk = logfc, q_bulk = q) %>%
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

ggsave_default("dge/mycn_bulk_vs_sc")



## LogFC correlation ----

#' Plot a Heatmap showing correlations of log fold changes.
#'
#' @param data DGE results.
#' @param filename Name of output file.
#'
#' @return A Heatmap object.
plot_logfc_correlation_heatmap <- function(data, filename = NULL) {
  corr_mat <- 
    data %>% 
    select(gene, cluster_id, contrast, logFC) %>%
    extract(contrast, "contrast") %>% 
    unite(cluster_id, contrast, col = "cluster_group") %>% 
    pivot_wider(names_from = "cluster_group", values_from = "logFC") %>%
    select(!gene) %>%
    cor(use = "pairwise.complete.obs")
  
  p <- Heatmap(
    corr_mat,
    col = circlize::colorRamp2(
      seq(0, max(corr_mat[lower.tri(corr_mat)]), length.out = 9),
      scico(9, palette = "davos", direction = -1),
    ),
    name = "correlation"
  )
  ggsave_default(filename, plot = p, width = 130, height = 110)
  p
}

plot_logfc_correlation_heatmap(dge_results,
                               filename = "dge/logfc_correlation_heatmap")


#' Scatter plots of gene expression log-flod change in two clusters.
#'
#' @param data DGE results.
#' @param x Cluster on x-axis.
#' @param y Cluster on y-axis.
#' @param filename 
#'
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_logfc_correlation_scatter <- function(data, x, y, filename = NULL) {
  p <-
    data %>% 
    select(gene, cluster_id, contrast, logFC) %>%
    extract(contrast, "contrast") %>% 
    unite(cluster_id, contrast, col = "cluster_group") %>% 
    pivot_wider(names_from = "cluster_group", values_from = "logFC") %>% 
    ggplot(aes({{x}}, {{y}})) +
    geom_point(alpha = .1) +
    geom_smooth(method = "lm") +
    coord_fixed() +
    theme_bw()
  
  ggsave_default(filename, height = 100)
  p
}

plot_logfc_correlation_scatter(dge_results, NK_III, E_II,
                               filename = "dge/logfc_correlation_example_1")
plot_logfc_correlation_scatter(dge_results, NK_III, T_III,
                               filename = "dge/logfc_correlation_example_2")



# Publication figures -----------------------------------------------------

## Figure 3c ----

plot_violin <- function(gene,
                        cell_type,
                        groups = c("I", "II", "III", "IV"),
                        direction = c("up", "down")) {
  direction <- match.arg(direction)
  
  barcodes <- 
    nb_metadata %>% 
    filter(
      cellont_cluster %in% used_clusters,
      cellont_abbr == {{cell_type}}
    ) %>% 
    pull(cell)
  
  logcounts(nb)[gene, barcodes, drop = FALSE] %>%
    t() %>%
    as.matrix() %>%
    magrittr::set_colnames("logexp") %>%
    as_tibble(rownames = "cell") %>%
    left_join(nb_metadata, by = "cell") %>%
    filter(group %in% {{groups}}) %>%
    mutate(
      sample = rename_patients(sample),
      group = rename_groups(group)
    ) %>%
    ggplot(aes(sample, logexp)) +
    geom_violin(
      size = BASE_LINE_SIZE,
      color = "gray60",
      scale = "width",
      width = 0.8
    ) +
    geom_quasirandom(
      aes(color = group),
      width = 0.4,
      bandwidth = 1,
      alpha = 0.5,
      show.legend = FALSE,
      size = 0.2,
      shape = 16
    ) +
    stat_summary(geom = "point", fun = mean, size = .2) +
    annotate(
      "text_npc",
      label = str_glue("{cell_type}, {gene} ({direction})"),
      npcx = 0.05,
      npcy = 0.95,
      size = BASE_TEXT_SIZE_MM,
      hjust = 0
    ) +
    xlab(NULL) +
    ylab(NULL) +
    scale_color_manual(values = GROUP_COLORS) +
    theme_nb(grid = FALSE)
}

wrap_plots(
  textGrob(
    "log-normalized expression",
    rot = 90,
    gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
  ),
  plot_violin("IRF9", "B", c("I", "II"), "up"),
  plot_violin("WDR74", "B", c("I", "III"), "up"),
  plot_violin("IFI44L", "SC", c("I", "IV"), "up"),
  plot_violin("SAP30", "SC", c("I", "II"), "down"),
  plot_violin("IRF9", "B", c("I", "III"), "up"),
  plot_violin("HIST1H1E", "SC", c("I", "IV"), "down"),
  design = "ABCD\nAEFG",  
  widths = c(.1, 9, 7, 10)
)
ggsave_publication("3c_exp_violin", width = 10, height = 6, type = "png")



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
    group_by(contrast, cluster)
  
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
      contrast != "tif",
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
      contrast = factor(contrast) %>% rename_contrast()
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
    ggplot(data_vis, aes(contrast, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    circle_significant +
    xlab("vs C (contrast)") +
    ylab(NULL) +
    scale_color_gsea(
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
    facet_wrap(vars(cluster), nrow = 1) +
    theme_nb(grid = FALSE) +
    theme(
      legend.box.just = "bottom",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = "top",
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(0, 1, -3, 1, "mm"),
      panel.spacing = unit(-.5, "pt"),
    )
  
  p
}

plot_gsea(gsea_results, "MSigDB_Hallmark_2020")
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
    dge_results %>% 
    filter(contrast != "tif") %>% 
    select(gene, cluster_id, contrast, logFC) %>%
    mutate(contrast = rename_contrast(contrast)) %>% 
    unite(cluster_id, contrast, col = "cluster_group", sep = "_") %>%
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
    "NPW",     0, -4,
    "NKAIN4",  0, -3,
    
    "CXorf57", 0,  3,
    "IL7",     0,  1, 
    "MUC15",   0, -3,
    "KRT19",   1, -1
  )
  
  data <- 
    dge_results_tumor %>% 
    filter(contrast == "II_vs_IV") %>%
    inner_join(
      read_csv("metadata/rifatbegovic2018_table_s5.csv", comment = "#")
    ) %>%
    rename(
      logfc_sc = logFC,
      q_sc = p_adj.loc,
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

dge_results_filtered_gsea %>% 
  arrange(contrast, cluster_id, desc(logFC)) %>% 
  mutate(
    contrast = rename_contrast_long(contrast),
    cluster_id = CELL_TYPE_ABBREVIATIONS[cluster_id]
  ) %>% 
  select(
    "Contrast" = contrast,
    "Cell type" = cluster_id,
    "Gene" = gene,
    "Log fold change" = logFC,
    "Adjusted p-value" = p_adj.loc
  ) %>%
  save_table("S3_dge", "Microenvironment")


## Table S4 ----

gsea_results %>% 
  arrange(db, contrast, cluster, desc(NES)) %>% 
  mutate(
    contrast = rename_contrast_long(contrast),
    cluster = CELL_TYPE_ABBREVIATIONS[cluster],
    leadingEdge = map_chr(leadingEdge, str_c, collapse = ", ")
  ) %>%
  select(
    "Database" = db,
    "Contrast" = contrast,
    "Cell type" = cluster,
    "Pathway" = pathway,
    "Normalized Enrichment Score" = NES,
    "Adjusted p-value" = padj,
    "Leading edge genes" = leadingEdge
  ) %>%
  save_table("S4_gsea", "Microenvironment")
