# @DEPI dge_pb_results.rds

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

dge <- readRDS("data_generated/dge_pb_results.rds")



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

plot_volcano(dge$results, "II_vs_I", filename = "dge_pb/volcano_II")
plot_volcano(dge$results, "III_vs_I", filename = "dge_pb/volcano_III")
plot_volcano(dge$results, "IV_vs_I", filename = "dge_pb/volcano_IV")



## Violin plot ----

#' Title
#'
#' @param data A DGE dataset.
#' @param contrast A contrast and ...
#' @param cluster ... cluster for which all samples should be plotted.
#' @param cds Cell dataset.
#' @param top_n Plot the top differentially expressed genes.
#' @param direction Plot up- or downregulated genes.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_violin <- function(data,
                        contrast,
                        cluster,
                        cds = dge$cds,
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
      cds[, cds$cluster_id == cluster &
            cds$group_id %in% c(group_left, group_right)],
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
#             filename = "dge_pb/violin_II_vs_I_c1_up")

dge$results_filtered %>% 
  distinct(contrast, cluster = as.character(cluster_id), direction) %>% 
  pwalk(
    function(contrast, cluster, direction) {
      info("Plotting {direction}regulated genes ",
           "in cluster {cluster} in contrast {contrast}")
      plot_violin(
        dge$results_filtered,
        top_n = Inf,
        contrast = contrast,
        cluster = cluster,
        direction = direction,
        filename = str_glue("dge_pb/violin_{contrast}_c{cluster}_{direction}")
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
    group_map(
      ~assays(pb_data)[[.y$cluster_id]][.x$gene, sample_data$sample]
    ) %>%
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
  dge$cds,
  assay = "logcounts",
  fun = "mean",
  by = c("cluster_id", "sample_id")
)

plot_pbheatmap(dge$results_filtered, pb_log,
               contrast = "II_vs_I", filename = "dge_pb/heatmap_II")
plot_pbheatmap(dge$results_filtered, pb_log,
               contrast = "III_vs_I", filename = "dge_pb/heatmap_III")
plot_pbheatmap(dge$results_filtered, pb_log,
               contrast = "IV_vs_I", filename = "dge_pb/heatmap_IV")



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
    filename <- str_glue("dge_pb/enrichr_{db}")
  }
  ggsave_default(filename, ...)
  p
}


plot_enrichr_dots(dge$enrichr,
                  db = "GO_Biological_Process_2018",
                  height = 400,
                  width = 400)

plot_enrichr_dots(dge$enrichr,
                  db = "GO_Cellular_Component_2018")

plot_enrichr_dots(dge$enrichr,
                  db = "GO_Molecular_Function_2018",
                  width = 420)

plot_enrichr_dots(dge$enrichr,
                  db = "KEGG_2019_Human")

plot_enrichr_dots(dge$enrichr,
                  db = "WikiPathways_2019_Human",
                  width = 420,
                  height = 250)

plot_enrichr_dots(dge$enrichr,
                  db = "MSigDB_Hallmark_2020",
                  width = 420)



## GSEA results ----

dge$gsea %>%
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
    filename <- str_glue("dge_pb/gsea_{db}")
  }
  ggsave_default(filename, width = 400, ...)
  p
}


plot_gsea_dots(dge$gsea,
               db = "MSigDB_Hallmark_2020")

plot_gsea_dots(dge$gsea,
               db = "GO_Biological_Process_2018",
               height = 500)

plot_gsea_dots(dge$gsea,
               db = "KEGG_2019_Human",
               height = 300)

plot_gsea_dots(dge$gsea,
               db = "WikiPathways_2019_Human",
               height = 300)

plot_gsea_dots(dge$gsea,
               db = "TRRUST_Transcription_Factors_2019",
               height = 300)



## DGE in tumor cluster ----

plot_volcano(
  dge$results_tumor,
  "II_vs_IV",
  filename = "dge_pb/volcano_tumor_II_vs_IV"
)

plot_violin(
  dge$results_tumor_filtered,
  contrast = "II_vs_IV",
  cluster = "NB",
  direction = "up",
  top_n = Inf,
  filename = "dge_pb/violin_tumor_II_vs_IV_up"
)
plot_violin(
  dge$results_tumor_filtered,
  contrast = "II_vs_IV",
  cluster = "NB",
  direction = "down",
  top_n = Inf,
  filename = "dge_pb/violin_tumor_II_vs_IV_down"
)

pb_log_tumor <- aggregateData(
  dge$cds_tumor,
  assay = "logcounts",
  fun = "mean",
  by = c("cluster_id", "sample_id")
)

plot_pbheatmap(
  dge$results_tumor_filtered,
  pb_log_tumor,
  contrast = "II_vs_IV",
  filename = "dge_pb/heatmap_tumor_II_vs_IV"
)



## Comparson of bulk and sc data ---- 

dge$results_tumor %>% 
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

ggsave_default("dge_pb/mycn_bulk_vs_sc")



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

plot_logfc_correlation_heatmap(dge$results,
                               filename = "dge_pb/logfc_correlation_heatmap")


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

plot_logfc_correlation_scatter(dge$results, NK_III, E_II,
                               filename = "dge_pb/logfc_correlation_example_1")
plot_logfc_correlation_scatter(dge$results, NK_III, T_III,
                               filename = "dge_pb/logfc_correlation_example_2")