# Find cluster markers, generate plots and export top markers to a CSV file.
#
# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds


library(monocle3)
library(scran)
library(scuttle)
library(ComplexHeatmap)
library(tidyverse)
library(patchwork)
library(scico)
source("common_functions.R")




# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts()

colData(nb) <-
  readRDS("data_generated/metadata.rds") %>%
  mutate(Size_Factor = colData(nb)$Size_Factor) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)



# Functions ---------------------------------------------------------------

#' Find conserved markers and combine results from several statistical tests.
#'
#' @param data A monocle dataset.
#' @param cluster_ids A vector of length `ncol(data)` containing cluster IDs,
#'   typically a metadata column.
#' @param pval_types Value for `pval.types` of `scran::findMarkers`.
#'
#' @return A named list with names equal to the values in `pval_types`. Each
#'   element is a data frame as returend by `scran::multiMarkerStats()`.
assemble_markers <- function(data,
                             cluster_ids,
                             pval_types = c("any", "some", "all")) {
  test_types <- c("t", "wilcox", "binom")
  
  message("Finding markers ...")
  map(
    pval_types,
    function(pval_type) {
      message("  pval.type=", pval_type)
      test_results <- map(
        test_types,
        function(test_type) {
          message("    test.type=", test_type)
          findMarkers(
            data,
            groups = cluster_ids,
            direction = "up",
            pval.type = pval_type,
            test.type = test_type
          )
        }
      ) %>% 
        set_names(test_types)
      multiMarkerStats(
        t = test_results$t,
        wilcox = test_results$wilcox,
        binom = test_results$binom
      )
    }
  ) %>% 
    set_names(pval_types)
}


#' Plot a DGE heatmap and dot plot for a selected ID of a cluster/partition.
#'
#' @param data A monocle dataset.
#' @param markers Results from `assemble_markers()`.
#' @param group Cluster/partition that should be plotted.
#' @param dotplot_groups colData column that contains all IDs.
#' @param n_genes NUmber of genes to plot.
#' @param pval_max Maximum allowed p value.
#' @param pval_type Selects value for `pval.types` of `scran::findMarkers` for
#'   which results should be plotted.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_conserved_markers <- function(data,
                                   markers,
                                   group,
                                   dotplot_groups = "cluster_50",
                                   n_genes = 30,
                                   pval_max = 1e-6,
                                   pval_type = c("any", "some", "all"),
                                   filename = NULL) {
  pval_type <- match.arg(pval_type)
  sort_col <- "t.summary.logFC"
  prefix <- "t.logFC."

  cluster_info <-
    read_csv("metadata/clusters.csv") %>%
    mutate(
      cell_type =
        as_factor(cell_type) %>%
        fct_relevel("T", "NK", "B", "M", "D", "E", "NB", "other")
    ) %>%
    arrange(cell_type, cluster) %>%
    mutate(cluster = as_factor(cluster) %>% fct_inorder())

  dge_df <-
    markers[[pval_type]][[group]] %>%
    as_tibble(rownames = "gene") %>%
    filter(p.value < pval_max) %>%
    slice_head(n = n_genes) %>%
    arrange(desc(!!sym(sort_col)))

  dge_matrix <-
    dge_df %>%
    select(gene, starts_with(prefix)) %>%
    rename_with(~str_replace(.x, prefix, "")) %>%
    mutate({{group}} := NA_real_) %>%
    column_to_rownames("gene") %>%
    as.matrix() %>%
    magrittr::extract(, levels(cluster_info$cluster))

  break_max <- quantile(dge_matrix, 0.95, na.rm = TRUE)

  heatmap <- Heatmap(
    dge_matrix,
    col = circlize::colorRamp2(
      seq(-break_max, break_max, length.out = 7),
      rev(scico(n = 7, palette = "roma"))
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "expression",
    top_annotation = HeatmapAnnotation(
      "cell type" =
        cluster_info %>%
        filter(cluster %in% colnames(dge_matrix)) %>%
        pull(cell_type),
      col = list("cell type" = c(
        "T" = "#1f78b4",
        NK = "#a6cee3",
        B = "#33a02c",
        M = "#ff7f00",
        NB = "#e31a1c",
        E = "#b15928",
        other = "#6a3d9a",
        D = "#fdbf6f"
      ))
    )
  )

  p <- wrap_plots(
    grid.grabExpr(draw(heatmap, merge_legend = TRUE)),
    plot_dots(
      logcounts(data),
      rownames(dge_matrix) %>% rev(),
      colData(data)[[dotplot_groups]] %>%
        fct_relevel(levels(cluster_info$cluster))
    )
  ) +
    plot_annotation(str_glue("Group {group} in {dotplot_groups}"))
  ggsave_default(filename)
  p
}



# Analysis ----------------------------------------------------------------

conmarkers_clusters <- assemble_markers(
  nb,
  colData(nb)$cluster_50,
  pval_types = "any"
)

plot_conserved_markers(nb, conmarkers_clusters, group = "9",
                       filename = "markers/conserved_markers_c9")

walk(
  levels(colData(nb)$cluster_50),
  function(cluster) {
    message("Plotting cluster ", cluster)
    plot_conserved_markers(
      nb,
      conmarkers_clusters,
      group = cluster,
      filename = str_glue("markers/cluster_{cluster}")
    )
    
    fs::dir_create("plots/markers/tables")
    conmarkers_clusters$any[[cluster]] %>% 
      as_tibble(rownames = "gene") %>% 
      filter(p.value < 1e-6) %>%
      select(gene, Top, p.value, FDR, t.summary.logFC) %>% 
      write_csv(str_glue("plots/markers/tables/cluster_{cluster}.csv"))
  }
)



