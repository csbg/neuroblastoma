# Find conserved markers.
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
  
  clusters <- tribble(
    ~id, ~type,
    "3",  "T",
    "5",  "T",
    "18", "T",
    "4",  "NK",
    "6",  "NK",
    "2",  "B",
    "9",  "B",
    "12", "B",
    "17", "B",
    "19", "B",
    "21", "B",
    "1",  "M",
    "15", "M",
    "16", "M",
    "22", "M",
    "14", "D",
    "13", "E",
    "7",  "other",
    "10", "other",
    "11", "other",
    "20", "other",
    "8",  "NB"
  )
  
  colData(data)[[dotplot_groups]] <- 
    colData(data)[[dotplot_groups]] %>% 
    fct_relevel(clusters$id)

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
    magrittr::extract(, clusters$id)

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
        clusters %>%
        filter(id %in% colnames(dge_matrix)) %>%
        pull(type),
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
    plot_genes_by_group(
      data,
      rownames(dge_matrix) %>% rev(),
      group_cells_by = dotplot_groups,
      ordering_type = "none"
    ) +
      scale_color_scico(
        name = "log(mean + 0.1)",
        palette = "acton",
        direction = -1
      ) +
      scale_radius(name = "percentage", range = c(0, 6)) +
      xlab(NULL) +
      ylab(NULL) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  ) +
    plot_layout(widths = c(1.5, 1)) +
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

plot_conserved_markers(nb, conmarkers_clusters, group = "8",
                       filename = "markers/conserved_markers_c8")

walk(
  levels(colData(nb)$cluster_50),
  function(cluster) {
    message("Plotting cluster ", cluster)
    plot_conserved_markers(
      nb,
      conmarkers_clusters,
      group = cluster,
      dotplot_groups = "cluster_50",
      filename = str_glue("markers/conserved_markers_c{cluster}")
    ) 
  }
)
