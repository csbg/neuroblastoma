# Find conserved markers.
#
# @DEPI rna_integrated_monocle.rds
# @DEPI metadata.rds


library(monocle3)
library(scran)
library(scuttle)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
source("common_functions.R")




# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_integrated_monocle.rds") %>% 
  logNormCounts()

nb@colData <-
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
                                   dotplot_groups = "cluster",
                                   n_genes = 30,
                                   pval_max = 1e-6,
                                   pval_type = c("any", "some", "all"),
                                   filename = NULL) {
  pval_type <- match.arg(pval_type)
  sort_col <- "t.summary.logFC"
  prefix <- "t.logFC."
  
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
    column_to_rownames("gene") %>%
    as.matrix()
  
  break_max <- quantile(dge_matrix, 0.95)

  p <- wrap_plots(
    pheatmap(
      dge_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      breaks = seq(-break_max, break_max, length.out = 101),
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      silent = TRUE
    )$gtable,
    plot_genes_by_group(
      data,
      rownames(dge_matrix) %>% rev(),
      group_cells_by = dotplot_groups,
      ordering_type = "none"
    )
  ) +
    plot_annotation(str_glue("Group {group} in {dotplot_groups}"))
  ggsave_default(filename)
  p
}



# Analysis ----------------------------------------------------------------

## Clusters ----

conmarkers_clusters <-
  assemble_markers(nb, colData(nb)$cluster_50, pval_types = "any")

plot_conserved_markers(nb, conmarkers_clusters,
                       group = "8", dotplot_groups = "cluster_50",
                       filename = "monocle/markers/conserved_markers_c_8")

levels(colData(nb)$cluster_50) %>% 
  walk(
    function(cluster) {
      message("Plotting cluster ", cluster)
      plot_conserved_markers(
        nb,
        conmarkers_clusters,
        group = cluster,
        dotplot_groups = "cluster_50",
        filename = str_glue("monocle/markers/conserved_markers/c_{cluster}")
      ) 
    }
  )


## Partitions ----

conmarkers_partitions <-
  assemble_markers(nb, colData(nb)$partition_50, pval_types = "any")

plot_conserved_markers(nb, conmarkers_partitions,
                       group = "1", dotplot_groups = "partition_50",
                       filename = "monocle/markers/conserved_markers_p_1")

levels(colData(nb)$partition_50) %>% 
  walk(
    function(partition) {
      message("Plotting partition ", partition)
      plot_conserved_markers(
        nb,
        conmarkers_partitions,
        group = partition,
        dotplot_groups = "partition_50",
        filename = str_glue("monocle/markers/conserved_markers/p_{partition}")
      ) 
    }
  )
