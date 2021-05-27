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
library(openxlsx)
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
#' @param id_col colData column that contains all IDs.
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
                                   id_col,
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
    mutate({{group}} := NA_real_) %>%
    column_to_rownames("gene") %>%
    as.matrix()
 
  break_max <- quantile(dge_matrix, 0.95, na.rm = TRUE)

  heatmap <- Heatmap(
    dge_matrix,
    col = circlize::colorRamp2(
      seq(-break_max, break_max, length.out = 7),
      rev(scico(n = 7, palette = "roma"))
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    name = "log FC"
  )

  p <-
    wrap_plots(
      grid.grabExpr(draw(heatmap, merge_legend = TRUE)),
      plot_dots(
        logcounts(data),
        rownames(dge_matrix) %>% rev(),
        colData(data) %>% as_tibble() %>% pull({{id_col}})
      ) +
        scale_x_discrete(labels = function(x) str_replace(x, " ", "\n"))
    ) +
    plot_annotation(
      str_glue("Cluster {group} in {rlang::as_string(ensym(id_col))}")
    )
  ggsave_default(filename)
  p
}



# Analysis ----------------------------------------------------------------

# conmarkers_clusters <- assemble_markers(
#   nb,
#   colData(nb)$cellont_cluster,
#   pval_types = "any"
# )

# test plotting function
plot_conserved_markers(nb, conmarkers_clusters,
                       group = "NB (8)", id_col = cellont_cluster)

# make all plots
walk(
  levels(colData(nb)$cellont_cluster),
  function(cluster) {
    message("Plotting cluster ", cluster)
    plot_conserved_markers(
      nb,
      conmarkers_clusters,
      group = cluster,
      id_col = cellont_cluster,
      filename = str_glue("markers/conserved/cluster_{cluster}")
    )
  }
)


# create tables
header_style <- createStyle(textDecoration = "bold")
wb <- createWorkbook()

walk(
  levels(colData(nb)$cellont_cluster),
  function(cluster) {
    table_data <-
      conmarkers_clusters$any[[cluster]] %>% 
      as_tibble(rownames = "gene") %>% 
      filter(p.value < 1e-6) %>%
      select(gene, top = Top, p_value = p.value, FDR, logFC = t.summary.logFC)
    addWorksheet(wb, cluster)
    writeData(wb, cluster, table_data, headerStyle = header_style)
    freezePane(wb, cluster, firstRow = TRUE)
  }
)

saveWorkbook(wb, "plots/markers/conserved/genes.xlsx", overwrite = TRUE)
