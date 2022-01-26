# Plot canonical cell type and neuroblastoma markers.
# Exports several plots to plots/markers.
#
# @DEPI metadata.rds
# @DEPI rna_decontaminated.rds

library(monocle3)
library(scuttle)
library(tidyverse)
library(scico)
library(ggnewscale)
library(patchwork)
library(ComplexHeatmap)

source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")

nb_metadata <- readRDS("data_generated/metadata.rds")

colData(nb) <-
  nb_metadata %>%
  mutate(Size_Factor = colData(nb)$Size_Factor) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

markers <- read_csv("metadata/cell_markers.csv", comment = "#")



# Canonical cell type markers ---------------------------------------------

#' Plot manually selected canonical cell type markers for each cluster.
#' Oder cluster by manually assigned cell types.
#'
#' @param cluster_col Column of the cell metadata that contains cluster ids.
#' @param filename Name of the output file.
#'
#' @return A ggplot2 object
plot_canonical_markers <- function(cluster_col, filename = NULL) {
  y_annotation_data <-
    markers %>%
    arrange(desc(row_number())) %>%
    mutate(r = row_number()) %>%
    group_by(label = cell_type) %>%
    summarise(
      yintercept = first(r) - 0.5,
      label_y = mean(r)
    )

  x_annotation_data <-
    tibble(level = levels(colData(nb)$cellont_cluster)) %>% 
    extract(level, into = "label", regex = "(\\w+)") %>%
    mutate(label = as_factor(label), r = row_number()) %>%
    group_by(label) %>%
    summarize(
      xmin = first(r) - 0.5,
      xmax = last(r) + 0.5,
      label_x = mean(r)
    ) %>%
    mutate(
      fill = case_when(
        row_number() %% 2 == 0 ~ "white",
        TRUE                   ~ "gray90"
      )
    )
  
  p <-
    plot_dots(
      logcounts(nb),
      rev(markers$gene),
      colData(nb) %>% as_tibble() %>% pull({{cluster_col}}),
      panel_annotation = x_annotation_data
    ) +
    scale_x_discrete(labels = function(x) str_replace(x, " ", "\n")) +
    geom_hline(
      yintercept = y_annotation_data$yintercept,
      linetype = "dashed",
      size = 0.25
    ) +
    geom_text(
      data = y_annotation_data,
      aes(
        x = nlevels(colData(nb)$cellont_cluster) + 2,
        y = label_y,
        label = label
      ),
      size = 4,
      hjust = 0
    ) +
    theme(legend.position = "left") +
    NULL

  ggsave_default(filename, height = 297, width = 350)
  p
}

plot_canonical_markers(cellont_cluster, filename = "markers/overview")



# NB markers --------------------------------------------------------------

#' Normalize expression values.
#' 
#' Values are first limited to the given lower and upper quantile (which are
#' calculated after exclusing zero values) and then normalized to the range [0,
#' 1].
#'
#' @param x Numeric vector.
#' @param lower_quantile Quantile for calculating the lower limit.
#' @param upper_quantile Quantile for calculating the upper limit.
#'
#' @return The quantile-normalized numeric vector.
normalize_expression <- function(x, lower_quantile, upper_quantile) {
  q <- quantile(x[x > 0], c(lower_quantile, upper_quantile))
  x_limited <- pmax(pmin(x, q[2]), q[1])
  (x_limited - min(x_limited)) / (max(x_limited) - min(x_limited))
}


#' Color single cells on a dimensional reduction plot by normalized expression.
#'
#' @param counts Count matrix whose rownames correspond to features.
#' @param x Vector of x-coordinates for each cell.
#' @param y Vector of y-coordinates for each cell.
#' @param features Genes to plot.
#' @param lower_quantile Cut low expression values at this quantile.
#' @param upper_quantile Cut high expression values at this quantile.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_features <- function(counts, x, y, features,
                          lower_quantile = 0.05, upper_quantile = 0.95,
                          filename = NULL) {
  expr_data <-
    counts[features, , drop = FALSE] %>%
    Matrix::t() %>%
    apply(
      2,
      normalize_expression,
      lower_quantile = lower_quantile,
      upper_quantile = upper_quantile
    ) %>%
    as_tibble()
  
  p <-
    tibble(x = {{x}}, y = {{y}}) %>% 
    bind_cols(expr_data) %>% 
    pivot_longer(!1:2, names_to = "feature", values_to = "expression") %>%
    arrange(expression) %>%
    ggplot(aes(x, y)) +
    geom_point(aes(color = expression), size = 0.1) +
    scale_color_scico("relative\nnormalized\nexpression", palette = "bamako") +
    coord_fixed() +
    facet_wrap(vars(feature)) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  ggsave_default(filename, width = 420, height = 297)
  p
}

nb_markers <-
  markers %>%
  filter(cell_type == "NB") %>%
  pull(gene)

plot_features(
  logcounts(nb),
  x = colData(nb)$umap_1_monocle,
  y = colData(nb)$umap_2_monocle,
  nb_markers,
  filename = "markers/neuroblastoma"
)



#' Plot markers for NB cells and cells with highest scores in a selected gene
#' signature. Group cells by cluster.
#'
#' @param signature_col Column with gene signature scores.
#' @param top_prop Fraction of cells with highest scores to select.
#' @param filename Name of the output file.
#'
#' @return A ggplot2 object
plot_nb_markers_of_high_signature_cells <- function(signature_col,
                                                    top_prop = 0.05,
                                                    filename = NULL) {
  selected_cells <- union(
    colData(nb) %>% 
      as_tibble(rownames = "cell") %>% 
      slice_max(prop = top_prop, order_by = {{signature_col}}) %>% 
      pull(cell),
    colData(nb) %>% 
      as_tibble(rownames = "cell") %>% 
      filter(cellont_abbr == "NB") %>% 
      pull(cell)
  )
  
  nb_subset <- nb[, selected_cells]
  
  nb_markers <-
    markers %>% 
    filter(cell_type == "NB") %>% 
    pull(gene)
  
  p <- plot_dots(
    logcounts(nb_subset),
    nb_markers,
    colData(nb_subset)$cluster_50
  )
  
  ggsave_default(filename, height = 100)
  p
}

plot_nb_markers_of_high_signature_cells(signature_mesenchymal,
                                        filename = "markers/nb_mesenchymal")
plot_nb_markers_of_high_signature_cells(signature_ncc_like,
                                        filename = "markers/nb_ncc_like")



# PanglaoDB ---------------------------------------------------------------

panglaodb <- 
  read_tsv("metadata/PanglaoDB_markers_27_Mar_2020.tsv.gz") %>% 
  set_names(colnames(.) %>% str_replace_all(" ", "_")) %>%
  filter(
    str_detect(species, "Hs"),
    canonical_marker == 1
  ) %>% 
  mutate(cell_type = str_replace(cell_type, "/", " or "))

panglaodb_cell_types <- 
  panglaodb %>%
  filter(
    organ %in% c("Immune system", "Blood", "Brain") |
    cell_type == "Hematopoietic stem cells"
  ) %>% 
  count(organ, cell_type)


#' Plot canonical markers as defined by the PanglaoDB database in a dotplot
#' and highlight manual cluster classifications.
#'
#' @param cell_type Cell type available in the database.
#' @param cluster_col Column of the cell metadata that contains cluster ids.
#' @param save_plot If `TRUE`, save the plot in an automatically derived folder.
#'
#' @return A ggplot object
plot_panglaodb_markers <- function(cell_type, cluster_col,
                                   folder = "panglaodb", save_plot = TRUE) {
  info("Plotting markers for {cell_type}")
  
  counts <- logcounts(nb)
  
  markers <-
    panglaodb %>%
    filter(cell_type == {{cell_type}}) %>%
    pull(official_gene_symbol) %>%
    unique() %>% 
    str_sort(decreasing = TRUE)
  
  organ <- 
    panglaodb %>%
    filter(cell_type == {{cell_type}}) %>%
    slice(1) %>% 
    pull(organ)
  
  known_markers <- intersect(markers, rownames(counts))
  missing_markers <- setdiff(markers, rownames(counts))
  if (length(missing_markers) > 0)
    missing_markers <- str_glue(
      "Not observed: {str_c(missing_markers, collapse = ', ')}"
    )
  else
    missing_markers <- waiver()

  figure_height <- max(6 * length(known_markers), 70)
  
  x_annotation_data <-
    tibble(level = levels(colData(nb)$cellont_cluster)) %>% 
    extract(level, into = "label", regex = "(\\w+)") %>%
    mutate(label = as_factor(label), r = row_number()) %>%
    group_by(label) %>%
    summarize(
      xmin = first(r) - 0.5,
      xmax = last(r) + 0.5,
      label_x = mean(r)
    ) %>%
    mutate(
      fill = case_when(
        row_number() %% 2 == 0 ~ "white",
        TRUE                   ~ "gray90"
      )
    )
  
  p <-
    plot_dots(
      counts,
      known_markers,
      colData(nb) %>% as_tibble() %>% pull({{cluster_col}}),
      panel_annotation = x_annotation_data
    ) +
    scale_x_discrete(labels = function(x) str_replace(x, " ", "\n")) +
    labs(
      title = cell_type,
      subtitle = "PanglaoDB canonical markers",
      caption = missing_markers
    )

  if (save_plot)
    ggsave_default(str_glue("markers/{folder}/{organ}/{cell_type}"),
                   height = figure_height, width = 250)
  p
}


plot_panglaodb_markers("B cells", cellont_cluster, folder = "panglaodb")

walk(
  panglaodb_cell_types$cell_type,
  plot_panglaodb_markers,
  cluster_50,
  folder = "panglaodb"
)



# Gene expression heatmaps ------------------------------------------------

plot_celltype_genes <- function(cell_type,
                                cells_per_patient = 100,
                                n_genes = 500) {
  metadata <-
    nb_metadata %>% 
    filter(cellont_abbr == {{cell_type}})
  
  nb_sub <- nb[, metadata$cell]
  nb_sub <- nb_sub[rowSums(counts(nb_sub) > 0) > 0, ]
  
  mat_metadata <- 
    metadata %>%
    group_by(group, sample) %>% 
    slice_sample(n = cells_per_patient) %>% 
    arrange(sample) %>% 
    mutate(group = rename_groups(group), sample = rename_patients(sample))
  
  mat <-
    nb_sub[
      nb_sub %>%
        logcounts() %>%
        rowSums() %>%
        sort(TRUE) %>%
        names() %>%
        magrittr::extract(1:n_genes),
      mat_metadata$cell
    ] %>% 
    logcounts() %>% 
    as.matrix()
  
  
  Heatmap(
    mat,
    col = circlize::colorRamp2(
      seq(0, quantile(mat, 0.95), length.out = 9),
      scico(9, palette = "oslo", direction = -1),
    ),
    name = "log expression",
    
    show_row_names = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    
    show_column_names = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    column_split = mat_metadata$sample,
    show_column_dend = FALSE,
    show_parent_dend_line = FALSE,
    
    top_annotation = HeatmapAnnotation(
      group = mat_metadata$group,
      col = list(
        group = GROUP_COLORS
      )
    )
  )
}

(p <- plot_celltype_genes("NK"))
ggsave_default("markers/genes_NK", plot = p)

(p <- plot_celltype_genes("M"))
ggsave_default("markers/genes_M", plot = p)



# NB markers in patients --------------------------------------------------

plot_dots(
  logcounts(nb),
  markers %>%
    filter(cell_type == "NB") %>%
    pull(gene),
  if_else(
    colData(nb)$cellont_abbr == "NB",
    as.character(colData(nb)$sample),
    as.character(colData(nb)$cellont_abbr)
  ) %>% 
    rename_patients() %>% 
    as_factor() %>% 
    fct_relevel(PATIENT_ORDER, names(CELL_TYPE_ABBREVIATIONS))
)
ggsave_default("markers/nb_markers_patients", height = 80)

