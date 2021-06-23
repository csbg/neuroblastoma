# Misc analyses and plots.
#
# @DEPI metadata.rds

library(ComplexHeatmap)
library(RColorBrewer)
library(ggpmisc)
library(patchwork)
library(tidyverse)
library(egg)
library(fs)
library(bluster)
library(ggbeeswarm)
library(latex2exp)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

nb_data <- readRDS("data_generated/metadata.rds")



# Clusters ----------------------------------------------------------------

#' Plot all clusters.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param label_direct If `TRUE`, print cluster labels at cluster mean
#' @param color_scale A ggplot2 color scale used for coloring the points.
#'   If `NULL`, use the default discrete color scale.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_clusters_all <- function(data, x, y, clusters, label_direct = TRUE,
                              color_scale = NULL, filename = NULL) {
  color_scale <-
    color_scale %||%
    scale_color_hue(guide = guide_legend(override.aes = list(size = 5)))
  
  cluster_labels <- 
    data %>% 
    group_by(label = {{clusters}}) %>% 
    summarise({{x}} := mean({{x}}), {{y}} := mean({{y}}))
  
  p <- 
    data %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_point(
      aes(color = {{clusters}}),
      size = .1,
      show.legend = !label_direct
    ) +
    {
      if (label_direct)
        geom_text(data = cluster_labels, aes(label = label))
    } +
    color_scale +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  ggsave_default(filename)
  p
}

plot_clusters_all(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                  filename = "seurat/clusters_all_umap_0.5")
plot_clusters_all(nb_data, tsne_1_seurat, tsne_2_seurat, cluster_0.5,
                  filename = "seurat/clusters_all_tsne_0.5")

plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, cluster_20,
                  filename = "monocle/clusters_all_umap_20")
plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
                  filename = "monocle/clusters_all_umap_50")
plot_clusters_all(nb_data, tsne_1_monocle, tsne_2_monocle, cluster_50,
                  filename = "monocle/clusters_all_tsne_50")
plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, partition_20,
                  filename = "monocle/partitions_all_umap_20")
plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, partition_50,
                  filename = "monocle/partitions_all_umap_50")
plot_clusters_all(nb_data, tsne_1_monocle, tsne_2_monocle, partition_50,
                  filename = "monocle/partitions_all_tsne_50")



#' Plot all clusters, facet by sample.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param sample Column with sample IDs.
#' @param nrow Number of rows.
#' @param show_legend If `TRUE`, include color legend.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_clusters_per_sample <- function(data, x, y, clusters, sample,
                                     nrow = NULL, show_legend = FALSE,
                                     filename = NULL) {
  p <- 
    data %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_point(
      aes(color = {{clusters}}),
      size = .01,
      shape = 20,
      show.legend = show_legend
    ) +
    facet_wrap(
      vars(group, {{sample}}),
      nrow = nrow,
      labeller = function(labels) label_value(labels, multi_line = FALSE)
    ) +
    scale_color_hue(guide = guide_legend(override.aes = list(size = 5))) +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    NULL
  
  ggsave_default(filename, width = 420, height = 297)
  p
}

plot_clusters_per_sample(nb_data, umap_1_monocle, umap_2_monocle,
                         cluster_50, sample,
                         nrow = 3, filename = "monocle/clusters_sample_umap_50")

plot_clusters_per_sample(nb_data %>% mutate(group2 = group),
                         umap_1_monocle, umap_2_monocle,
                         cluster_50, group2,
                         nrow = 2, filename = "monocle/clusters_groups_umap_50")


#' Highlight an individual cluster over all cells on the background (light
#' gray), facet by cluster.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param nrow Number of rows.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_clusters_highlight <- function(data, x, y, clusters,
                                    nrow = NULL, filename = NULL) {
  p <- 
    data %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_point(
      color = "gray90",
      data = data %>% select(!{{clusters}}),
      size = .01,
      shape = 20
    ) +
    geom_point(
      aes(color = {{clusters}}),
      size = .01,
      shape = 20,
      show.legend = FALSE
    ) +
    facet_wrap(vars({{clusters}}), nrow = nrow) +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  ggsave_default(filename, width = 420, height = 297)
  p
}

plot_clusters_highlight(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
                        nrow = 5,
                        filename = "monocle/clusters_highlight_umap_50")



#' Highlight an individual cluster, generate one plot per cluster.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param folder Folder where individual plots will be saved (filenames
#'   correspond to cluster IDs).
#'
#' @return `NULL`
plot_clusters_selected <- function(data, x, y, clusters, folder = NULL) {
  walk(
    data %>% pull({{clusters}}) %>% levels(),
    function(cluster) {
      message("Plotting cluster ", cluster)
      data_hl <-
        data %>%
        filter({{clusters}} == {{cluster}})
      
      n_cells <- nrow(data_hl)
      
      p <- 
        data %>%
        ggplot(aes({{x}}, {{y}})) +
        geom_point(
          color = "gray90",
          size = .01,
          shape = 20
        ) +
        geom_point(
          data = data_hl,
          size = .01,
          shape = 20
        ) +
        annotate(
          "text_npc", npcx = 0.05, npcy = 0.95, size = 6,
          label = str_glue("cluster {cluster}, {n_cells} cells")
        ) +
        coord_fixed() +
        theme_classic() +
        NULL
      
      if (!is.null(folder))
        ggsave_default(str_glue("{folder}/{cluster}"),
                       width = 200, height = 200, crop = FALSE)
    }
  )
}

plot_clusters_selected(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
                       folder = "monocle/clusters_selected_umap_50")



# Cluster diagnostics -----------------------------------------------------

#' Plot silhouettes and neighborhood purity as calculated by `bluster`.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_diagnostics <- function(data, x, y, clusters, filename = NULL) {
  mat <- 
    nb_data %>%
    select({{x}}, {{y}}) %>%
    as.matrix()
  
  clusters <-
    data %>%
    pull({{clusters}})
  
  p1 <-
    approxSilhouette(mat, clusters) %>% 
    as_tibble(rownames = "cell") %>% 
    mutate(
      cluster = clusters,
      closest = if_else(width > 0, cluster, other)
    ) %>% 
    ggplot(aes(cluster, width, color = closest)) +
    geom_quasirandom(size = 0.1, show.legend = FALSE) +
    geom_hline(yintercept = 0) +
    ggtitle("Silhouette")
  
  p2 <-
    neighborPurity(mat, clusters) %>% 
    as_tibble(rownames = "cell") %>% 
    mutate(
      cluster = clusters,
      maximum = as_factor(maximum)
    ) %>% 
    ggplot(aes(cluster, purity, color = maximum)) +
    geom_quasirandom(alpha = .1, show.legend = FALSE) +
    ggtitle("Neighborhood purity")
  
  p <- wrap_plots(p1, p2)
  ggsave_default(filename)
  p
}

plot_cluster_diagnostics(nb_data, umap_1_monocle, umap_2_monocle, cluster_20,
                         filename = "monocle/cluster_diagnostics_20")

plot_cluster_diagnostics(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
                         filename = "monocle/cluster_diagnostics_50")

plot_cluster_diagnostics(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                         filename = "seurat/cluster_diagnostics_0.5")



# Samples -----------------------------------------------------------------

plot_clusters_all(nb_data %>% arrange(umap_1_seurat),
                  umap_1_seurat, umap_2_seurat, sample,
                  label_direct = FALSE, filename = "seurat/samples_all_umap")
plot_clusters_all(nb_data %>% arrange(umap_1_monocle),
                  umap_1_monocle, umap_2_monocle, sample,
                  label_direct = FALSE, filename = "monocle/samples_all_umap")
plot_clusters_all(nb_data,
                  umap_1_unaligned, umap_2_unaligned, sample,
                  label_direct = FALSE,
                  filename = "monocle/samples_all_umap_unaligned")



# Subclustering -----------------------------------------------------------

#' Plot all sub clusters.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param superclusters Column with supercluster IDs.
#' @param subclusters Column with subcluster IDs.
#' @param label_direct If `TRUE`, print cluster labels at cluster mean.
#' @param color_scale A ggplot2 color scale used for coloring the points.
#'   If `NULL`, use the default discrete color scale.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_subclusters <- function(data, x, y, superclusters, subclusters,
                             label_direct = TRUE, color_scale = NULL,
                             filename = NULL) {
  cluster_labels <-
    data %>%
    group_by({{superclusters}}, label = {{subclusters}}) %>%
    summarise({{x}} := mean({{x}}), {{y}} := mean({{y}}))
  
  color_scale <-
    color_scale %||%
    scale_color_hue(guide = guide_legend(override.aes = list(size = 5)))
  
  p <-
    data %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_point(
      aes(color = {{subclusters}}),
      size = .01,
      shape = 20,
      show.legend = !label_direct
    ) +
    {
      if (label_direct)
        geom_text(data = cluster_labels, aes(label = label), size = 3)
    } +
    color_scale +
    facet_wrap(vars({{superclusters}}), ncol = 6) +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  ggsave_default(filename, width = 420, height = 297)
  p
}


plot_subclusters(nb_data, umap_1_subcluster, umap_2_subcluster,
                 cluster_50, subcluster_20,
                 filename = "monocle/subclusters_50_20")

nb_data %>% 
  mutate(
    cluster_label =
      str_glue("{cluster_50} ({cellont_name})") %>%
      as_factor() %>% 
      fct_reorder(as.integer(cluster_50))
  ) %>% 
  plot_subclusters(umap_1_subcluster, umap_2_subcluster,
                   cluster_label, subcluster_20,
                   filename = "monocle/subclusters_50_20")



# Cell types --------------------------------------------------------------

#' Preprocess cell types: Convert to factor, order by frequency, possibly lump.
#'
#' @param data Metadata.
#' @param ref Reference cell type dataset (hpca, blueprint, dice, dmap, monaco).
#' @param label Label granularity (broad fine).
#' @param prop Lump cell types with a lower relative abundance.
#' @param n Lump cell types except for the most frequent ones. If `prop` and `n`
#'   are `NULL` (default value), do nothing.
#'
#' @return A data frame with a new column `cell_type`.
preprocess_celltypes <- function(data, ref, label, prop = NULL, n = NULL) {
  if (!is.null(prop)) {
    if (!is.null(n))
      stop("Values were specified for both prop and n.")
    lump <- partial(fct_lump_prop, prop = prop)
  } else if (!is.null(n)) {
    lump <- partial(fct_lump_n, n = n)
  } else {
    lump <- function(x) x
  }
  
  cell_type_column <- rlang::sym(str_glue("cell_type_{ref}_{label}"))
  
  data %>% 
    mutate(
      cell_type =
        as_factor(!!cell_type_column) %>%
        fct_infreq() %>% 
        lump() %>% 
        fct_explicit_na("Unknown")
    )
}


#' Plot a UMAP with cell types highlighted.
#'
#' @param data Metadata.
#' @param ref   |
#' @param label |
#' @param prop  ↓
#' @param n Passed to `preprocess_celltypes()`.
#' @param ... Passed to `plot_clusters_all()`.
#'
#' @return A ggplot object.
plot_celltypes_all <- function(data, ref, label,
                               prop = NULL, n = NULL, ...) {
  colors <- list(
    T4 = "#1f78b4",
    T8 = "#62a3cb",
    NK = "#a6cee3",
    B = "#33a02c",
    B_prec = "#b2df8a",
    mono = "#ff7f00",
    neuron = "#e31a1c",
    ery = "#b15928",
    hsc = "#6a3d9a",
    gran = "#cab2d6",
    den = "#fdbf6f",
    other = "black",
    na = "gray80"
  )
  
  cell_type_colors <- c(
    T_cells = colors$T4,
    "T cells" = colors$T4,
    "CD4+ T cells" = colors$T4,
    "CD4+ T-cells" = colors$T4,
    "T cells, CD4+" = colors$T4,
    "CD8+ T cells" = colors$T8,
    "CD8+ T-cells" = colors$T8,
    "T cells, CD8+" = colors$T8,
    
    NK_cell = colors$NK,
    "NK cells" = colors$NK,
    
    B_cell = colors$B,
    "B-cells" = colors$B,
    "B cells" = colors$B,
    
    "Pro-B_cell_CD34+" = colors$B_prec,
    "Pre-B_cell_CD34-" = colors$B_prec,
    
    Monocyte = colors$mono,
    Monocytes = colors$mono,
    
    Neurons = colors$neuron,
    
    Erythroblast = colors$ery,
    Erythrocytes = colors$ery,
    "Erythroid cells" = colors$ery,
    
    GMP = colors$hsc,
    CMP = colors$hsc,
    "Pro-Myelocyte" = colors$hsc,
    BM = colors$hsc,
    HSC = colors$hsc,
    HSCs = colors$hsc,
    MEPs = colors$hsc,
    Progenitors = colors$hsc,
    
    Granulocytes = colors$gran,
    Basophils = colors$gran,
    Neutrophils = colors$gran,
    
    "Dendritic cells" = colors$den,
    
    Other = colors$other,
    Unknown = colors$na
  )
  
  data %>% 
    preprocess_celltypes(ref, label, prop, n) %>% 
    plot_clusters_all(
      umap_1_monocle,
      umap_2_monocle,
      cell_type,
      label_direct = FALSE,
      color_scale = scale_color_manual(
        values = cell_type_colors,
        guide = guide_legend(override.aes = list(size = 5))
      ),
      ...
    )
}

plot_celltypes_all(nb_data, "hpca", "broad", prop = 0.01,
                   filename = "cell_types/celltype_all_hpca")
plot_celltypes_all(nb_data, "blueprint", "broad", prop = 0.01,
                   filename = "cell_types/celltype_all_blueprint")
plot_celltypes_all(nb_data, "dice", "broad",
                   filename = "cell_types/celltype_all_dice")
plot_celltypes_all(nb_data, "dmap", "broad", prop = 0.01,
                   filename = "cell_types/celltype_all_dmap")
plot_celltypes_all(nb_data, "monaco", "broad",
                   filename = "cell_types/celltype_all_monaco")



#' Make subplots, each of which highlights a particular cell type.
#'
#' @param data Metadata.
#' @param ref   |
#' @param label |
#' @param prop  ↓
#' @param n Passed to `preprocess_celltypes()`.
#' @param ... Passed to `plot_clusters_all()`.
#'
#' @return A ggplot object.
plot_celltypes_highlight <- function(data, ref, label,
                                     prop = NULL, n = NULL, ...) {
  data %>% 
    preprocess_celltypes(ref, label, prop, n) %>% 
    plot_clusters_highlight(
      umap_1_monocle,
      umap_2_monocle,
      cell_type,
      ...
    )
}

plot_celltypes_highlight(nb_data, "hpca", "broad", prop = 0.01,
                         filename = "cell_types/celltype_highlight_hpca")
plot_celltypes_highlight(nb_data, "blueprint", "broad", prop = 0.01,
                         filename = "cell_types/celltype_highlight_blueprint")
plot_celltypes_highlight(nb_data, "dice", "broad",
                         filename = "cell_types/celltype_highlight_dice")
plot_celltypes_highlight(nb_data, "dmap", "broad", prop = 0.01,
                         filename = "cell_types/celltype_highlight_dmap")
plot_celltypes_highlight(nb_data, "monaco", "broad",
                         filename = "cell_types/celltype_highlight_monaco")



# Cell types vs clusters --------------------------------------------------

#' Plot a heatmap of cell type vs cluster.
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param ref   |
#' @param label |
#' @param prop  ↓
#' @param n Passed to `preprocess_celltypes()`.
#' @param filename Name of output file.
#' @param ... Passed to `Heatmap()`.
#'
#' @return A ggplot object.
plot_cvt_heatmap <- function(data, clusters, ref, label,
                             prop = NULL, n = NULL, filename = NULL, ...) {
  
  mat <-
    data %>% 
    preprocess_celltypes(ref, label, prop, n) %>% 
    mutate(cluster = {{clusters}}) %>%
    count(cluster, cell_type, .drop = FALSE) %>%
    group_by(cluster) %>%
    mutate(n_rel = n / sum(n)) %>%
    select(!n) %>%
    ungroup() %>%
    pivot_wider(names_from = "cell_type", values_from = "n_rel") %>%
    column_to_rownames("cluster") %>%
    as.matrix() %>%
    t() %>%
    replace_na(0)

  p <- Heatmap(
    mat,
    col = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
    show_heatmap_legend = FALSE,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_rows = FALSE,
    column_names_rot = 0
  )

  ggsave_default(filename, plot = p)
  p
}

plot_cvt_heatmap(nb_data, cluster_50, "hpca", "broad",
                 filename = "cell_types/cvt_heatmap_hpca_broad")
plot_cvt_heatmap(nb_data, cluster_50, "hpca", "fine", n = 35,
                 filename = "cell_types/cvt_heatmap_hpca_fine")

plot_cvt_heatmap(nb_data, cluster_50, "blueprint", "broad",
                 filename = "cell_types/cvt_heatmap_blueprint_broad")
plot_cvt_heatmap(nb_data, cluster_50, "blueprint", "fine", n = 25,
                 filename = "cell_types/cvt_heatmap_blueprint_fine")

plot_cvt_heatmap(nb_data, cluster_50, "dice", "broad",
                 filename = "cell_types/cvt_heatmap_dice_broad")
plot_cvt_heatmap(nb_data, cluster_50, "dice", "fine",
                 filename = "cell_types/cvt_heatmap_dice_fine")

plot_cvt_heatmap(nb_data, cluster_50, "dmap", "broad",
                 filename = "cell_types/cvt_heatmap_dmap_broad")
plot_cvt_heatmap(nb_data, cluster_50, "dmap", "fine", n = 25,
                 filename = "cell_types/cvt_heatmap_dmap_fine")

plot_cvt_heatmap(nb_data, cluster_50, "monaco", "broad",
                 filename = "cell_types/cvt_heatmap_monaco_broad")
plot_cvt_heatmap(nb_data, cluster_50, "monaco", "fine",
                 filename = "cell_types/cvt_heatmap_monaco_fine")




#' For a given cell type reference dataset, label, and cluster, plot a bar chart
#' that counts the most frequent cell types.
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param selected_cluster Cluster whose cell types should be counted.
#' @param ref Reference cell type dataset.
#' @param label Cell type label.
#' @param max_bars Maximum number of bars in each subplot. Collapse additional
#'   cell types (except "Unknown") into type "Other".
#'
#' @return A ggplot object.
plot_cvt_bar_single <- function(data, clusters, selected_cluster, ref, label,
                                max_bars = 8L) {
  cell_types <-
    data %>%
    preprocess_celltypes(ref, label) %>% 
    select(cell_type, cluster = {{clusters}}) %>% 
    filter(cluster == {{selected_cluster}}) %>%
    pull(cell_type) %>% 
    fct_drop() %>% 
    fct_infreq()
  
  if (nlevels(cell_types) > max_bars) {
    if ("Unknown" %in% levels(cell_types))
      cell_types <- fct_relevel(cell_types, "Unknown", after = 0)

    kept_levels <- levels(cell_types)[1:max_bars - 1]
    cell_types <- fct_relabel(
      cell_types,
      ~case_when(. %in% kept_levels ~ ., TRUE ~ "Other")
    )
  }
  
  p <-
    fct_count(cell_types, prop = TRUE) %>%
    mutate(
      f =
        f %>%
        fct_reorder(n) %>%
        fct_relevel("Other") %>%
        fct_relevel("Unknown")
    ) %>%
    ggplot(aes(f, n)) +
    geom_col(aes(fill = p), show.legend = FALSE) +
    xlab("") +
    ylab("") +
    scale_fill_distiller(
      palette = "YlOrRd",
      direction = 1,
      limits = c(0, 1)
    ) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 270, vjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  p
}

plot_cvt_bar_single(nb_data, cluster_50, "3", "hpca", "broad")



#' Count cell types in a given cluster and plot results in bar charts for all
#' cell type reference datasets.
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param selected_cluster Cluster whose cell types should be counted.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cvt_bar_cluster <- function(data, clusters, selected_cluster,
                                 filename = NULL) {
  references <-
    colnames(data) %>%
    str_match("cell_type_(.+)_broad") %>% 
    magrittr::extract(, 2) %>% 
    magrittr::extract(!is.na(.))
  
  p <- 
    list(ref = references, label = c("broad", "fine")) %>%
    cross_df() %>%
    pmap(
      function(ref, label) {
        plot_cvt_bar_single(
          data,
          clusters = {{clusters}},
          selected_cluster = selected_cluster,
          ref = ref,
          label = label,
          max_bars = if (label == "broad") 8L else 15L
        ) +
          {if (label == "broad") xlab(ref) else NULL} +
          {if (ref == references[length(references)]) ylab(label) else NULL}
      }
    ) %>% 
    wrap_plots(ncol = 2, byrow = FALSE) +
    plot_annotation(
      title = str_glue("Cell types in cluster {selected_cluster}")
    )
  
  ggsave_default(filename, height = 297, width = 210)
  p
}

plot_cvt_bar_cluster(nb_data, cluster_50, "3")

walk(
  levels(nb_data$cluster_50),
  ~plot_cvt_bar_cluster(
    data = nb_data,
    clusters = cluster_50,
    selected_cluster = .,
    filename = str_glue("cell_types/cvt_bar_cluster_{.}")
  )
)



#' For each cluster, plot a bar chart that counts the most frequent cell types.
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param ref Reference cell type dataset.
#' @param label Cell type label.
#' @param max_bars Maximum number of bars in each subplot. Collapse additional
#'   cell types (except "Unknown") into type "Other".
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cvt_bar_dataset <- function(data, clusters, ref, label,
                                 max_bars = 8L, filename = NULL) {
  
  p <-
    levels(pull(data, {{clusters}})) %>%
    map(
      plot_cvt_bar_single,
      data = data, 
      clusters = {{clusters}},
      ref = ref,
      label = label,
      max_bars = max_bars
    ) %>%
    wrap_plots()

  ggsave_default(filename, width = 420, height = 297)
  p
}

colnames(nb_data) %>%
  str_match("cell_type_(.+)_(.+)") %>%
  magrittr::set_colnames(c("match", "ref", "label")) %>% 
  as_tibble() %>% 
  filter(!is.na(ref)) %>% 
  pwalk(
    function(match, ref, label) {
      info("Plotting {ref}, {label}")
      plot_cvt_bar_dataset(
        nb_data,
        cluster_50,
        ref,
        label,
        filename = str_glue("cell_types/cvt_bar_dataset_{ref}_{label}")
      )
    }
  )



#' Count cell types in all subclusters of a given cluster and plot (a, right)
#' results in bar charts for all cell type reference datasets. Also plot (b,
#' upper left) a UMAP with subclusters highlighted, and (c, lower left) up to
#' four UMAPs, each of which highlights cells from a different group.
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param subclusters Column with subcluster IDs.
#' @param selected_cluster Cluster whose cell types should be counted.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cvt_bar_subcluster <- function(data, clusters, subclusters,
                                    selected_cluster,
                                    filename = NULL) {
  subcluster_data <- 
    data %>% 
    filter({{clusters}} == {{selected_cluster}})
  
  ref <-
    colnames(data) %>%
    str_match("cell_type_(.+)_fine") %>% 
    magrittr::extract(, 2) %>% 
    magrittr::extract(!is.na(.))
  
  subcluster_ids <-
    subcluster_data %>% 
    pull({{subclusters}}) %>% 
    unique() %>% 
    str_sort()
  
  p1 <- 
    plot_subclusters(
      subcluster_data,
      umap_1_subcluster,
      umap_2_subcluster,
      {{clusters}},
      {{subclusters}}
    ) +
    theme(strip.text = element_blank())
  
  p2 <- 
    ggplot(subcluster_data, aes(umap_1_subcluster, umap_2_subcluster)) +
    geom_point(
      data = subcluster_data %>% select(!group),
      color = "gray90",
      size = 0.01,
      shape = 20
    ) +
    geom_point(color = "black", size = 0.01, shape = 20) +
    facet_wrap(vars(group), ncol = 1) +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  p3 <-
    list(ref = ref, subcluster = subcluster_ids) %>%
    cross_df() %>%
    pmap(
      function(ref, subcluster) {
        plot_cvt_bar_single(
          subcluster_data,
          clusters = {{subclusters}},
          selected_cluster = subcluster,
          ref = ref,
          label = "fine",
          max_bars = 15L
        ) +
          {if (subcluster == "a") xlab(ref) else NULL} +
          {if (ref == ref[length(ref)]) ylab(subcluster) else NULL}
      }
    ) %>%
    wrap_plots(ncol = length(subcluster_ids), byrow = FALSE)
  
  p <-
    wrap_plots(
      p1, p2, p3,
      widths = c(1, length(subcluster_ids)),
      heights = c(1, 4),
      design = "AC\nBC"
    ) +
    plot_annotation(
      title = str_glue(
        "Cell types in subclusters of cluster {selected_cluster}"
      )
    )

  ggsave_default(filename, height = 297,
                 width = 100 * (length(subcluster_ids) + 1))
  p
}

# plot_cvt_bar_subcluster(nb_data, cluster_50, subcluster_50, "1",
#                         filename = "scvt")

walk(
  levels(nb_data$cluster_50),
  ~plot_cvt_bar_subcluster(
    data = nb_data,
    clusters = cluster_50,
    subclusters = subcluster_20,
    selected_cluster = .,
    filename = str_glue("cell_types/cvt_bar_subcluster_{.}")
  )
)



# Cell type abundance -----------------------------------------------------

nb_data %>%
  group_by(group, sample, cellont_name) %>% 
  summarise(n = n()) %>% 
  mutate(n_rel = n / sum(n)) %>% 
  ungroup() %>% 
  filter(cellont_name == "natural killer cell") %>% 
  mutate(sample = fct_rev(sample)) %>% 
  ggplot(aes(sample, n_rel)) +
  geom_col(aes(fill = group)) +
  coord_flip() +
  ylab("fraction of total cells") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  NULL
ggsave_default("cell_types/abundances_NK", width = 100, height = 80)


#' Determine the enrichment of a cell type in a sample via Fisher's exact test.
#'
#' @param sample Sample.
#' @param cellont_abbr Cell type.
#' @param ref_groups Groups that should be used for comparison.
#'
#' @return A named list with elements `odds_ratio` and `p_val`.
do_fisher_test <- function(sample,
                           cellont_abbr,
                           ref_groups = c("I", "II", "III", "IV")) {
  res <- 
    nb_data %>%
    filter(group %in% ref_groups | sample == .env$sample) %>%
    mutate(
      sample = fct_collapse(
        sample,
        yes = .env$sample,
        other_level = "no"
      ),
      cell_type = fct_collapse(
        cellont_abbr,
        yes = .env$cellont_abbr,
        other_level = "no"
      )
    ) %>%
    count(sample, cell_type) %>%
    pivot_wider(names_from = cell_type, values_from = n) %>%
    column_to_rownames("sample") %>%
    as.matrix() %>%
    fisher.test()
    
  list(
    odds_ratio = res$estimate,
    p_val = res$p.value
  )
}

# do_fisher_test("2016_4503", "B")  # for testing

cell_type_enrichment <-
  nb_data %>%
  distinct(
    group,
    sample = as.character(sample),
    cell_type = as.character(cellont_abbr)
  ) %>%
  filter(group != "I", cell_type != "NB") %>%
  rowwise() %>%
  mutate(fisher = list(do_fisher_test(sample, cell_type, ref_groups = "I"))) %>%
  # filter(cell_type != "NB") %>%
  # rowwise() %>% 
  # mutate(fisher = list(do_fisher_test(sample, cell_type))) %>% 
  unnest_wider(fisher) %>% 
  mutate(p_adj = p.adjust(p_val, method = "BH"))

cell_type_enrichment 



#' Plot a dotplot for cell type enrichment.
#'
#' @param data Results from `do_fisher_test()`.
#' @param p_lim Cap for -log10 p-value.
#' @param or_lim Cap for absolute log2 odds ratio.
#'
#' @return A ggplot object.
plot_enrichment <- function(data, p_lim = 20, or_lim = 2) {
  data %>% 
    mutate(
      log_p_adj = pmin(-log10(p_adj), p_lim),
      log_odds_ratio = log2(odds_ratio),
      log_odds_ratio = case_when(
        log_odds_ratio > or_lim  ~ or_lim,
        log_odds_ratio < -or_lim ~ -or_lim,
        TRUE                     ~ log_odds_ratio
      )
    ) %>% 
    ggplot(aes(cell_type, sample)) + 
    geom_point(aes(color = log_odds_ratio, size = log_p_adj), shape = 16) + 
    scale_color_distiller(
      name = TeX("log_2 odds ratio"),
      palette = "PiYG",
      direction = 1,
      limits = c(-or_lim, or_lim),
      breaks = c(-or_lim, 0, or_lim),
      labels = c(
        str_glue("-{or_lim} or lower"),
        "0",
        str_glue("{or_lim} or higher")
      )
    ) +
    scale_radius(
      name = TeX("-log_{10} p_{adj}"),
      range = c(1, 6),
      labels = function(x) c(x[-length(x)], paste(x[length(x)], "or higher"))
    ) +
    xlab("cell type") +
    facet_grid(vars(group), scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      strip.background = element_blank(),
      strip.text.y.right = element_text(face = "bold", angle = 0)
    ) +
    NULL
}

cell_type_enrichment %>% 
  mutate(
    cell_type =
      as_factor(cell_type) %>%
      fct_relevel(
        levels(nb_data$cellont_cluster) %>%
          str_extract("\\w+") %>% 
          unique()
      )
  ) %>% 
  plot_enrichment(or_lim = 3)
# ggsave_default("cell_types/abundances_enriched", width = 120, height = 130)
ggsave_default("cell_types/abundances_enriched", width = 120, height = 100)



# Unified cell labels -----------------------------------------------------

#' Plot cell types with cluster IDs.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param label_by Column with cluster IDs.
#' @param color_by Column with cell types.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_unified_cell_label <- function(data, x, y, label_by, color_by,
                                    filename = NULL) {
  color_scale <- scale_color_manual(
    name = "cell type",
    values = list(
      "T cell" = "#1f78b4",
      # T8 = "#62a3cb",
      "natural killer cell" = "#a6cee3",
      "B cell" = "#33a02c",
      # B_prec = "#b2df8a",
      "monocyte" = "#ff7f00",
      "neuron" = "#e31a1c",
      "erythroid lineage cell" = "#b15928",
      "hematopoietic precursor cell" = "#6a3d9a",
      # gran = "#cab2d6",
      "plasmacytoid dendritic cell" = "#fdbf6f",
      "other" = "black",
      na = "gray80"
    ),
    guide = guide_legend(override.aes = list(size = 5))
  )
  
  cluster_labels <- 
    data %>% 
    group_by(label = {{label_by}}) %>% 
    summarise({{x}} := mean({{x}}), {{y}} := mean({{y}}))
  
  p <- 
    data %>%
    ggplot(aes({{x}}, {{y}})) +
    geom_point(
      aes(color = {{color_by}}),
      size = .1,
      show.legend = TRUE
    ) +
    geom_text(data = cluster_labels, aes(label = label)) +
    color_scale +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  ggsave_default(filename)
  p
}

plot_unified_cell_label(
  nb_data, umap_1_monocle, umap_2_monocle,
  label_by = cluster_50, color_by = cellont_name,
  filename = "cell_types/unified"
)



# Cluster sizes -----------------------------------------------------------

#' Plot a heatmap of cluster size vs sample (top) and group (bottom).
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param angle_col Angle of column labels.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_size <- function(data, clusters,
                              angle_col = 90, filename = NULL) {
  pivot_and_scale <- function(x) {
    x %>%
      mutate(n_rel = n / sum(n)) %>%
      pivot_wider(
        names_from = "sample",
        values_from = "n_rel",
        id_cols = "cluster"
      ) %>%
      column_to_rownames("cluster") %>%
      as.matrix() %>%
      t()
  }
  
  # scale columnwise for samples and groups
  size_mat <-
    rbind(
      data %>%
        count(group, sample, cluster = {{clusters}}) %>%
        group_by(group, sample) %>%
        pivot_and_scale(),
      data %>%
        count(group, cluster = {{clusters}}) %>%
        group_by(group) %>%
        mutate(sample = group) %>%
        pivot_and_scale()
    )

  annotation_row <-
    bind_rows(
      count(data, group, sample),
      count(data, group)  
    ) %>% 
    mutate(
      split = case_when(
        is.na(sample) ~ "groups", TRUE ~ as.character(group)
      ) %>% 
        as_factor() %>% 
        fct_relevel("groups", after = Inf)
    )

  p <- Heatmap(
    size_mat,
    col = colorRampPalette(brewer.pal(9, "Greens"))(100),
    column_names_rot = angle_col,
    left_annotation = rowAnnotation(
      group = annotation_row$group,
      col = list(
        group = c(
          I = "#1b9e77",
          II = "#d95f02",
          III = "#7570b3",
          IV = "#e7298a"
        )
      ),
      show_legend = FALSE
    ),
    row_split = annotation_row$split,
    cluster_row_slices = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(title = "relative\nabundance")
  )

  ggsave_default(filename, plot = p)
  p
}

plot_cluster_size(nb_data, cluster_50, angle_col = 0,
                  filename = "monocle/cluster_size_vs_samples")

nb_data %>%
  preprocess_celltypes("hpca", "broad", prop = 0.01) %>% 
  plot_cluster_size(cell_type,
                    filename = "cell_types/celltype_hpca_broad_count_vs_samples")

nb_data %>%
  preprocess_celltypes("hpca", "fine", prop = 0.01) %>% 
  plot_cluster_size(cell_type,
                    filename = "cell_types/celltype_hpca_fine_count_vs_samples")



# Putative NB cells -------------------------------------------------------

#' Highlight neurons per group, include a table with clusters that contain
#' neurons.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param cell_types Column with cell types.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_neurons <- function(data, x, y, clusters, cell_types, filename = NULL) {
  p <- wrap_plots(
    data %>%
      ggplot(aes({{x}}, {{y}})) +
      geom_point(
        color = "gray90",
        size = .01,
        shape = 20
      ) +
      geom_point(
        data = data %>%
          filter({{cell_types}} == "Neurons"),
        size = .01,
        shape = 20
      ) +
      coord_fixed() +
      facet_wrap(vars(group)) +
      theme_classic() +
      theme(strip.background = element_blank()) +
      NULL,
    data %>% 
      count(cell_type = {{cell_types}}, {{clusters}}) %>% 
      group_by(cell_type) %>% 
      mutate(n_rel = n / sum(n) * 100) %>% 
      ungroup() %>% 
      filter(cell_type == "Neurons") %>%
      arrange(desc(n)) %>%
      select(!cell_type) %>% 
      tableGrob(rows = NULL)
  ) +
    plot_annotation("Distribution of neurons in clusters")
  
  ggsave_default(filename)
  p
}

plot_neurons(nb_data, umap_1_seurat, umap_2_seurat,
             cluster_0.5, cell_type_hpca_broad,
             filename = "seurat/groupwise_abundance_neurons_0.5")
plot_neurons(nb_data, umap_1_monocle, umap_2_monocle,
             cluster_50, cell_type_hpca_broad,
             filename = "monocle/groupwise_abundance_neurons_50")



# NB gene signatures ------------------------------------------------------

#' Plot gene signature scores, facet by cluster or cell type.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param ncol Number of columns in the plot.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_gene_program <- function(data, x, y, clusters,
                              ncol = 6, filename = NULL) {
  p <- 
    data %>% 
    ggplot(aes({{x}}, {{y}})) +
    geom_vline(xintercept = 0, size = .25) +
    geom_hline(yintercept = 0, size = .25) +
    geom_point(aes(color = factor({{clusters}})), size = .1, show.legend = FALSE) +
    coord_fixed() +
    facet_wrap(vars({{clusters}}), ncol = ncol) +
    theme(panel.grid = element_blank()) +
    NULL
  ggsave_default(filename)
  p
}

plot_gene_program(nb_data,
                  signature_adrenergic, signature_mesenchymal,
                  cluster_50,
                  filename = "gene_programs/gene_programs_am_clusters")
plot_gene_program(nb_data,
                  signature_noradrenergic, signature_ncc_like,
                  cluster_50,
                  filename = "gene_programs/gene_programs_nn_clusters")

nb_data %>% 
  preprocess_celltypes("hpca", "broad", prop = 0.01) %>% 
  plot_gene_program(signature_adrenergic, signature_mesenchymal, cell_type,
                    ncol = 5, filename = "gene_programs/gene_programs_am_hpca")
nb_data %>% 
  preprocess_celltypes("hpca", "broad", prop = 0.01) %>% 
  plot_gene_program( signature_noradrenergic, signature_ncc_like, cell_type,
                  ncol = 4, filename = "gene_programs/gene_programs_nn_hpca")



# Quality control ---------------------------------------------------------

plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, percent_mito,
                  label_direct = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "qc/qc_mtgene_umap")
plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, bcds_score,
                  label_direct = FALSE,
                  color_scale = scale_color_viridis_c(limits = c(0, 1)),
                  filename = "qc/qc_doublets_umap")



# Infiltration rate -------------------------------------------------------

#' Plot tumor infiltration rate as determined by FACS or barcode counting.
#'
#' @param data Metadata.
#' @param filter_col If this column contains ...
#' @param filter_values ... any of these values, the respective barcode
#'                      is regarded as tumor cell.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_infiltration_rate <- function(data, filter_col, filter_values,
                                   filename = NULL) {
  tif_facs <-
    read_csv("metadata/sample_groups.csv", comment = "#") %>%
    filter(!is.na(facs_alive)) %>% 
    mutate(tif_facs = facs_tumor / facs_alive) %>% 
    select(group, sample, tif_facs)
  
  tif_data <-
    data %>%
    mutate(is_tumor = !is.na(match({{filter_col}}, {{filter_values}}))) %>%
    group_by(group, sample) %>% 
    summarise(tif_sc = sum(is_tumor) / n())
  
  infiltration_rates <-
    left_join(
      tif_facs,
      tif_data,
      by = c("group", "sample")
    )

  p <- 
    infiltration_rates %>%
    pivot_longer(
      starts_with("tif"),
      names_to = "method",
      names_prefix = "tif_",
      values_to = "tif"
    ) %>%
    ggplot(aes(method, tif, color = group)) +
    geom_point(show.legend = FALSE) +
    geom_line(aes(group = sample), show.legend = FALSE) +
    geom_text(
      data =
        infiltration_rates %>%
        group_by(group) %>%
        summarise(across(starts_with("tif"), mean)) %>%
        pivot_longer(
          starts_with("tif"),
          names_to = "method",
          names_prefix = "tif_",
          values_to = "tif"
        ) %>%
        mutate(label = sprintf("%.1f", tif * 100), tif = 0.6),
      aes(label = label),
      color = "black"
    ) +
    facet_wrap(vars(group)) +
    ggtitle(
      "Tumor infiltration rates",
      str_glue("{quo_name(enquo(filter_col))} contains ",
               "{{{str_c(filter_values, collapse = ',')}}}")
    ) +
    theme_bw() +
    ylab("tumor infiltration rate") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
  NULL
  
  ggsave_default(filename, height = 160, width = 80)
  p
}

plot_infiltration_rate(nb_data, cell_type_hpca_broad, "Neurons",
                       filename = "monocle/tif_neurons")
plot_infiltration_rate(nb_data, cluster_50, "8",
                       filename = "monocle/tif_clusters")


nb_data %>% 
  mutate(
    cells = case_when(
      cellont_cluster == "NB (8)" ~ "tumor",
      TRUE ~ "other"
    )
  ) %>%
  pivot_longer(
    starts_with("signature"),
    names_to = "signature",
    names_prefix = "signature_"
  ) %>% 
  ggplot(aes(cells, value)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(aes(fill = cells), scale = "width") +
  geom_boxplot(
    aes(group = cells),
    width = 0.075,
    outlier.shape = NA,
    coef = 0,
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_manual(values = c("gray80", "#e41a1c")) +
  coord_flip() +
  facet_grid(vars(signature), switch = "y") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "gray95"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )
ggsave_default("monocle/tif_violin", width = 150, height = 100)



#' Plot a UMAP, highlight cells with the hightest score in a given signature.
#'
#' @param data Metadata.
#' @param signature Column with gene signature scores. 
#' @param top_prop Proportion of cells to highlight.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_signature_highlight <- function(data, signature, top_prop = 0.05,
                                     filename = NULL) {
  p <- 
    data %>%
    ggplot(aes(umap_1_monocle, umap_2_monocle)) +
    geom_point(color = "gray90", size = 0.1) +
    geom_point(
      data =
        nb_data %>% 
        slice_max({{signature}}, prop = top_prop) %>% 
        arrange({{signature}}),
      aes(color = {{signature}}),
      size = 0.1
    ) +
    scale_color_gradient(low = "#0c2c84", high = "#e7298a") +
    coord_fixed() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    NULL
  
  ggsave_default(filename)
  p
}

plot_signature_highlight(nb_data, signature_mesenchymal,
                         filename = "gene_programs/umap_mesenchymal")
plot_signature_highlight(nb_data, signature_ncc_like,
                         filename = "gene_programs/umap_ncc_like")




# Publication figures -----------------------------------------------------

## Figure 1b ----

plot_umap <- function() {
  adjust_positions <- tribble(
    ~label, ~dx, ~dy,
    "13", 1.5, -0.5,
    "14", -1, -1,
    "15", 0, 1,
    "16", 0, 0.5,
    "17", 1, 0,
    "18", -1, 0,
    "19", 1, 0,
    "20", 1, 0,
    "21", 0, -1
  )
  
  cluster_labels <- 
    nb_data %>% 
    group_by(label = cluster_50) %>% 
    summarise(
      umap_1_monocle = mean(umap_1_monocle),
      umap_2_monocle = mean(umap_2_monocle)
    ) %>% 
    left_join(adjust_positions, by = "label") %>%
    mutate(
      label = as_factor(label),
      umap_1_monocle = umap_1_monocle + replace_na(dx, 0),
      umap_2_monocle = umap_2_monocle + replace_na(dy, 0)
    )
   
  p <- 
    nb_data %>%
    ggplot(aes(umap_1_monocle, umap_2_monocle)) +
    geom_point(
      aes(color = cellont_abbr),
      size = .001,
      shape = 16
    ) +
    geom_text(
      data = cluster_labels,
      aes(label = label),
      size = BASE_TEXT_SIZE_MM
    ) +
    scale_x_continuous("UMAP1", breaks = c(-10, 0, 10)) +
    scale_y_continuous("UMAP2", breaks = c(-10, 0, 10)) +
    scale_color_manual(
      name = "cell type",
      values = CELL_TYPE_COLORS,
      guide = guide_legend(override.aes = list(size = 1))
    ) +
    coord_fixed() +
    theme_nb(grid = FALSE) +
    theme(
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = c(.9, .2)
    )
  
  p
}

plot_umap()
ggsave_publication("1b_umap_dataset", type = "png", width = 6, height = 6)


## Figure 1c ----

nb_data %>% 
  mutate(
    cells = case_when(
      cellont_cluster == "NB (8)" ~ "tumor",
      TRUE ~ "other"
    )
  ) %>%
  pivot_longer(
    starts_with("signature"),
    names_to = "signature",
    names_prefix = "signature_"
  ) %>% 
  mutate(signature = recode(signature, ncc_like = "NCC-like")) %>% 
  ggplot(aes(cells, value)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    size = BASE_LINE_SIZE
  ) +
  geom_violin(
    aes(fill = cells),
    scale = "width",
    size = BASE_LINE_SIZE
  ) +
  stat_summary(geom = "point", fun = mean, size = .1) +
  scale_fill_manual(values = c("gray80", CELL_TYPE_COLORS["NB"])) +
  xlab(NULL) +
  ylab("gene signature score") +
  coord_flip() +
  facet_grid(vars(signature)) +
  theme_nb(grid = FALSE) +
  theme(legend.position = "none")

ggsave_publication("1c_gene_programs", width = 5, height = 6)



## Figure 3a ----

plot_celltype_dots <- function(data, p_lim = 20, or_lim = 3) {
  data %>% 
    filter(cell_type != "other") %>% 
    mutate(
      cell_type =
        as_factor(cell_type) %>%
        fct_relevel(
          levels(nb_data$cellont_cluster) %>%
            str_extract("\\w+") %>% 
            unique()
        ),
      log_p_adj = pmin(-log10(p_adj), p_lim),
      log_odds_ratio = log2(odds_ratio),
      log_odds_ratio = case_when(
        log_odds_ratio > or_lim  ~ or_lim,
        log_odds_ratio < -or_lim ~ -or_lim,
        TRUE                     ~ log_odds_ratio
      ),
      sample =
        sample %>% 
        factor(levels = levels(nb_data$sample)) %>% 
        rename_patients() %>% 
        fct_rev(),
      group = rename_groups(group)
    ) %>% 
    ggplot(aes(cell_type, sample)) + 
    geom_point(aes(color = log_odds_ratio, size = log_p_adj), shape = 16) + 
    xlab("cell type") +
    ylab("patient") +
    scale_color_gsea(
      name = TeX("log_2 odds ratio"),
      limits = c(-or_lim, or_lim),
      breaks = c(-or_lim, 0, or_lim),
      labels = c(
        str_glue("-{or_lim} or lower"),
        "0",
        str_glue("{or_lim} or higher")
      ),
      guide = guide_colorbar(ticks = FALSE)
    ) +
    scale_radius(
      name = TeX("-log_{10} p_{adj}"),
      range = c(0.25, 3),
      breaks = c(0, 10, 20),
      labels = function(x) c(x[-length(x)], paste(x[length(x)], "or higher"))
    ) +
    facet_grid(vars(group), scales = "free_y", space = "free_y") +
    theme_nb(grid = FALSE) +
    theme(
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      legend.spacing = unit(5, "mm"),
      legend.margin = margin(0, 0, 0, 0, "mm"),
      panel.spacing = unit(-.5, "pt"),
      strip.text.y.right = element_text(angle = 0)
    )
}

plot_celltype_dots(cell_type_enrichment)
ggsave_publication("3a_cell_type_abundances", width = 5, height = 4)


## Figure 3b ----

nb_data %>%
  group_by(group, sample, cellont_name) %>% 
  summarise(n = n()) %>% 
  mutate(n_rel = n / sum(n)) %>% 
  ungroup() %>% 
  filter(cellont_name == "natural killer cell") %>% 
  mutate(
    sample = rename_patients(sample) %>% fct_rev(),
    group = rename_groups(group)
  ) %>% 
  ggplot(aes(sample, n_rel)) +
  geom_col(aes(fill = group), show.legend = FALSE) +
  annotate(
    "text",
    x = 15,
    y = .21,
    label = "abundance of\nNK cells",
    size = BASE_TEXT_SIZE_MM,
    hjust = 0
  ) +
  xlab("patient") +
  scale_y_continuous(
    "fraction of total cells",
    expand = expansion()
  ) +
  scale_fill_manual(values = GROUP_COLORS) +
  coord_flip() +
  theme_nb(grid = FALSE) +
  theme(
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "mm"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.ontop = TRUE,
    panel.grid.major.x = element_line(color = "white", size = BASE_LINE_SIZE),
  )
ggsave_publication("3b_abundances_NK", width = 5, height = 4)


## Figure S1a ----

set.seed(1L)
nb_data %>%
  select(
    cell, sample,
    umap_1_monocle, umap_2_monocle,
    umap_1_unaligned, umap_2_unaligned
  ) %>% 
  pivot_longer(
    starts_with("umap"),
    names_to = c("coord", "aligned"),
    names_pattern = "(umap_\\d)_(.+)"
  ) %>% 
  pivot_wider(names_from = coord) %>% 
  mutate(
    aligned =
      as_factor(aligned) %>%
      fct_recode(integrated = "monocle", original = "unaligned") %>% 
      fct_rev(),
    sample = rename_patients(sample)
  ) %>% 
  slice_sample(prop = 1) %>% 
  ggplot(aes(umap_1, umap_2)) +
  geom_point(
    aes(color = sample),
    size = .001,
    shape = 16
  ) +
  scale_x_continuous("UMAP1", breaks = c(-10, 0, 10)) +
  scale_y_continuous("UMAP2", breaks = c(-10, 0, 10)) +
  scale_color_manual(
    name = "patient",
    values = PATIENT_COLORS,
    guide = guide_legend(override.aes = list(size = 1), ncol = 2)
  ) +
  coord_fixed() +
  facet_wrap(vars(aligned)) +
  theme_nb(grid = FALSE) +
  theme(
    legend.title.align = 1,
    legend.key.height = unit(1, "mm"),
    legend.key.width = unit(1, "mm"),
    legend.position = c(.92, .24)
  )
ggsave_publication("S1a_umap_integration", type = "png", width = 9, height = 5)


## Figure S1b ----

plot_infiltration_rate <- function(show_mean = FALSE) {
  tif_facs <-
    read_csv("metadata/sample_groups.csv", comment = "#") %>%
    filter(!is.na(facs_alive)) %>% 
    mutate(tif_facs = facs_tumor / facs_alive) %>% 
    select(group, sample, tif_facs)
  
  tif_data <-
    nb_data %>%
    group_by(group, sample) %>% 
    summarise(tif_sc = sum(cluster_50 == "8") / n())
  
  infiltration_rates <-
    left_join(
      tif_facs,
      tif_data,
      by = c("group", "sample")
    ) %>% 
    mutate(group = rename_groups(group) %>% fct_relevel("C", "M", "A", "S"))
  
  if (show_mean) {
    geom_text_mean_tif <- list(
      geom_text(
        data =
          infiltration_rates %>%
          group_by(group) %>%
          summarise(across(starts_with("tif"), mean)) %>%
          pivot_longer(
            starts_with("tif"),
            names_to = "method",
            names_prefix = "tif_",
            values_to = "tif"
          ) %>%
          mutate(
            label = sprintf("%.1f", tif * 100),
            tif = 0.6,
            method = recode(method, facs = "FACS", sc = "scRNA-seq")
          ),
        aes(label = label),
        color = "black",
        size = BASE_TEXT_SIZE_MM
      ),
      geom_text(
        data = tibble(
          label = "mean over\npatients",
          group = factor("C", levels = c("C", "M", "A", "S"))
        ),
        aes(label = label, x = 1.5, y = 0.52),
        color = "black",
        size = BASE_TEXT_SIZE_MM
      )
    )
  } else {
    geom_text_mean_tif <- NULL
  }
  
  p <- 
    infiltration_rates %>%
    pivot_longer(
      starts_with("tif"),
      names_to = "method",
      names_prefix = "tif_",
      values_to = "tif"
    ) %>%
    mutate(method = recode(method, facs = "FACS", sc = "scRNA-seq")) %>% 
    ggplot(aes(method, tif, color = group)) +
    geom_point(show.legend = FALSE) +
    geom_line(aes(group = sample), show.legend = FALSE) +
    geom_text_mean_tif +
    xlab(NULL) +
    scale_y_continuous(
      name = "tumor infiltration rate",
      limits = c(0, 0.6)
    ) +
    scale_color_manual(values = GROUP_COLORS) +
    facet_wrap(vars(group), nrow = 1) +
    theme_nb(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "grey92", size = BASE_LINE_SIZE)
    )
  
  p
}

plot_infiltration_rate()
ggsave_publication("S1b_tif", width = 9, height = 4)
