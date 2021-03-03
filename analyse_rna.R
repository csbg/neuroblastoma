library(pheatmap)
library(RColorBrewer)
library(ggpmisc)
library(patchwork)
library(tidyverse)
library(ggalluvial)
library(egg)
library(fs)
library(clustree)
library(bluster)
library(ggbeeswarm)
source("common_functions.R")



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
  
  res <-
    enquo(clusters) %>% 
    rlang::as_name() %>%
    str_sub(start = -3L)
  
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

plot_clusters_all(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.2,
                  filename = "seurat/clusters_all_umap_0.2")
plot_clusters_all(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                  filename = "seurat/clusters_all_umap_0.5")
plot_clusters_all(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.8,
                  filename = "seurat/clusters_all_umap_0.8")
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

plot_clusters_per_sample(nb_data, umap_1_seurat, umap_2_seurat,
                         cluster_0.5, sample,
                         nrow = 3, filename = "seurat/clusters_sample_umap_0.5")

plot_clusters_per_sample(nb_data %>% mutate(group2 = group),
                         umap_1_seurat, umap_2_seurat,
                         cluster_0.5, group2,
                         nrow = 2, filename = "seurat/clusters_groups_umap_0.5")

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

plot_clusters_highlight(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                        nrow = 4,
                        filename = "seurat/clusters_highlight_umap_0.5")
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

plot_clusters_selected(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                       folder = "seurat/clusters_selected_umap_0.5")
plot_clusters_selected(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
                       folder = "monocle/clusters_selected_umap_50")



# Cluster correspondence --------------------------------------------------

#' Create alluvial plot that shows how clusters change with resolution.
#'
#' @param data Metadata.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_corr_resolution <- function(data, filename = NULL) {
  cluster_flow <- 
    data %>%
    count(
      low = cluster_0.2,
      mid = cluster_0.5,
      high = cluster_0.8
    ) %>% 
    mutate(
      low = fct_relabel(low, ~str_glue("L-{.}")),
      mid = fct_relabel(mid, ~str_glue("M-{.}")),
      high = fct_relabel(high, ~str_glue("H-{.}"))
    )
  
  p <- 
    cluster_flow %>% 
    ggplot(aes(axis1 = low, axis2 = mid, axis3 = high,  y = n)) +
    stat_alluvium(aes(fill = mid), show.legend = FALSE, reverse = FALSE) +
    stat_stratum(reverse = FALSE) +
    stat_stratum(
      geom = "text",
      aes(label = str_sub(after_stat(stratum), start = 3L)),
      reverse = FALSE,
      size = 3
    ) +
    scale_x_discrete(
      "clustering resolution",
      limits = c("low", "mid", "high"),
      expand = expansion(add = .2)
    ) +
    scale_y_continuous(
      "cumulative number of cells",
      expand = expansion()
    ) +
    coord_flip() +
    theme_classic() +
    theme(axis.line.x = element_blank()) +
    NULL
  
  ggsave_default(filename, width = 400, height = 200)
  p
}

plot_cluster_corr_resolution(
  nb_data,
  filename = "seurat/cluster_correspondence_resolution"
)


# alternative visualization using clustree
nb_data %>% 
  select(starts_with("cluster_0")) %>% 
  clustree(prefix = "cluster_", edge_arrow = FALSE)
ggsave_default("seurat/clustree")

nb_data %>% 
  select(cluster_20, cluster_50) %>% 
  clustree(prefix = "cluster_", edge_arrow = FALSE)
ggsave_default("monocle/clustree")

nb_data %>% 
  select(cluster_0.5, cluster_50) %>% 
  clustree(prefix = "cluster_", edge_arrow = FALSE)
ggsave_default("clustree_seurat_monocle")



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

plot_cluster_diagnostics(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.2,
                         filename = "seurat/cluster_diagnostics_0.2")

plot_cluster_diagnostics(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                         filename = "seurat/cluster_diagnostics_0.5")

plot_cluster_diagnostics(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.8,
                         filename = "seurat/cluster_diagnostics_0.8")



# Samples -----------------------------------------------------------------

plot_clusters_all(nb_data %>% arrange(umap_1_seurat),
                  umap_1_seurat, umap_2_seurat, sample,
                  label_direct = FALSE, filename = "seurat/samples_all_umap")
plot_clusters_all(nb_data %>% arrange(umap_1_monocle),
                  umap_1_monocle, umap_2_monocle, sample,
                  label_direct = FALSE, filename = "monocle/samples_all_umap")



# Cell types --------------------------------------------------------------

plot_clusters_all(nb_data,
                  umap_1_seurat, umap_2_seurat,
                  cell_type_broad_lumped,
                  label_direct = FALSE,
                  filename = "seurat/celltype_all_umap")

plot_clusters_per_sample(nb_data,
                         umap_1_seurat, umap_2_seurat,
                         cell_type_broad_lumped,
                         sample,
                         nrow = 3, show_legend = TRUE,
                         filename = "seurat/celltype_sample_umap")

plot_clusters_per_sample(nb_data %>% mutate(group2 = group),
                         umap_1_seurat, umap_2_seurat,
                         cell_type_broad_lumped,
                         group2,
                         nrow = 2, show_legend = TRUE,
                         filename = "seurat/celltype_group_umap")

plot_clusters_highlight(nb_data,
                        umap_1_seurat, umap_2_seurat,
                        cell_type_broad_lumped,
                        nrow = 3, filename = "seurat/celltype_highlight_umap")



plot_clusters_all(nb_data,
                  umap_1_monocle, umap_2_monocle,
                  cell_type_broad_lumped,
                  label_direct = FALSE,
                  color_scale = scale_color_manual(
                    values = c(
                      T_cell = "#1f78b4",
                      NK_cell = "#a6cee3",
                      B_cell = "#33a02c",
                      "Pro-B_cell_CD34+" = "#b2df8a",
                      Monocyte = "#ff7f00",
                      Neurons = "#e31a1c",
                      Erythroblast = "#b15928",
                      "Pre-B_call_CD34-" = "#6a3d9a",
                      GMP = "#6a3d9a",
                      CMP = "#6a3d9a",
                      "Pro-Myelocyte" = "#6a3d9a",
                      Other = "#6a3d9a",
                      "NA" = "gray80"
                    ),
                    na.value = "gray80"
                  ),
                  filename = "monocle/celltype_all_umap")

plot_clusters_per_sample(nb_data,
                         umap_1_monocle, umap_2_monocle,
                         cell_type_broad_lumped,
                         sample,
                         nrow = 3, show_legend = TRUE,
                         filename = "monocle/celltype_sample_umap")

plot_clusters_per_sample(nb_data %>% mutate(group2 = group),
                         umap_1_monocle, umap_2_monocle,
                         cell_type_broad_lumped,
                         group2,
                         nrow = 2, show_legend = TRUE,
                         filename = "monocle/celltype_group_umap")

plot_clusters_highlight(nb_data,
                        umap_1_monocle, umap_2_monocle,
                        cell_type_broad_lumped,
                        nrow = 3, filename = "monocle/celltype_highlight_umap")



# Cell types vs clusters --------------------------------------------------

#' Plot a heatmap of cell type vs cluster.
#'
#' @param data Metadata.
#' @param cell_types Column with cell types.
#' @param clusters Column with cluster IDs.
#' @param lump_n Preserve the `lump_n` most common cell types, lump the
#'   remaining ones as "other".
#' @param sample Only include this sample.
#' @param cluster_cols If `TRUE`, cluster the columns (i.e., Seurat clusters).
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cvt_heatmap <- function(data, cell_types, clusters,
                             lump_n = NULL, sample = NULL, cluster_cols = TRUE,
                             filename = NULL) {
  if (is.null(sample))
    sample <- "."
  if (is.null(lump_n))
    lump_n <- Inf
  
  p <-
    data %>% 
    mutate(
      cell_type = fct_lump_n({{cell_types}}, lump_n),
      cluster = {{clusters}}
    ) %>%
    filter(sample %>% str_detect({{sample}})) %>%     
    count(cluster, cell_type, .drop = FALSE) %>% 
    group_by(cluster) %>%
    mutate(n_rel = n / sum(n)) %>%
    select(!n) %>%
    ungroup() %>%
    pivot_wider(names_from = "cell_type", values_from = "n_rel") %>%
    column_to_rownames("cluster") %>%
    as.matrix() %>%
    t() %>%
    replace_na(0) %>%
    pheatmap(
      color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
      legend = FALSE,
      border_color = "white",
      cluster_cols = cluster_cols,
      cluster_rows = FALSE,
      angle_col = "0",
      fontsize = 12
    )
  set_last_plot(p)

  ggsave_default(filename)
  invisible(p)
}


plot_cvt_heatmap(nb_data, cell_type_broad, cluster_0.5,
                 filename = "seurat/cvt_heatmap_broad_0.5")
plot_cvt_heatmap(nb_data, cell_type_fine, cluster_0.5,
                 lump_n = 35, filename = "seurat/cvt_heatmap_fine_0.5")

plot_cvt_heatmap(nb_data, cell_type_broad, cluster_50,
                 filename = "monocle/cvt_heatmap_broad_50")
plot_cvt_heatmap(nb_data, cell_type_fine, cluster_50,
                 lump_n = 35, filename = "monocle/cvt_heatmap_fine_50")



#' For each cluster, plot a bar chart that counts the most frequent cell types.
#'
#' @param data Metadata.
#' @param cell_types Column with cell types.
#' @param clusters Column with cluster IDs.
#' @param lump_n Preserve the `lump_n` most common cell types, lump the
#'   remaining ones as "other".
#' @param save_subplots If `TRUE`, also save each subplot into a separate file.
#' @param filename Name of output file.
#' @param ... further arguments passed to `ggsave_default()`
#'
#' @return A ggplot object.
plot_cvt_bar <- function(data, cell_types, clusters,
                         lump_n = 10, save_subplots = FALSE,
                         filename = NULL, ...) {
  pdata <- 
    data %>% 
    transmute(
      cell_type = {{cell_types}},
      cluster = {{clusters}}
    )
  
  plot_single <- function(cluster) {
    ps <- 
      pdata %>% 
      filter(cluster == {{cluster}}) %>%
      mutate(cell_type = fct_lump_n(cell_type, lump_n)) %>%
      count(cell_type) %>%
      mutate(
        cell_type = cell_type %>% 
          fct_reorder(n) %>%
          fct_relevel("Other"),
        n_rel = n / sum(n) * 100
      ) %>% 
      ggplot(aes(cell_type, n)) +
      geom_col(aes(fill = n_rel), show.legend = FALSE) +
      annotate("text_npc", npcx = 0.5, npcy = 0.9, label = cluster) +
      xlab("") +
      ylab("") +
      scale_fill_distiller(
        palette = "YlOrRd",
        direction = 1,
        limits = c(0, 100)
      ) +
      coord_flip() +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 270, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
      ) +
      NULL
    if (save_subplots) {
      message("Plotting ", cluster)
      ggsave_default(
        str_glue("{filename}/{cluster}"),
        plot = set_panel_size(
          ps,
          width = unit(30, "mm"),
          height = unit(30, "mm")
        )
      ) 
    }
    ps
  }
  
  p <-
    levels(pdata$cluster) %>% 
    map(plot_single) %>% 
    wrap_plots()
  set_last_plot(p)
  
  ggsave_default(filename, ...)
  p
}


plot_cvt_bar(nb_data, cell_type_broad, cluster_0.5,
             lump_n = 5, filename = "seurat/cvt_bar_broad_0.5")
plot_cvt_bar(nb_data, cell_type_fine, cluster_0.5,
             lump_n = 8, filename = "seurat/cvt_bar_fine_0.5",
             width = 420, height = 297)

plot_cvt_bar(nb_data, cell_type_broad, cluster_50,
             lump_n = 5, width = 420, height = 297,
             filename = "monocle/cvt_bar_broad_cluster_50")
plot_cvt_bar(nb_data, cell_type_broad, partition_50,
             lump_n = 5,
             filename = "monocle/cvt_bar_broad_partition_50")



#' For each group, plot a bar chart that counts the most frequent cell types.
#' Only include selected clusters.
#'
#' @param data Metadata.
#' @param cell_types Column with cell types.
#' @param clusters Column with cluster IDs.
#' @param selected_clusters Vector of selected cluster IDs.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cvt_group <- function(data, cell_types, clusters, selected_clusters,
                           filename = NULL) {
  p <- 
    map(
      levels(data$group),
      function(group) {
        data_plot <- 
          data %>%
          filter(
            {{clusters}} %in% {{selected_clusters}},
            group == {{group}}
          ) %>%
          count(cell_type = {{cell_types}}) %>%
          mutate(
            cell_type = fct_reorder(cell_type, n),
            n_rel = n / sum(n) * 100
          )
        
        ggplot(data_plot, aes(cell_type, n)) +
          geom_col(aes(fill = n_rel), show.legend = FALSE) +
          annotate(
            "text_npc",
            npcx = 0.5,
            npcy = 0.9,
            label = str_glue("group {group}\n{sum(data_plot$n)} cells")
          ) +
          xlab("") +
          ylab("") +
          scale_fill_distiller(
            palette = "YlOrRd",
            direction = 1,
            limits = c(0, 100)
          ) +
          coord_flip() +
          theme_classic() +
          theme(
            axis.text.x = element_text(angle = 270, vjust = 0.5),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold")
          ) +
          NULL
      }
    ) %>% 
    wrap_plots() +
    plot_annotation(
      str_glue("Cell types in cluster ",
               "{str_c(selected_clusters, collapse = ', ')}, ",
               "split by patient groups")
    )
  
  ggsave_default(filename, width = 200, height = 150)
  p
}



# Cluster sizes -----------------------------------------------------------

#' Plot a heatmap of cluster size vs sample (top) and group (bottom).
#'
#' @param data Metadata.
#' @param clusters Column with cluster IDs.
#' @param angle_col Angle of pheatmap column labels.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_size <- function(data, clusters,
                              angle_col = "90", filename = NULL) {
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
    data %>%
    distinct(group, sample) %>%
    column_to_rownames("sample")

  gaps_row <-
    data %>%
    group_by(group) %>%
    summarise(n = n_distinct(sample)) %>%
    pull(n) %>%
    cumsum()

  p <- pheatmap(
    size_mat,
    color = colorRampPalette(brewer.pal(9, "Greens"))(100),
    border_color = "white",
    angle_col = angle_col,
    annotation_row = annotation_row,
    gaps_row = gaps_row,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Relative cluster sizes in samples and groups"
  )
  set_last_plot(p)

  ggsave_default(filename)
  invisible(p)
}

plot_cluster_size(nb_data, cluster_0.5,
                  angle_col = "0",
                  filename = "seurat/cluster_size_vs_samples")

plot_cluster_size(nb_data, cluster_50,
                  angle_col = "0",
                  filename = "monocle/cluster_size_vs_samples")

nb_data %>%
  mutate(cell_type = fct_explicit_na(cell_type_broad_lumped)) %>% 
  plot_cluster_size(cell_type, 
                    filename = "celltype_broad_count_vs_samples")

nb_data %>%
  mutate(cell_type = fct_explicit_na(cell_type_fine_lumped)) %>% 
  plot_cluster_size(cell_type,
                    filename = "celltype_fine_count_vs_samples")



# Putative NB cells -------------------------------------------------------

#' Highlight neurons per group, include a table with clusters that contain
#' neurons.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_neurons <- function(data, x, y, clusters, filename = NULL) {
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
          filter(cell_type_broad == "Neurons"),
        size = .01,
        shape = 20
      ) +
      coord_fixed() +
      facet_wrap(vars(group)) +
      theme_classic() +
      theme(strip.background = element_blank()) +
      NULL,
    data %>% 
      count(cell_type = cell_type_broad, {{clusters}}) %>% 
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

plot_neurons(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
             filename = "seurat/groupwise_abundance_neurons_0.5")
plot_neurons(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.8,
             filename = "seurat/groupwise_abundance_neurons_0.8")
plot_neurons(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
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
                  cluster_0.5,
                  filename = "seurat/gene_programs_am_clusters")
plot_gene_program(nb_data,
                  signature_noradrenergic, signature_ncc_like,
                  cluster_0.5,
                  filename = "seurat/gene_programs_nn_clusters")

plot_gene_program(nb_data,
                  signature_adrenergic, signature_mesenchymal,
                  cluster_50,
                  filename = "monocle/gene_programs_am_clusters")
plot_gene_program(nb_data,
                  signature_noradrenergic, signature_ncc_like,
                  cluster_50,
                  filename = "monocle/gene_programs_nn_clusters")

plot_gene_program(nb_data,
                  signature_adrenergic, signature_mesenchymal,
                  cell_type_broad_lumped,
                  ncol = 5, filename = "gene_programs_am_ctb")
plot_gene_program(nb_data,
                  signature_noradrenergic, signature_ncc_like,
                  cell_type_broad_lumped,
                  ncol = 4, filename = "gene_programs_nn_ctb")



# Quality control ---------------------------------------------------------

## Mitochondrial genes ----

plot_clusters_all(nb_data, umap_1_seurat, umap_2_seurat, percent_mt,
                  label_direct = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "seurat/qc_mtgene_umap")
plot_clusters_all(nb_data, umap_1_monocle, umap_2_monocle, percent_mt,
                  label_direct = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "monocle/qc_mtgene_umap")



## Doublets ----

#' Plot doublet scores calculated by scds.
#' 
#' For each type of score, plot (1) a UMAP with points colored by score,
#' (2) a histogram of scores, and (3) a boxplot of score vs cluster.
#'
#' @param data Metadata.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_doublet_scores <- function(data, x, y, clusters, filename = NULL) {
  subplot_histogram <- function(doublet_score) {
    ggplot(nb_data, aes({{doublet_score}})) +
      geom_histogram(bins = 100)
  }
  
  subplot_boxplot <- function(doublet_score) {
    data %>% 
      mutate(
        cluster = fct_reorder(
          {{clusters}}, {{doublet_score}}, .desc = TRUE
        )
      ) %>% 
      ggplot(aes(cluster, {{doublet_score}})) +
      geom_boxplot(
        aes(fill = {{clusters}}),
        outlier.shape = NA,
        show.legend = FALSE
      ) +
      NULL
  }
  
  p <- wrap_plots(
    plot_clusters_all(data %>% arrange(cxds_score),
                      {{x}}, {{y}}, cxds_score,
                      label_direct = FALSE,
                      color_scale = scale_color_viridis_c()),
    plot_clusters_all(data %>% arrange(bcds_score),
                      {{x}}, {{y}}, bcds_score,
                      label_direct = FALSE,
                      color_scale = scale_color_viridis_c()),
    plot_clusters_all(data %>% arrange(hybrid_score),
                      {{x}}, {{y}}, hybrid_score,
                      label_direct = FALSE,
                      color_scale = scale_color_viridis_c()),
    subplot_histogram(cxds_score),
    subplot_histogram(bcds_score),
    subplot_histogram(hybrid_score),
    subplot_boxplot(cxds_score),
    subplot_boxplot(bcds_score),
    subplot_boxplot(hybrid_score),
    nrow = 3
  )

  ggsave_default(filename, width = 420, height = 300)
  p
}

plot_doublet_scores(nb_data, umap_1_seurat, umap_2_seurat, cluster_0.5,
                    filename = "seurat/qc_doublet_scores")
plot_doublet_scores(nb_data, umap_1_monocle, umap_2_monocle, cluster_50,
                    filename = "monocle/qc_doublet_scores")



#' TODO: currently, subclustering is not done
#' # Subclustering -----------------------------------------------------------
#' 
#' #' Plot all sub clusters.
#' #'
#' #' @param data Metadata.
#' #' @param x Column with x-axis data.
#' #' @param y Column with y-axis data.
#' #' @param superclusters Column with supercluster IDs.
#' #' @param subclusters Column with subcluster IDs.
#' #' @param label_direct If `TRUE`, print cluster labels at cluster mean.
#' #' @param color_scale A ggplot2 color scale used for coloring the points.
#' #'   If `NULL`, use the default discrete color scale.
#' #' @param filename Name of output file.
#' #'
#' #' @return A ggplot object.
#' plot_subclusters <- function(data, x, y, superclusters, subclusters,
#'                              label_direct = TRUE, color_scale = NULL,
#'                              filename = NULL) {
#'   cluster_labels <-
#'     data %>%
#'     group_by({{superclusters}}, label = {{subclusters}}) %>%
#'     summarise({{x}} := mean({{x}}), {{y}} := mean({{y}}))
#' 
#'   color_scale <-
#'     color_scale %||%
#'     scale_color_hue(guide = guide_legend(override.aes = list(size = 5)))
#' 
#'   p <-
#'     data %>%
#'     ggplot(aes({{x}}, {{y}})) +
#'     geom_point(
#'       aes(color = {{subclusters}}),
#'       size = .01,
#'       shape = 20,
#'       show.legend = !label_direct
#'     ) +
#'     {
#'       if (label_direct)
#'         geom_text(data = cluster_labels, aes(label = label), size = 3)
#'     } +
#'     color_scale +
#'     facet_wrap(vars({{superclusters}})) +
#'     coord_fixed() +
#'     theme_classic() +
#'     theme(
#'       strip.background = element_blank(),
#'       strip.text = element_text(face = "bold")
#'     ) +
#'     NULL
#' 
#'   ggsave_default(filename, width = 420, height = 297)
#'   p
#' }
#' 
#' 
#' #' For each subcluster, plot a bar chart that counts the most frequent cell
#' #' types. All barch charts that belong to a single supercluster are printed in
#' #' one row.
#' #'
#' #' @param data Metadata.
#' #' @param cell_types Column with cell types.
#' #' @param superclusters Column with supercluster IDs.
#' #' @param subclusters Column with subcluster IDs.
#' #' @param lump_n Preserve the `lump_n` most common cell types, lump the
#' #'   remaining ones as "other".
#' #' @param filename Name of output file.
#' #'
#' #' @return A ggplot object.
#' plot_scvt_bar <- function(data, cell_types, superclusters, subclusters,
#'                           lump_n = 10, filename = NULL) {
#'   pdata <-
#'     data %>%
#'     transmute(
#'       cell_type = {{cell_types}},
#'       supercluster = {{superclusters}},
#'       subcluster = {{subclusters}}
#'     )
#' 
#'   plot_single <- function(supercluster, subcluster) {
#'     bar_data <-
#'       pdata %>%
#'       filter(supercluster == {{supercluster}}, subcluster == {{subcluster}})
#' 
#'     if (nrow(bar_data) == 0)
#'       return(plot_spacer())
#' 
#'     bar_data %>%
#'       mutate(cell_type = fct_lump_n(cell_type, lump_n)) %>%
#'       count(cell_type) %>%
#'       mutate(
#'         cell_type = cell_type %>%
#'           fct_reorder(n) %>%
#'           fct_relevel("Other"),
#'         n_rel = n / sum(n) * 100
#'       ) %>%
#'       ggplot(aes(cell_type, n)) +
#'       geom_col(aes(fill = n_rel), show.legend = FALSE) +
#'       annotate("text_npc", npcx = 0.1, npcy = 0.9,
#'                label = str_glue("{supercluster}-{subcluster}")) +
#'       xlab("") +
#'       ylab("") +
#'       scale_fill_distiller(
#'         palette = "YlOrRd",
#'         direction = 1,
#'         limits = c(0, 100)
#'       ) +
#'       coord_flip() +
#'       theme_classic() +
#'       theme(
#'         axis.text.x = element_text(angle = 270, vjust = 0.5),
#'         strip.background = element_blank(),
#'         strip.text = element_text(face = "bold")
#'       ) +
#'       NULL
#'   }
#' 
#'   p <-
#'     list(
#'       supercluster = as.integer(levels(pdata$supercluster)),
#'       subcluster = as.integer(levels(pdata$subcluster))
#'     ) %>%
#'     cross_df() %>%
#'     arrange(supercluster) %>%
#'     pmap(plot_single) %>%
#'     wrap_plots(nrow = nlevels(pdata$supercluster))
#'   set_last_plot(p)
#' 
#'   ggsave_default(filename, width = 420, height = 840, crop = FALSE)
#'   p
#' }
#' 
#' 
#' ## of clusters ------------------------------------------------------------
#' 
#' plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
#'                  integrated_snn_res.0.5, subcluster_mid_0.2,
#'                  filename = "subcluster_mid_0.2")
#' 
#' plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
#'                  integrated_snn_res.0.5, cell_type_broad_lumped,
#'                  label_direct = FALSE, filename = "subcluster_mid_ctb")
#' 
#' plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
#'                  integrated_snn_res.0.5, percent.mt,
#'                  color_scale = scale_color_viridis_c(),
#'                  label_direct = FALSE, filename = "subcluster_mid_mtgenes")
#' 
#' plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
#'                  integrated_snn_res.0.5, nFeature_SCT,
#'                  color_scale = scale_color_viridis_c(),
#'                  label_direct = FALSE, filename = "subcluster_mid_nFeature")
#' 
#' plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
#'                  integrated_snn_res.0.5, nCount_SCT,
#'                  color_scale = scale_color_viridis_c(),
#'                  label_direct = FALSE, filename = "subcluster_mid_nCount")
#' 
#' plot_scvt_bar(nb_data, cell_type_fine,
#'               integrated_snn_res.0.5, subcluster_mid_0.2,
#'               lump_n = 4, filename = "subcluster_mid_bars")
#' 
#' 
#' 
#' ## of cell types ----------------------------------------------------------
#' 
#' plot_subclusters(nb_data, subcluster_ctb_UMAP_1, subcluster_ctb_UMAP_2,
#'                  cell_type_broad, subcluster_ctb_0.2,
#'                  filename = "subcluster_ctb_0.2")
#' 
#' plot_subclusters(nb_data, subcluster_ctb_UMAP_1, subcluster_ctb_UMAP_2,
#'                  cell_type_broad, integrated_snn_res.0.5,
#'                  label_direct = FALSE, filename = "subcluster_ctb_mid")
#' 
#' 
#' 
#' # Refined clusters --------------------------------------------------------
#' 
#' plot_clusters_all(nb_data, UMAP_1, UMAP_2, refined_cluster,
#'                   show_resolution = FALSE,
#'                   filename = "clusters_all_UMAP_refined")
#' 
#' plot_clusters_selected(nb_data, UMAP_1, UMAP_2, refined_cluster,
#'                        folder = "clusters_highlighted_UMAP_refined")
#' 
#' 
#' 
#' # Cluster splitting -------------------------------------------------------
#' 
#' plot_splitcluster_bar <- function(data, limit = 10, filename = NULL) {
#'   plot_data <- 
#'     data %>% 
#'     group_by(integrated_snn_res.0.5, cell_type_broad_lumped) %>% 
#'     summarise(count = n()) %>%
#'     mutate(prop = count / sum(count) * 100, is_selected = prop > limit)
#'   
#'   p <-
#'     ggplot(plot_data, aes(cell_type_broad_lumped, prop)) +
#'     geom_hline(yintercept = limit, linetype = "dashed", size = 0.25) +
#'     geom_col(aes(fill = is_selected), show.legend = FALSE) +
#'     xlab(NULL) +
#'     ylab("Relative abundance") +
#'     scale_fill_manual(values = c("TRUE" = "#006d2c", "FALSE" = "#a1d99b")) +
#'     coord_flip() +
#'     facet_wrap(vars(integrated_snn_res.0.5), nrow = 4) +
#'     labs(caption = str_glue("Split clusters by cell type ",
#'                             "with abundance >{limit}% ",
#'                             "→ {sum(plot_data$is_selected)} subclusters")) +
#'     theme(panel.grid = element_blank()) +
#'     NULL
#'   ggsave_default(filename)
#'   p
#' }
#' 
#' plot_splitcluster_bar(nb_data, limit = 10)
#' plot_splitcluster_bar(nb_data, limit = 20, filename = "splitcluster_celltypes")
#' 
#' 
#' nb_data %>% 
#'   inner_join(
#'     nb_data %>% 
#'       group_by(integrated_snn_res.0.5, cell_type_broad_lumped) %>% 
#'       summarise(count = n()) %>%
#'       filter(count / sum(count) > 0.2) %>% 
#'       mutate(
#'         split_cluster = str_c(
#'           integrated_snn_res.0.5,
#'           case_when(
#'             n() > 1 ~ letters[row_number()],
#'             TRUE ~ ""
#'           )
#'         )
#'       ) %>% 
#'       select(!count)
#'   ) %>% 
#'   plot_clusters_all(UMAP_1, UMAP_2, integrated_snn_res.0.5,
#'                     show_resolution = FALSE,
#'                     filename = "clusters_all_UMAP_0.5_split")



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

plot_infiltration_rate(nb_data, cell_type_broad, "Neurons",
                       filename = "tif_neurons")
plot_infiltration_rate(nb_data, cluster_0.5, c("6", "10"),
                       filename = "seurat/tif_clusters")
plot_infiltration_rate(nb_data, cluster_50, "8",
                       filename = "monocle/tif_clusters")
