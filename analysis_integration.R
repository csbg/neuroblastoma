library(pheatmap)
library(RColorBrewer)
library(ggpmisc)
library(patchwork)
library(tidyverse)
library(ggalluvial)
library(egg)
library(fs)


#' Save a plot with default settings.
#'
#' @param filename output filename, stored as PNG in folder ./plots; any missing
#'   folders are created. If NULL, do not create output.
#' @param width figure width in mm (size is not limited)
#' @param height figure height in mm
#' @param crop If `TRUE`, crop the generated plot.
#' @param ... additional arguments passed to ggsave().
#'
#' @return The final file name, invisibly.
ggsave_default <- function(filename, width = 297, height = 210,
                           crop = TRUE, ...) {
  if (is.null(filename))
    return()
  
  filename <- str_glue("plots/{filename}.png")
  filename %>% path_dir() %>% dir_create()
  
  ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
         width = width, height = height, ...)
  
  if (crop) knitr::plot_crop(filename)
  
  invisible(filename)
}



# Load data ---------------------------------------------------------------

nb_data <- readRDS("data_generated/all_datasets_current/nb_metadata.rds")



# Clusters ----------------------------------------------------------------

#' Plot all clusters.
#'
#' @param data Data extracted from a Seurat object.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param label_direct If `TRUE`, print cluster labels at cluster mean
#' @param show_resolution If `TRUE`, print clustering resolution, which is
#'   derived from the cluster column name.
#' @param color_scale A ggplot2 color scale used for coloring the points.
#'   If `NULL`, use the default discrete color scale.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_clusters_all <- function(data, x, y, clusters, label_direct = TRUE,
                              show_resolution = TRUE, color_scale = NULL,
                              filename = NULL) {
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
    {
      if (show_resolution)
        annotate(
          "text_npc", npcx = 0.05, npcy = 1, size = 6,
          label = str_glue("resolution {res}")
        )
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

plot_clusters_all(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.2,
                  filename = "clusters_all_UMAP_0.2")
plot_clusters_all(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.5,
                  filename = "clusters_all_UMAP_0.5")
plot_clusters_all(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.8,
                  filename = "clusters_all_UMAP_0.8")
plot_clusters_all(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.5,
                  filename = "clusters_all_tSNE_0.5")



#' Plot all clusters, facet by sample.
#'
#' @param data Data extracted from a Seurat object.
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

plot_clusters_per_sample(nb_data, UMAP_1, UMAP_2,
                         integrated_snn_res.0.5, sample,
                         nrow = 3, filename = "clusters_sample_UMAP_0.5")



#' Highlight an individual cluster over all cells on the background (light
#' gray), facet by cluster.
#'
#' @param data Data extracted from a Seurat object.
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
}

plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.2,
                        nrow = 3, filename = "clusters_hl_UMAP_0.2")
plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.5,
                        nrow = 4, filename = "clusters_hl_UMAP_0.5")
plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.8,
                        nrow = 5, filename = "clusters_hl_UMAP_0.8")



#' Highlight an individual cluster, generate one plot per cluster.
#'
#' @param data Data extracted from a Seurat object.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param folder Folder where individual plots will be saved (filenames
#'   correspond to cluster IDs).
#'
#' @return `NULL`
plot_clusters_selected <- function(data, x, y, clusters, folder = NULL) {
  res <-
    enquo(clusters) %>% 
    rlang::as_name() %>%
    str_sub(start = -3L)
  
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
          label = str_glue("cluster {cluster}, resolution {res}, n = {n_cells}")
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

plot_clusters_selected(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.2,
                       folder = "clusters_highlighted_UMAP_0.2")
plot_clusters_selected(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.5,
                       folder = "clusters_highlighted_UMAP_0.5")
plot_clusters_selected(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.8,
                       folder = "clusters_highlighted_UMAP_0.8")



# Clustering resolution ---------------------------------------------------

#' Create alluvial plot that shows how clusters change with resolution.
#'
#' @param data Data extracted from a Seurat object.
#' @param other_prop Fractions of low-res clusters that are smaller than this
#'   value will flow to a high-res cluster "other".
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_changes <- function(data, other_prop = 0.1, filename = NULL) {
  cluster_flow <- 
    data %>%
    count(
      low = integrated_snn_res.0.2,
      mid = integrated_snn_res.0.5,
      high = integrated_snn_res.0.8
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

plot_cluster_changes(nb_data, filename = "cluster_change")



# Samples -----------------------------------------------------------------

plot_clusters_all(nb_data, UMAP_1, UMAP_2, sample, label_direct = FALSE,
                  show_resolution = FALSE, filename = "samples_all_UMAP")
plot_clusters_all(nb_data, tSNE_1, tSNE_2, sample, label_direct = FALSE,
                  show_resolution = FALSE, filename = "samples_all_tSNE")



# Cell types --------------------------------------------------------------

plot_clusters_all(nb_data, UMAP_1, UMAP_2, cell_type_broad_lumped,
                  label_direct = FALSE, show_resolution = FALSE,
                  filename = "celltype_broad_all_UMAP")

plot_clusters_per_sample(nb_data, UMAP_1, UMAP_2, cell_type_broad_lumped,
                         sample, nrow = 3, show_legend = TRUE,
                         filename = "celltype_broad_sample_UMAP")

plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, cell_type_broad_lumped,
                        nrow = 3, filename = "celltype_broad_hl_UMAP")


plot_clusters_all(nb_data, UMAP_1, UMAP_2, cell_type_fine_lumped,
                  label_direct = FALSE, show_resolution = FALSE,
                  filename = "celltype_fine_all_UMAP")

plot_clusters_per_sample(nb_data, UMAP_1, UMAP_2, cell_type_fine_lumped,
                         sample, nrow = 3, show_legend = TRUE,
                         filename = "celltype_fine_sample_UMAP")

plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, cell_type_fine_lumped,
                        nrow = 4, filename = "celltype_fine_hl_UMAP")



# Cell types vs clusters --------------------------------------------------

#' Plot a heatmap of cell type vs cluster.
#'
#' @param data Data extracted from a Seurat object.
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
      fontsize = 15
    ) %>%
    {.}
  set_last_plot(p)

  ggsave_default(filename)
  invisible(p)
}


plot_cvt_heatmap(nb_data, cell_type_broad, integrated_snn_res.0.5,
                 filename = "cvt_heatmap_broad_0.5")
plot_cvt_heatmap(nb_data, cell_type_fine, integrated_snn_res.0.5,
                 lump_n = 35, filename = "cvt_heatmap_fine_0.5")

walk(
  levels(nb_data$sample),
  function(sample) {
    message("Plotting sample ", sample)
    plot_cvt_heatmap(nb_data, cell_type_broad, integrated_snn_res.0.5,
                     sample = sample,
                     filename = str_glue("cvt_heatmap_broad/{sample}"))
    plot_cvt_heatmap(nb_data, cell_type, integrated_snn_res.0.5,
                     sample = sample, lump_n = 35,
                     filename = str_glue("cvt_heatmap_fine/{sample}"))
  }
)


#' For each cluster, plot a bar chart that counts the most frequent cell types.
#'
#' @param data Data extracted from a Seurat object.
#' @param cell_types Column with cell types.
#' @param clusters Column with cluster IDs.
#' @param lump_n Preserve the `lump_n` most common cell types, lump the
#'   remaining ones as "other".
#' @param filename Name of output file.
#' @param ... further arguments passed to `ggsave_default()`
#'
#' @return A ggplot object.
plot_cvt_bar <- function(data, cell_types, clusters,
                         lump_n = 10, filename = NULL, ...) {
  pdata <- 
    data %>% 
    transmute(
      cell_type = {{cell_types}},
      cluster = {{clusters}}
    )
  
  plot_single <- function(cluster) {
    message("Plotting ", cluster)
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
    ggsave_default(
      str_glue("{filename}/{cluster}"),
      plot = set_panel_size(
        ps,
        width = unit(30, "mm"),
        height = unit(30, "mm")
      )
    )
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


plot_cvt_bar(nb_data, cell_type_broad, integrated_snn_res.0.2,
             lump_n = 5, filename = "cvt_bar_broad_0.2")
plot_cvt_bar(nb_data, cell_type_fine, integrated_snn_res.0.2,
             lump_n = 8, filename = "cvt_bar_fine_0.2",
             width = 420, height = 297)
plot_cvt_bar(nb_data, cell_type_broad, integrated_snn_res.0.5,
             lump_n = 5, filename = "cvt_bar_broad_0.5")
plot_cvt_bar(nb_data, cell_type_fine, integrated_snn_res.0.5,
             lump_n = 8, filename = "cvt_bar_fine_0.5",
             width = 420, height = 297)
plot_cvt_bar(nb_data, cell_type_broad, integrated_snn_res.0.8,
             lump_n = 5, filename = "cvt_bar_broad_0.8")
plot_cvt_bar(nb_data, cell_type_fine, integrated_snn_res.0.8,
             lump_n = 8, filename = "cvt_bar_fine_0.8",
             width = 600, height = 297)


#' For each group, plot a bar chart that counts the most frequent cell types.
#' Only include selected clusters.
#'
#' @param data Data extracted from a Seurat object.
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

plot_cvt_group(nb_data, cell_type_broad_lumped, integrated_snn_res.0.5,
               selected_clusters = 5,
               filename = "cvt_bar_broad_0.5_groupwise_cluster_5")
plot_cvt_group(nb_data, cell_type_broad_lumped, integrated_snn_res.0.5,
               selected_clusters = 9,
               filename = "cvt_bar_broad_0.5_groupwise_cluster_9")



# Cluster sizes -----------------------------------------------------------

#' Plot a heatmap of cluster size vs sample (top) and group (bottom).
#' 
#' Columns are scaled separately within the top and bottom part.
#'
#' @param data Data extracted from a Seurat object.
#' @param clusters Column with cluster IDs.
#' @param angle_col Angle of pheatmap column labels.
#' @param center_color_scale If `TRUE`, center color scale at zero.
#' @param color_palette Color palette to use.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_size <- function(data, clusters,
                              angle_col = "90", center_color_scale = TRUE,
                              color_palette = brewer.pal(11, "RdYlBu"),
                              filename = NULL) {
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
      t() %>%
      scale()
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
  
  if (center_color_scale) {
    size_mat_max <- max(abs(range(size_mat, na.rm = TRUE)))
    breaks <- seq(from = -size_mat_max, to = size_mat_max, length.out = 101)  
  } else {
    breaks <- NA
  }
  
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

  pheatmap(
    size_mat,
    color = colorRampPalette(color_palette)(100),
    breaks = breaks,
    border_color = "white",
    angle_col = angle_col,
    annotation_row = annotation_row,
    gaps_row = gaps_row,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Relative cluster sizes in samples and groups (scaled columnwise)"
  ) %>%
  set_last_plot()

  ggsave_default(filename)
}

plot_cluster_size(nb_data, integrated_snn_res.0.5, angle_col = "0",
                  center_color_scale = FALSE,
                  color_palette = brewer.pal(9, "Greens"),
                  filename = "cluster_size_vs_samples")

nb_data %>%
  mutate(cell_type = fct_explicit_na(cell_type_broad_lumped)) %>% 
  plot_cluster_size(cell_type, center_color_scale = FALSE,
                    color_palette = brewer.pal(9, "Greens"),
                    filename = "celltype_broad_count_vs_samples")

nb_data %>%
  mutate(cell_type = fct_explicit_na(cell_type_fine_lumped)) %>% 
  plot_cluster_size(cell_type, center_color_scale = FALSE,
                    color_palette = brewer.pal(9, "Greens"),
                    filename = "celltype_fine_count_vs_samples")
  


# Putative NB cells -------------------------------------------------------

nb_data %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(
    color = "gray90",
    size = .01,
    shape = 20
  ) +
  geom_point(
    data = nb_data %>%
      filter(cell_type_broad == "Neurons"),
    size = .01,
    shape = 20
  ) +
  coord_fixed() +
  facet_wrap(vars(group)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  NULL
ggsave_default("groupwise_abundance_neurons")

nb_data %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(
    color = "gray90",
    size = .01,
    shape = 20
  ) +
  geom_point(
    aes(color = refined_cluster),
    data = nb_data %>%
      filter(refined_cluster %in% c("5a", "9a", "5b", "9b", "20b")),
    size = .01,
    shape = 20
  ) +
  coord_fixed() +
  facet_wrap(vars(group)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  NULL
ggsave_default("groupwise_abundance_cluster5_9")



# NB gene signatures ------------------------------------------------------

#' Plot gene signature scores, facet by cluster or cell type.
#'
#' @param data Data extracted from a Seurat object.
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

plot_gene_program(nb_data, adrenergic, mesenchymal, integrated_snn_res.0.5,
                  filename = "gene_programs_am_clusters")
plot_gene_program(nb_data, noradrenergic, ncc_like, integrated_snn_res.0.5,
                  filename = "gene_programs_nn_clusters")
plot_gene_program(nb_data, adrenergic, mesenchymal, cell_type_broad_lumped,
                  ncol = 5, filename = "gene_programs_am_ctb")
plot_gene_program(nb_data, noradrenergic, ncc_like, cell_type_broad_lumped,
                  ncol = 5, filename = "gene_programs_nn_ctb")



#' Plot gene signature scores, facet by group and selected clusters/cell types.
#'
#' @param data Data extracted from a Seurat object.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param selected_clusters IDs of clusters that should be plotted.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_gene_program_cvg <- function(data, x, y, clusters, selected_clusters,
                                  filename = NULL) {
  p <-
    data %>% 
    filter({{clusters}} %in% selected_clusters) %>% 
    ggplot(aes({{x}}, {{y}})) +
    geom_vline(xintercept = 0, size = .25) +
    geom_hline(yintercept = 0, size = .25) +
    geom_point(aes(color = group), size = .1, show.legend = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    coord_fixed() +
    facet_grid(vars({{clusters}}), vars(group)) +
    theme(panel.grid = element_blank()) +
    NULL
  ggsave_default(filename)
  p
}

plot_gene_program_cvg(nb_data, adrenergic, mesenchymal, integrated_snn_res.0.5,
                      c(5, 9, 10, 17),
                      filename = "gene_programs_am_clusters_vs_group")
plot_gene_program_cvg(nb_data, noradrenergic, ncc_like, integrated_snn_res.0.5,
                      c(5, 9, 10, 17),
                      filename = "gene_programs_nn_clusters_vs_group")
plot_gene_program_cvg(nb_data, adrenergic, mesenchymal, cell_type_broad_lumped,
                      c("T_cell", "Neurons"),
                      filename = "gene_programs_am_ctb_vs_groups")
plot_gene_program_cvg(nb_data, noradrenergic, ncc_like, cell_type_broad_lumped,
                      c("T_cell", "Neurons"),
                      filename = "gene_programs_nn_ctb_vs_groups")



# Quality control ---------------------------------------------------------

plot_clusters_all(nb_data, UMAP_1, UMAP_2, percent.mt,
                  label_direct = FALSE, show_resolution = FALSE,
                  color_scale = scale_color_viridis_c(),
                  filename = "qc_mtgene_umap")

ggplot(nb_data, aes(integrated_snn_res.0.5, percent.mt)) +
  geom_violin(aes(fill = integrated_snn_res.0.5), show.legend = FALSE) +
  geom_jitter(alpha = .1) +
  xlab("Cluster") +
  ylab("% mitochondrial genes") +
  theme_classic()
ggsave_default("qc_mtgene_per_cluster", width = 200, height = 150)



# Subclustering -----------------------------------------------------------

#' Plot all sub clusters.
#'
#' @param data Data extracted from a Seurat object.
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
    facet_wrap(vars({{superclusters}})) +
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


#' For each subcluster, plot a bar chart that counts the most frequent cell
#' types. All barch charts that belong to a single supercluster are printed in
#' one row.
#'
#' @param data Data extracted from a Seurat object.
#' @param cell_types Column with cell types.
#' @param superclusters Column with supercluster IDs.
#' @param subclusters Column with subcluster IDs.
#' @param lump_n Preserve the `lump_n` most common cell types, lump the
#'   remaining ones as "other".
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_scvt_bar <- function(data, cell_types, superclusters, subclusters,
                          lump_n = 10, filename = NULL) {
  pdata <- 
    data %>% 
    transmute(
      cell_type = {{cell_types}},
      supercluster = {{superclusters}},
      subcluster = {{subclusters}}
    )
  
  plot_single <- function(supercluster, subcluster) {
    bar_data <- 
      pdata %>% 
      filter(supercluster == {{supercluster}}, subcluster == {{subcluster}})
    
    if (nrow(bar_data) == 0)
      return(plot_spacer())
    
    bar_data %>%   
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
      annotate("text_npc", npcx = 0.1, npcy = 0.9,
               label = str_glue("{supercluster}-{subcluster}")) +
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
  
  p <- 
    list(
      supercluster = as.integer(levels(pdata$supercluster)),
      subcluster = as.integer(levels(pdata$subcluster))
    ) %>% 
    cross_df() %>% 
    arrange(supercluster) %>% 
    pmap(plot_single) %>%
    wrap_plots(nrow = nlevels(pdata$supercluster))
  set_last_plot(p)
  
  ggsave_default(filename, width = 420, height = 840, crop = FALSE)
  p
}


## of clusters ------------------------------------------------------------

plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
                 integrated_snn_res.0.5, subcluster_mid_0.2,
                 filename = "subcluster_mid_0.2")

plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
                 integrated_snn_res.0.5, cell_type_broad_lumped,
                 label_direct = FALSE, filename = "subcluster_mid_ctb")

plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
                 integrated_snn_res.0.5, percent.mt,
                 color_scale = scale_color_viridis_c(),
                 label_direct = FALSE, filename = "subcluster_mid_mtgenes")

plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
                 integrated_snn_res.0.5, nFeature_SCT,
                 color_scale = scale_color_viridis_c(),
                 label_direct = FALSE, filename = "subcluster_mid_nFeature")

plot_subclusters(nb_data, subcluster_mid_UMAP_1, subcluster_mid_UMAP_2,
                 integrated_snn_res.0.5, nCount_SCT,
                 color_scale = scale_color_viridis_c(),
                 label_direct = FALSE, filename = "subcluster_mid_nCount")

plot_scvt_bar(nb_data, cell_type_fine,
              integrated_snn_res.0.5, subcluster_mid_0.2,
              lump_n = 4, filename = "subcluster_mid_bars")



## of cell types ----------------------------------------------------------

plot_subclusters(nb_data, subcluster_ctb_UMAP_1, subcluster_ctb_UMAP_2,
                 cell_type_broad, subcluster_ctb_0.2,
                 filename = "subcluster_ctb_0.2")

plot_subclusters(nb_data, subcluster_ctb_UMAP_1, subcluster_ctb_UMAP_2,
                 cell_type_broad, integrated_snn_res.0.5,
                 label_direct = FALSE, filename = "subcluster_ctb_mid")



# Refined clusters --------------------------------------------------------
  
plot_clusters_all(nb_data, UMAP_1, UMAP_2, refined_cluster,
                  show_resolution = FALSE,
                  filename = "clusters_all_UMAP_refined")

plot_clusters_selected(nb_data, UMAP_1, UMAP_2, refined_cluster,
                       folder = "clusters_highlighted_UMAP_refined")




# Infiltration rate -------------------------------------------------------

#' Plot tumor infiltration rate as determined by FACS or barcode counting.
#'
#' @param data Data extracted from a Seurat object.
#' @param filter_col If this column contains ...
#' @param filter_values ... any of these values, the respective barcode
#'                      is regarded as tumor cell.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_infiltration_rate <- function(data, filter_col, filter_values,
                                   filename = NULL) {
  tif_facs <-
    read_csv("data_raw/metadata/sample_groups.csv", comment = "#") %>%
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
plot_infiltration_rate(nb_data, refined_cluster, c("5a", "9a", "20b"),
                       filename = "tif_clusters")
