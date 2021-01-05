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

load_singler_details <- function(folder) {
  df <- 
    read_csv(str_glue("{folder}/nb_singler_details.csv")) %>%
    rowwise() %>% 
    mutate(median_score = median(c_across(starts_with("scores.")))) %>% 
    ungroup() %>% 
    mutate(
      diff_next = tuning.scores.first - tuning.scores.second,
      label_name = str_glue("scores.{make.names(labels)}")
    ) %>% 
    pivot_longer(
      starts_with("scores."),
      names_to = "scores",
      values_to = "label_score"
    ) %>% 
    filter(scores == label_name) %>% 
    mutate(delta_score = label_score - median_score) %>% 
    select(!c(label_name, scores))
  
  left_join(
    df,
    df %>% 
      group_by(labels) %>% 
      summarise(
        median_delta_score = median(delta_score),
        mad_delta_score = mad(delta_score)
      ),
    by = "labels"
  ) %>% 
  mutate(z_score = (delta_score - median_delta_score) / mad_delta_score)
}

load_data <- function(folder) {
  files <- c(
    "nb_clusters_0.2.csv",
    "nb_clusters_0.5.csv",
    "nb_clusters_0.8.csv",
    "nb_tsne.csv",
    "nb_umap.csv",
    "nb_singler.csv"
  )
  
  nb_groups <-
    read_csv(
      "data_raw/metadata/sample_groups.csv",
      col_types = "cf",
      comment = "#"
    ) %>%
    mutate(group = fct_relevel(group, "I", "II", "III", "IV")) %>%
    arrange(group, sample) %>%
    mutate(sample = as_factor(sample), sample_id = as.integer(sample))
  
  str_glue("{folder}/{files}") %>% 
    map(read_csv) %>% 
    reduce(left_join, by = c("cell", "sample")) %>%
    left_join(nb_groups, by = "sample") %>%
    extract(
      cell_type,
      into = "cell_type_broad",
      regex = "([^:]*)",
      remove = FALSE
    ) %>%
    mutate(
      sample = as_factor(sample) %>% fct_reorder(sample_id),
      integrated_snn_res.0.2 = as_factor(integrated_snn_res.0.2) %>%
        fct_inseq(),
      integrated_snn_res.0.5 = as_factor(integrated_snn_res.0.5) %>%
        fct_inseq(),
      integrated_snn_res.0.8 = as_factor(integrated_snn_res.0.8) %>%
        fct_inseq(),
      cell_type = as_factor(cell_type) %>% fct_infreq(),
      cell_type_broad = as_factor(cell_type_broad) %>% fct_infreq()
    ) %>%
    select(!sample_id)
}

add_filtered_cell_types <- function(df_seurat,
                                    df_singler,
                                    min_z_score = -3,
                                    min_delta_score = -Inf,
                                    min_diff_next = 0) {
  df_singler <- 
    df_singler %>% 
    mutate(
      labels_pruned = case_when(
        z_score >= min_z_score &
          delta_score >= min_delta_score &
          diff_next > min_diff_next
          ~ labels,
        TRUE
          ~ NA_character_
      )
    ) %>% 
    select(cell, sample, cell_type_filtered = labels_pruned)
  
  df_seurat %>% 
    left_join(df_singler, by = c("cell", "sample")) %>% 
    extract(
      cell_type_filtered,
      into = "cell_type_filtered_broad",
      regex = "([^:]*)",
      remove = FALSE
    ) %>% 
    mutate(
      cell_type_filtered = as_factor(cell_type_filtered) %>%
        fct_infreq(),
      cell_type_filtered_broad = as_factor(cell_type_filtered_broad) %>%
        fct_infreq()
    )
}


singler_data <- load_singler_details("data_generated/all_datasets_current")
nb_data <-
  load_data("data_generated/all_datasets_current") %>% 
  add_filtered_cell_types(singler_data)



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
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_clusters_all <- function(data, x, y, clusters, label_direct = TRUE,
                              show_resolution = TRUE, filename = NULL) {
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
    scale_color_hue(guide = guide_legend(override.aes = list(size = 5))) +
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
        ggsave_default(str_glue("{folder}/{cluster}"))
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



# Cell types (SingleR) ----------------------------------------------------

# broad cell types
nb_data_ctb <- 
  nb_data %>% 
  mutate(
    cell_type_filtered_broad = fct_lump_prop(cell_type_filtered_broad, 0.01)
  )

plot_clusters_all(nb_data_ctb, UMAP_1, UMAP_2, cell_type_filtered_broad,
                  label_direct = FALSE, show_resolution = FALSE,
                  filename = "celltype_broad_all_UMAP")

plot_clusters_per_sample(nb_data_ctb, UMAP_1, UMAP_2, cell_type_filtered_broad,
                         sample, nrow = 3, show_legend = TRUE,
                         filename = "celltype_broad_sample_UMAP")

plot_clusters_highlight(nb_data_ctb, UMAP_1, UMAP_2, cell_type_filtered_broad,
                        nrow = 3, filename = "celltype_broad_hl_UMAP")


# fine cell types
nb_data_ctf <- 
  nb_data %>% 
  mutate(cell_type_filtered = fct_lump_n(cell_type_filtered, 24))

plot_clusters_all(nb_data_ctf, UMAP_1, UMAP_2, cell_type_filtered,
                  label_direct = FALSE, show_resolution = FALSE,
                  filename = "celltype_fine_all_UMAP")

plot_clusters_per_sample(nb_data_ctf, UMAP_1, UMAP_2, cell_type_filtered,
                         sample, nrow = 3, show_legend = TRUE,
                         filename = "celltype_fine_sample_UMAP")

plot_clusters_highlight(nb_data_ctf, UMAP_1, UMAP_2, cell_type_filtered,
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


plot_cvt_heatmap(nb_data, cell_type_filtered_broad, integrated_snn_res.0.5,
                 filename = "cvt_heatmap_broad_0.5")
plot_cvt_heatmap(nb_data, cell_type_filtered, integrated_snn_res.0.5,
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


plot_cvt_bar(nb_data, cell_type_filtered_broad, integrated_snn_res.0.2,
             lump_n = 5, filename = "cvt_bar_broad_0.2")
plot_cvt_bar(nb_data, cell_type_filtered, integrated_snn_res.0.2,
             lump_n = 8, filename = "cvt_bar_fine_0.2",
             width = 420, height = 297)
plot_cvt_bar(nb_data, cell_type_filtered_broad, integrated_snn_res.0.5,
             lump_n = 5, filename = "cvt_bar_broad_0.5")
plot_cvt_bar(nb_data, cell_type_filtered, integrated_snn_res.0.5,
             lump_n = 8, filename = "cvt_bar_fine_0.5",
             width = 420, height = 297)
plot_cvt_bar(nb_data, cell_type_filtered_broad, integrated_snn_res.0.8,
             lump_n = 5, filename = "cvt_bar_broad_0.8")
plot_cvt_bar(nb_data, cell_type_filtered, integrated_snn_res.0.8,
             lump_n = 8, filename = "cvt_bar_fine_0.8",
             width = 600, height = 297)




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

nb_data_ctb %>%
  mutate(cell_type = fct_explicit_na(cell_type_filtered_broad)) %>% 
  plot_cluster_size(cell_type, center_color_scale = FALSE,
                    color_palette = brewer.pal(9, "Greens"),
                    filename = "celltype_broad_count_vs_samples")

nb_data_ctf %>%
  mutate(cell_type = fct_explicit_na(cell_type_filtered)) %>% 
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
    data = nb_data %>%
      filter(integrated_snn_res.0.5 %in% c(4, 10)),
    size = .01,
    shape = 20
  ) +
  coord_fixed() +
  facet_wrap(vars(group)) +
  theme_classic() +
  theme(strip.background = element_blank()) +
  NULL
ggsave_default("groupwise_abundance_cluster4_10")
