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
    )
}


singler_data <- load_singler_details("data_generated/all_datasets_sc")
nb_data <-
  load_data("data_generated/all_datasets_sc") %>% 
  add_filtered_cell_types(singler_data)



# Clusters ----------------------------------------------------------------

#' Plot all clusters.
#'
#' @param data Data extracted from a Seurat object.
#' @param x Column with x-axis data.
#' @param y Column with y-axis data.
#' @param clusters Column with cluster IDs.
#' @param label_direct If `TRUE`, print cluster labels at cluster mean
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_clusters_all <- function(data, x, y, clusters, label_direct = TRUE,
                              filename = NULL) {
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
    {if (label_direct) geom_text(data = cluster_labels, aes(label = label))} +
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

plot_clusters_all(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.5,
                  filename = "clusters_all_UMAP_0.5")
plot_clusters_all(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.8,
                  filename = "clusters_all_UMAP_0.8")
# plot_clusters_all(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.5,
#                   filename = "clusters_all_tSNE_0.5")
# plot_clusters_all(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.8,
#                   filename = "clusters_all_tSNE_0.8")



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
    facet_wrap(vars({{sample}}), nrow = nrow) +
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
# plot_clusters_per_sample(nb_data, UMAP_1, UMAP_2,
#                          integrated_snn_res.0.8, sample,
#                          nrow = 3, filename = "clusters_sample_UMAP_0.8")
# plot_clusters_per_sample(nb_data, tSNE_1, tSNE_2,
#                          integrated_snn_res.0.5, sample,
#                          nrow = 3, filename = "clusters_sample_tSNE_0.5")
# plot_clusters_per_sample(nb_data, tSNE_1, tSNE_2,
#                          integrated_snn_res.0.8, sample,
#                          nrow = 3, filename = "clusters_sample_tSNE_0.8")



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

plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.5,
                        nrow = 4, filename = "clusters_hl_UMAP_0.5")
plot_clusters_highlight(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.8,
                        nrow = 5, filename = "clusters_hl_UMAP_0.8")
# plot_clusters_highlight(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.5,
#                         nrow = 4, filename = "clusters_hl_tSNE_0.5")
# plot_clusters_highlight(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.8,
#                         nrow = 5, filename = "clusters_hl_tSNE_0.8")



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

plot_clusters_selected(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.5,
                       folder = "clusters_highlighted_UMAP_0.5")
plot_clusters_selected(nb_data, UMAP_1, UMAP_2, integrated_snn_res.0.8,
                       folder = "clusters_highlighted_UMAP_0.8")
# plot_clusters_selected(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.5,
#                        folder = "clusters_highlighted_tSNE_0.5")
# plot_clusters_selected(nb_data, tSNE_1, tSNE_2, integrated_snn_res.0.8,
#                        folder = "clusters_highlighted_tSNE_0.8")



# Clustering resolution ---------------------------------------------------

#' Create alluvial plot that shows how clusters change upon increasing resolution.
#'
#' @param data Data extracted from a Seurat object.
#' @param clusters Selection of low-resolution clusters to plot.
#' @param other_prop Fractions of low-res clusters that are smaller than this
#'   value will flow to a high-res cluster "other".
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cluster_changes <- function(data, clusters,
                                 other_prop = 0.1, filename = NULL) {
  cluster_flow <- 
    data %>%
    count(low = integrated_snn_res.0.5, high = integrated_snn_res.0.8) %>% 
    mutate(high = as.character(high)) %>%
    group_by(low) %>% 
    mutate(
      prop = n / sum(n),
      high = case_when(
        prop >= other_prop ~ high,
        TRUE ~ "other"
      ) %>% as_factor()
    )
  
  p <- 
    cluster_flow %>% 
    filter(low %in% clusters) %>% 
    ggplot(aes(axis1 = low, axis2 = high, y = n)) +
    stat_alluvium(aes(fill = low), show.legend = FALSE) +
    stat_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(
      "clustering resolution",
      limits = c("low", "high"),
      expand = expansion(add = .2)
    ) +
    scale_y_continuous(
      "cumulative number of cells",
      expand = expansion()
    ) +
    theme_classic() +
    theme(axis.line.x = element_blank()) +
    NULL
  
  ggsave_default(filename, width = 210, height = 297)
  p
}

plot_cluster_changes(nb_data, as.character(0:11),
                     filename = "cluster_change_a")
plot_cluster_changes(nb_data, as.character(12:23),
                     filename = "cluster_change_b")



# Cell types (SingleR) ----------------------------------------------------

# broad cell types
nb_data_ctb <- 
  nb_data %>% 
  mutate(cell_type_broad = fct_lump_prop(cell_type_broad, 0.01))

plot_clusters_all(nb_data_ctb, UMAP_1, UMAP_2, cell_type_broad,
                  label_direct = FALSE, filename = "celltype_broad_all_UMAP")
# plot_clusters_all(nb_data_ctb, tSNE_1, tSNE_2, cell_type_broad,
#                   label_direct = FALSE, filename = "celltype_broad_all_tSNE")

plot_clusters_per_sample(nb_data_ctb, UMAP_1, UMAP_2, cell_type_broad, sample,
                         nrow = 3, show_legend = TRUE,
                         filename = "celltype_broad_sample_UMAP")
# plot_clusters_per_sample(nb_data_ctb, tSNE_1, tSNE_2, cell_type_broad, sample,
#                          nrow = 3, show_legend = TRUE,
#                          filename = "celltype_broad_sample_tSNE")

plot_clusters_highlight(nb_data_ctb, UMAP_1, UMAP_2, cell_type_broad,
                        nrow = 3, filename = "celltype_broad_hl_UMAP")
# plot_clusters_highlight(nb_data_ctb, tSNE_1, tSNE_2, cell_type_broad,
#                         nrow = 3, filename = "celltype_broad_hl_tSNE")

# plot_clusters_selected(nb_data_ctb, UMAP_1, UMAP_2, cell_type_broad,
#                        folder = "cell_type_broad_highlighted_UMAP")
# plot_clusters_selected(nb_data_ctb, tSNE_1, tSNE_2, cell_type_broad,
#                        folder = "cell_type_broad_highlighted_tSNE")

# fine cell types
nb_data_ctf <- 
  nb_data %>% 
  mutate(cell_type = fct_lump_n(cell_type, 24))

plot_clusters_all(nb_data_ctf, UMAP_1, UMAP_2, cell_type,
                  label_direct = FALSE, filename = "celltype_fine_all_UMAP")
# plot_clusters_all(nb_data_ctf, tSNE_1, tSNE_2, cell_type,
#                   label_direct = FALSE, filename = "celltype_fine_all_tSNE")

plot_clusters_per_sample(nb_data_ctf, UMAP_1, UMAP_2, cell_type, sample,
                         nrow = 3, show_legend = TRUE,
                         filename = "celltype_fine_sample_UMAP")
# plot_clusters_per_sample(nb_data_ctf, tSNE_1, tSNE_2, cell_type, sample,
#                          nrow = 3, show_legend = TRUE,
#                          filename = "celltype_fine_sample_tSNE")

plot_clusters_highlight(nb_data_ctf, UMAP_1, UMAP_2, cell_type,
                        nrow = 4, filename = "celltype_fine_hl_UMAP")
# plot_clusters_highlight(nb_data_ctf, tSNE_1, tSNE_2, cell_type,
#                         nrow = 4, filename = "celltype_fine_hl_tSNE")

# plot_clusters_selected(nb_data_ctf, UMAP_1, UMAP_2, cell_type,
#                        folder = "cell_type_fine_highlighted_UMAP")
# plot_clusters_selected(nb_data_ctf, tSNE_1, tSNE_2, cell_type,
#                        folder = "cell_type_fine_highlighted_tSNE")



# Cell types vs clusters --------------------------------------------------

#' Plot a heatmap of cell type vs cluster.
#'
#' @param data Data extracted from a Seurat object.
#' @param cell_types Column with cell types.
#' @param clusters Column with cluster IDs.
#' @param lump_n Preserve the `lump_n` most common cell types, lump the
#'   remaining ones as "other".
#' @param sample Only include this sample.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cvt_heatmap <- function(data, cell_types, clusters,
                             lump_n = NULL, sample = NULL, filename = NULL) {
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
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      angle_col = "0"
    ) %>%
    {.}
  set_last_plot(p)

  ggsave_default(filename)
  invisible(p)
}


plot_cvt_heatmap(nb_data, cell_type_broad, integrated_snn_res.0.5,
                 filename = "cvt_heatmap_broad_0.5")
plot_cvt_heatmap(nb_data, cell_type, integrated_snn_res.0.5,
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
      scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 100)) +
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
      plot = set_panel_size(ps, width = unit(30, "mm"), height = unit(30, "mm"))
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

plot_cvt_bar(nb_data, cell_type_broad, integrated_snn_res.0.5,
             lump_n = 5, filename = "cvt_bar_broad_0.5")
plot_cvt_bar(nb_data, cell_type, integrated_snn_res.0.5,
             lump_n = 8, filename = "cvt_bar_fine_0.5",
             width = 420, height = 297)
plot_cvt_bar(nb_data, cell_type_broad, integrated_snn_res.0.8,
             lump_n = 5, filename = "cvt_bar_broad_0.8")
plot_cvt_bar(nb_data, cell_type, integrated_snn_res.0.8,
             lump_n = 8, filename = "cvt_bar_fine_0.8",
             width = 600, height = 297)
