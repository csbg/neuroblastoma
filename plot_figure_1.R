# @DEPI infercnv_output

library(scuttle)
library(monocle3)
library(GenomeInfoDb)
library(infercnv)
library(RColorBrewer)
library(scico)
library(ComplexHeatmap)
library(tidyverse)
library(fs)
source("common_functions.R")
source("styling.R")

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(1, "mm"),
  show_parent_dend_line = FALSE,
  message = FALSE
)



# Loader functions --------------------------------------------------------

#' Load several infercnv result files that are required for plotting.
#'
#' @param folder Folder where infercnv stored results.
#'
#' @return A list with elements `final_obj` and `regions_data`.
read_infercnv_data <- function(folder) {
  final_obj <- readRDS(path_join(c(folder, "run.final.infercnv_obj")))
  regions_data <-
    dir_ls(folder, glob = "*Pnorm*regions*") %>% 
    map_dfr(read_tsv, .id = "file") %>% 
    extract(file, into = "prob", regex = "Pnorm_(.+)\\.pred", convert = TRUE)
    
  list(
    final_obj = final_obj,
    regions_data = regions_data
  )
}



#' Load SNP array CNV data.
#'
#' @param folder Folder with SNP array data, i.e., BEDGRAPH and BED files with
#'   matching names.
#' @param new_sample_names Optional named vector of new sample names, passed to
#'   `dplyr::recode()`.
#'
#' @return A list with elements `logrr_data`, `cnv_regions`, and
#'   `chromosome_size`. `logrr_data` is a data frame with columns sample,
#'   chr(omosome), start, end, logrr (raw log R ratio), and smoothed (smoothed
#'   ratio). `cnv_regions` is a data frame with columns sample, chr, start, end,
#'   type, and copy_number. `chromosome_size` is a dataframe with columns chr
#'   and end.
read_snparray_data <- function(folder, new_sample_names = NULL) {
  # chromosome sizes from GenomeInfoDb
  chromosome_size <- 
    Seqinfo(genome = "GRCh38.p13") %>%
    seqlengths() %>% 
    enframe("chr", "end") %>% 
    filter(str_detect(chr, "^(\\d\\d?|X|Y)$")) %>%
    mutate(chr = as_factor(chr) %>% fct_inseq())
  
  # manually curated data
  cnv_regions <-
    dir_ls(folder, glob = "*.bed") %>% 
    map_dfr(
      read_tsv,
      col_names = FALSE,
      col_types = cols(X2 = "i", X3 = "i", X13 = "d", .default = "c"),
      .id = "sample"
    ) %>%
    transmute(
      sample = sample %>% path_file() %>% path_ext_remove(),
      chr = X1 %>% str_sub(4) %>% factor(levels = chromosome_size$chr),
      start = X2,
      end = X3,
      type = X4,
      copy_number = X13,
      ploidy = replace_na(X21, "2") %>% as.numeric(),
      delta_copy_number = copy_number - ploidy
    )
  
  # raw SNP array signal
  logrr_data <-
    dir_ls(folder, glob = "*.bedgraph") %>% 
    map_dfr(
      read_tsv,
      col_names = c("chr", "start", "end", "logrr", "smoothed"),
      col_types = "ciidd",
      .id = "sample"
    ) %>%
    mutate(
      sample =
        sample %>%
        path_file() %>%
        path_ext_remove() %>%
        path_ext_remove(),
      chr = factor(chr, levels = chromosome_size$chr)
    )
  
  # recode sample names
  if (!is.null(new_sample_names)) {
    cnv_regions <- 
      cnv_regions %>% 
      mutate(sample = recode(sample, !!!new_sample_names))
    logrr_data <- 
      logrr_data %>% 
      mutate(sample = recode(sample, !!!new_sample_names))
  }
  
  # check for missing data
  samples1 <- unique(cnv_regions$sample)
  samples2 <- unique(logrr_data$sample)
  common_samples <- intersect(samples1, samples2)
  all_samples <- union(samples1, samples2)
  if (length(common_samples) != length(all_samples))
    warn("Data missing for some samples")
  info("Loaded samples: {str_c(common_samples, collapse = ', ')}")
    
  list(
    logrr_data = logrr_data,
    cnv_regions = cnv_regions,
    chromosome_size = chromosome_size
  )
}


# Load data ---------------------------------------------------------------

nb_metadata <- readRDS("data_generated/metadata.rds")

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")
colData(nb) <-
  nb_metadata %>%
  mutate(Size_Factor = colData(nb)$Size_Factor) %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

selected_markers <- read_csv("metadata/cell_markers.csv", comment = "#")

infercnv_data <- read_infercnv_data("data_generated/infercnv_output")

snparray_sample_names <- 
  read_csv("metadata/sample_groups.csv", comment = "#") %>% 
  select(snp_array_id, sample) %>% 
  filter(!is.na(snp_array_id)) %>% 
  deframe()
snparray_data <- read_snparray_data("data_raw/snp_array", snparray_sample_names)



# Figures -----------------------------------------------------------------

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(1, "mm")
)


## 1a ----

# no coded needed



## 1b ----

plot_umap <- function() {
  adjust_positions <- tribble(
    ~label, ~dx, ~dy,
    "13", 1.5, -0.5,
    "14", -1, -1,
    "15", 0, 1,
    "16", 0, 0.5,
    "17", 0, 1,
    "18", -1.2, 0,
    "19", 1, 0,
    "20", 1, 0,
    "21", -.5, 1
  )
  
  cluster_labels <- 
    nb_metadata %>% 
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
    nb_metadata %>%
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
      legend.position = c(.9, .22)
    )
  
  p
}

plot_umap()
ggsave_publication("1b_umap_dataset", type = "png", width = 5, height = 5)



## 1c ----

plot_celltype_heatmap <- function(clusters = 1:21, body_width = 150,
                                  collapse = TRUE) {
  # generate matrix of cell type abundances
  make_matrix <- function(ref) {
    cell_type_column <- rlang::sym(str_glue("cell_type_{ref}_broad"))
    
    nb_metadata %>% 
      filter(cluster_50 %in% {{clusters}}) %>% 
      mutate(
        cell_type =
          as_factor(!!cell_type_column) %>%
          fct_infreq() %>% 
          fct_explicit_na("Unknown") %>% 
          fct_relabel(~str_c(ref, .x, sep = "_"))
      ) %>% 
      count(cluster = cluster_50, cell_type) %>%
      group_by(cluster) %>%
      mutate(n_rel = n / sum(n)) %>%
      select(!n) %>%
      ungroup() %>%
      arrange(cell_type) %>% 
      pivot_wider(names_from = "cell_type", values_from = "n_rel") %>%
      arrange(cluster) %>% 
      column_to_rownames("cluster") %>%
      as.matrix() %>%
      replace_na(0)
  }
  
  mat <- 
    map(
      c("blueprint", "hpca", "dice", "dmap", "monaco"),
      make_matrix
    ) %>% 
    reduce(cbind)
  
  # set up column metadata
  col_metadata <-
    tibble(colname = colnames(mat)) %>% 
    left_join(
      read_csv("metadata/celldex_celltypes.csv", comment = "#"),
      by = "colname"
    ) %>% 
    separate(
      colname,
      into = c("ref", "cell_type"),
      extra = "merge",
      remove = FALSE
    ) %>% 
    mutate(
      ref =
        as_factor(ref) %>% 
        fct_recode(
          "Human Primary Cell Atlas" = "hpca",
          "Blueprint/ENCODE" = "blueprint",
          "DICE" = "dice",
          "Novershtern" = "dmap",
          "Monaco" = "monaco"
        ),
      abbr = factor(abbr, levels = names(CELL_TYPE_COLORS))
    ) %>% 
    group_by(ref) %>% 
    arrange(abbr, .by_group = TRUE) %>% 
    ungroup()
  
  mat <- mat[, col_metadata$colname]
  
  if (collapse) {
    mat <-
      mat %>% 
      t() %>% 
      as_tibble(rownames = "colname") %>% 
      left_join(col_metadata, by = "colname") %>% 
      group_by(ref, abbr) %>% 
      summarise(across(where(is.numeric), sum)) %>%
      ungroup() %>%
      unite(ref, abbr, col = "ref_abbr") %>% 
      column_to_rownames("ref_abbr") %>%
      as.matrix() %>%
      t()
    
    col_metadata <- 
      col_metadata %>% 
      distinct(ref, abbr) %>% 
      mutate(cell_type = abbr)
  }
  
  colnames(mat) <- col_metadata$cell_type
  
  # set up row metadata
  row_metadata <- 
    tibble(cluster = levels(nb_metadata$cellont_cluster)) %>% 
    extract(
      cluster,
      into = c("cell_type", "cluster"),
      regex = "(\\w+) \\((\\d+)",
      convert = TRUE
    ) %>%
    filter(cluster %in% {{clusters}}) %>% 
    mutate(
      cell_type =
        cell_type %>% 
        factor(levels = names(CELL_TYPE_COLORS)) %>% 
        fct_drop()
    ) %>%
    arrange(cluster) %>%
    pull(cell_type)
  
  # draw heatmap  
  ht_opt(
    simple_anno_size = unit(1.5, "mm"),
    COLUMN_ANNO_PADDING = unit(1, "pt"),
    DENDROGRAM_PADDING = unit(1, "pt"),
    HEATMAP_LEGEND_PADDING = unit(1, "mm"),
    ROW_ANNO_PADDING = unit(1, "pt"),
    TITLE_PADDING = unit(1, "mm")
  )
  
  set.seed(2)
  Heatmap(
    mat,
    col = colorRampPalette(brewer.pal(9, "YlOrBr"))(100),
    
    heatmap_legend_param = list(
      at = c(0, 1),
      border = FALSE,
      grid_width = unit(2, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      legend_height = unit(15, "mm"),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    name = "relative\nabundance",
    
    row_title = "cluster and assigned cell type",
    row_title_side = "right",
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    row_dend_width = unit(3, "mm"),
    row_km = 5,
    
    column_split = col_metadata$ref,
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    cluster_columns = FALSE,
    show_parent_dend_line = FALSE,
    
    width = unit(body_width, "mm"),
    height = unit((body_width - 4) / ncol(mat) * nrow(mat) + 4, "mm"),
    
    right_annotation = rowAnnotation(
      cell_type = row_metadata,
      col = list(cell_type = CELL_TYPE_COLORS),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      annotation_legend_param = list(
        cell_type = list(
          title = "cell type",
          grid_height = unit(2, "mm"),
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    ),
    
    bottom_annotation = HeatmapAnnotation(
      cell_type = col_metadata$abbr,
      col = list(cell_type = CELL_TYPE_COLORS),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      annotation_legend_param = list(
        cell_type = list(
          title = "cell type",
          grid_height = unit(2, "mm"),
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    )
  ) %>%
    draw(
      gap = unit(50, "mm"),
      column_title = "cell type in reference dataset",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    )
}

(p <- plot_celltype_heatmap(body_width = 80))
ggsave_publication("1c_cell_type_classification",
                   plot = p, width = 11, height = 6)



## 1d ----

plot_cm_dots <- function(counts, features, groups,
                         min_exp = -2.5, max_exp = 2.5,
                         panel_annotation = NULL) {
  scale_and_limit <- function(x) {
    scale(x)[,1] %>% 
      pmax(min_exp) %>% 
      pmin(max_exp)
  }
  
  known_features <- intersect(features, rownames(counts))
  missing_features <- setdiff(features, rownames(counts))
  if (length(missing_features) > 0)
    warn(
      "The following requested features are missing: ",
      "{str_c(missing_features, collapse = ', ')}"
    )
  
  vis_data <- 
    counts[known_features, , drop = FALSE] %>% 
    Matrix::t() %>% 
    as.matrix() %>% 
    as_tibble(rownames = "cell") %>% 
    group_by(id = groups) %>% 
    summarise(
      across(
        where(is.numeric),
        list(
          avg_exp = ~mean(expm1(.)),
          pct_exp = ~length(.[. > 0]) / length(.) * 100
        ),
        .names = "{.col}__{.fn}"
      )
    ) %>% 
    mutate(across(ends_with("avg_exp"), scale_and_limit)) %>%
    pivot_longer(
      !id,
      names_to = c("feature", ".value"),
      names_pattern = "(.+)__(.+)"
    ) %>% 
    mutate(feature = factor(feature, levels = features))
  
  if (!is.null(panel_annotation)) {
    panel_bg <- list(
      geom_point(color = "white"),  # initialize discrete coordinate system
      geom_rect(
        data = panel_annotation,
        aes(xmin = xmin, xmax = xmax,
            ymin = 0.5, ymax = nlevels(vis_data$feature) + 0.5,
            fill = fill),
        show.legend = FALSE,
        inherit.aes = FALSE,
      ),
      scale_fill_identity()
    )
  } else {
    panel_bg <- NULL
  }
  
  color_limits <- c(min(vis_data$avg_exp), max(vis_data$avg_exp))
  
  ggplot(vis_data, aes(id, feature)) +
    panel_bg +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_x_discrete(
      name = "cluster / assigned cell type",
      labels = function(x) str_replace(x, "(.+) \\((.+)\\)", "\\2\n\\1"),
      expand = expansion(add = 0.5)
    ) +
    scale_y_discrete(NULL, expand = expansion(add = 0.5)) +
    scale_color_dotplot(
      "scaled average expression",
      guide = guide_colorbar(
        barheight = unit(2, "mm"),
        barwidth = unit(15, "mm"),
        label.position = "top",
        title.vjust = 0.1,
        ticks = FALSE
      ),
      breaks = color_limits,
      labels = function(x) round(x, 1)
    ) +
    scale_radius("% expressed", range = c(0, 2.5)) +
    theme_nb(grid = FALSE) +
    theme(panel.grid = element_blank()) +
    theme(
      legend.box.just = "bottom",
      legend.key.height = unit(1, "mm"),
      legend.key.width = unit(1, "mm"),
      legend.position = "bottom",
      legend.spacing = unit(0, "mm"),
      legend.margin = margin(-2, 1, 0, 1, "mm"),
      panel.border = element_rect(color = "black", size = .25)
    )
}


plot_canonical_markers <- function() {
  # remove other cells
  nb <- nb[, colData(nb)$cellont_abbr != "other"]
  colData(nb)$cellont_abbr <- fct_drop(colData(nb)$cellont_abbr)
  colData(nb)$cellont_cluster <- fct_drop(colData(nb)$cellont_cluster)
  
  y_annotation_data <-
    selected_markers %>%
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
  
  n_col <- nlevels(colData(nb)$cellont_cluster)
  
  p <-
    plot_cm_dots(
      logcounts(nb),
      rev(selected_markers$gene),
      colData(nb)$cellont_cluster,
      panel_annotation = x_annotation_data
    ) +
    geom_hline(
      yintercept = y_annotation_data$yintercept,
      linetype = "dashed",
      size = 0.25
    ) +
    geom_text(
      data = y_annotation_data,
      aes(
        x = n_col + 1,
        y = label_y,
        label = label
      ),
      size = BASE_TEXT_SIZE_MM,
      hjust = 0
    ) +
    coord_fixed(
      xlim = c(1, n_col),
      clip = "off"
    )
  p
}

plot_canonical_markers()
ggsave_publication("1d_markers", height = 12, width = 7)



## 1e ----

plot_cnv_data_comparison <- function(selected_samples,
                                     probability = 0.5,
                                     logrr_prop = 0.01) {
  # labels and colors for copy numbers
  cn_metadata <-  # colorbrewer 11-class BrBG 
    tribble(
      ~delta_copy_number, ~label,        ~color,
      -2, "complete loss",               "#01665EFF",
      -1, "loss of one copy",            "#80CDC1FF",
      0, "neutral",                      "#F7F7F7FF",
      1, "gain of one copy",             "#DFC27DFF",
      2, "gain of two copies",           "#BF812DFF",
      3, "gain of >2 copies",            "#8C510AFF",
      38, "amplification",               "#543005FF"
    ) %>%
    mutate(label = as_factor(label))
  cn_color <- circlize::colorRamp2(
    cn_metadata$delta_copy_number,
    cn_metadata$color
  )
  
  # prepare data for plotting
  plot_data_regions <-
    infercnv_data$regions_data %>% 
    filter(near(prob, probability)) %>% 
    extract(
      cell_group_name,
      into = c("sample", "group"),
      regex = "malignant_(.*?)_([IV]+)"
    ) %>% 
    transmute(
      sample = sample,
      type = "sc",
      chr = str_sub(chr, 4) %>% factor(levels = snparray_data$chromosome_size$chr),
      start = as.integer(start),
      end = as.integer(end),
      delta_copy_number = as.numeric(state) - 3,
      ymin = 0,
      ymax = 1
    ) %>% 
    bind_rows(
      snparray_data$cnv_regions %>% 
        select(!type:ploidy) %>%
        mutate(type = "snp", ymin = 1, ymax = 2)
    ) %>% 
    mutate(fill = cn_color(delta_copy_number)) %>% 
    filter(chr != "Y")
  # return(plot_data_regions)
  
  plot_data_logrr <-
    snparray_data$logrr_data %>% 
    filter(chr != "Y") %>% 
    group_by(sample) %>% 
    slice_sample(prop = logrr_prop)
  # return(plot_data_logrr)
  
  # filter for selected or common samples
  plot_data_regions <-
    plot_data_regions %>%
    filter(sample %in% selected_samples) %>% 
    mutate(
      sample = sample %>% factor(level = selected_samples) %>% rename_patients()
    )
  plot_data_logrr <-
    plot_data_logrr %>%
    filter(sample %in% selected_samples) %>% 
    mutate(
      sample = sample %>% factor(level = selected_samples) %>% rename_patients()
    )
  
  # make plot
  p <-
    ggplot(NULL, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
    geom_rect(
      data = snparray_data$chromosome_size %>% filter(chr != "Y"),
      aes(xmin = 0, ymin = -0.1, ymax = 0)
    ) +
    geom_rect(
      data = plot_data_regions,
      aes(fill = fill),
      key_glyph = "rect"
    ) +
    geom_line(
      data = plot_data_logrr,
      inherit.aes = FALSE,
      aes(x = start, y = smoothed + 1.5, color = "log R ratio"),
      size = BASE_LINE_SIZE,
      alpha = .5,
      key_glyph = "timeseries"
    ) +
    geom_hline(yintercept = 1, size = BASE_LINE_SIZE) +
    scale_x_continuous(
      name = "chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = NULL,
      breaks = c(0.5, 1.5),
      labels = c("scRNA-seq", "SNP array"),
      expand = c(0, 0),
      sec.axis = sec_axis(
        name = "log R ratio",
        trans = ~. - 1.5,
        breaks = c(-0.5, 0, .5)
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = "black"
    ) +
    scale_fill_identity(
      name = NULL,
      guide = guide_legend(ncol = 1),
      breaks = cn_metadata$color,
      labels = cn_metadata$label
    ) +
    coord_cartesian(ylim = c(0, 2)) +
    facet_grid(
      cols = vars(chr),
      rows = vars(sample),
      space = "free_x",
      scales = "free_x",
      switch = "both"
    ) +
    theme_nb(grid = FALSE) +
    theme(
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.length.x = unit(0, "mm"),
      axis.text.x = element_blank(),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      legend.margin = margin(0, 0, -3, -1, "mm"),
      panel.spacing.x = unit(-.5, "pt"),
      panel.spacing.y = unit(1, "mm"),
      panel.background = element_rect(color = NA, fill = cn_color(0)),
      strip.placement = "outside",
      strip.text.x = element_text(margin = margin()),
      strip.text.y.left = element_text(angle = 0, margin = margin())
    )
  
  p
}

plot_cnv_data_comparison(c("2016_4503", "2005_1702", "2006_2684"))
ggsave_publication("1e_cnv_comparison", width = 18, height = 6)



## S1a ----

plot_umap_unintegrated <- function() {
  set.seed(1L)
  nb_metadata %>%
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
}

plot_umap_unintegrated()
ggsave_publication("S1a_umap_integration", type = "png", width = 9, height = 5)



## S1b ----

plot_infiltration_rate <- function(show_mean = FALSE) {
  tif_facs <-
    read_csv("metadata/sample_groups.csv", comment = "#") %>%
    filter(!is.na(facs_alive)) %>% 
    mutate(tif_facs = facs_tumor / facs_alive) %>% 
    select(group, sample, tif_facs)
  
  tif_data <-
    nb_metadata %>%
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
    ggplot(aes(tif, method, color = group)) +
    geom_point(show.legend = FALSE, size = .5) +
    geom_line(
      aes(group = sample),
      size = BASE_LINE_SIZE,
      show.legend = FALSE
    ) +
    geom_text_mean_tif +
    ylab(NULL) +
    scale_x_continuous(
      name = "tumor infiltration rate",
      limits = c(0, 0.6),
      expand = expansion(mult = c(0.03, 0.01))
    ) +
    scale_color_manual(values = GROUP_COLORS) +
    # coord_flip() +
    facet_grid(vars(group)) +
    theme_nb(grid = FALSE) +
    theme(
      panel.grid.major.x = element_line(
        color = "grey92",
        size = BASE_LINE_SIZE
      ),
      panel.border = element_blank(),
      panel.spacing = unit(1, "mm"),
      strip.text.y = element_text(angle = 0)
    )
  
  p
}

plot_infiltration_rate()
ggsave_publication("S1b_tif", width = 5, height = 4)



## S1c ----

plot_resexp_tumor <- function(cells_per_sample = 50L) {
  # cell metadata
  cell_metadata <- 
    infercnv_data$final_obj@tumor_subclusters$subclusters %>%
    enframe() %>%
    filter(str_starts(name, "malignant")) %>%
    unnest_longer(value) %>%
    select(value) %>%
    unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>%
    left_join(nb_metadata %>% select(cell, sample, group), by = "cell") %>%
    mutate(
      sample = rename_patients(sample),
      group = rename_groups(group)
    )
  
  # subset n cells per sample
  set.seed(1L)
  cell_metadata <- 
    cell_metadata %>% 
    group_by(sample) %>% 
    slice_sample(n = cells_per_sample) %>% 
    ungroup()
  
  # column splits, i.e., vector of chromosme names
  col_split <-
    infercnv_data$final_obj@gene_order$chr %>%
    str_sub(4) %>%
    as_factor()
  
  # matrix with residual expression, rows ordered, Y chromosome removed
  cnv_mat <-
    infercnv_data$final_obj@expr.data %>% 
    magrittr::extract(col_split != "Y", cell_metadata$cell_index) %>% 
    t()
  col_split <- col_split[col_split != "Y"]
  
  # draw heatmap
  Heatmap(
    cnv_mat,
    name = "residual\nexpression",
    col = circlize::colorRamp2(
      breaks = seq(0.85, 1.15, length.out = 7),
      colors = rev(brewer.pal(7, "RdBu"))
    ),
    heatmap_legend_param = list(
      at = c(0.85, 1.15),
      border = FALSE,
      grid_width = unit(1.5, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    row_split = cell_metadata$sample,
    row_gap = unit(0, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    
    cluster_columns = FALSE,
    column_split = col_split,
    cluster_column_slices = FALSE,
    column_gap = unit(0, "mm"),
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    border = TRUE,
    border_gp = gpar(lwd = 0.5),
    
    left_annotation = rowAnnotation(
      group = cell_metadata$group,
      col = list(group = GROUP_COLORS[c("M", "A", "S")]),
      annotation_name_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      annotation_name_side = "bottom",
      annotation_legend_param = list(
        group = list(
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    ),
    
    use_raster = FALSE,
  ) %>%
    draw(
      gap = unit(50, "mm"),
      row_title = "patient",
      row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      column_title = "chromosome",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      padding = unit(c(0, 0, 0, 0), "mm")
    ) 
}

(p <- plot_resexp_tumor())
ggsave_publication("S1c_resexp_tumor", plot = p,
                   type = "png", width = 11, height = 6)



## S1d ----

plot_resexp_marrow <- function(cells_per_type = 50L) {
  # cell metadata
  cell_metadata <- 
    infercnv_data$final_obj@tumor_subclusters$subclusters %>%
    enframe(name = "row_split") %>% 
    filter(!str_starts(row_split, "malignant")) %>% 
    unnest_longer(value) %>%
    select(value) %>%
    unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>%
    left_join(
      nb_metadata %>%
        select(cell, sample, group, cellont_abbr), 
      by = "cell"
    ) %>%
    filter(cellont_abbr != "other") %>% 
    mutate(
      sample = rename_patients(sample),
      group = rename_groups(group)
    )
  
  # subset n cells per sample
  set.seed(1L)
  cell_metadata <- 
    cell_metadata %>% 
    group_by(cellont_abbr) %>% 
    slice_sample(n = cells_per_type) %>% 
    ungroup()
  
  # column splits, i.e., vector of chromosme names
  col_split <-
    infercnv_data$final_obj@gene_order$chr %>%
    str_sub(4) %>%
    as_factor()
  
  # matrix with residual expression, rows ordered, Y chromosome removed
  cnv_mat <-
    infercnv_data$final_obj@expr.data %>% 
    magrittr::extract(col_split != "Y", cell_metadata$cell_index) %>% 
    t()
  col_split <- col_split[col_split != "Y"]
  
  # draw heatmap
  Heatmap(
    cnv_mat,
    name = "residual\nexpression",
    col = circlize::colorRamp2(
      breaks = seq(0.85, 1.15, length.out = 7),
      colors = rev(brewer.pal(7, "RdBu"))
    ),
    heatmap_legend_param = list(
      at = c(0.85, 1.15),
      border = FALSE,
      grid_width = unit(1.5, "mm"),
      labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ),
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    row_split = cell_metadata$cellont_abbr,
    row_gap = unit(0, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    row_dend_gp = gpar(lwd = 0.5),
    
    cluster_columns = FALSE,
    column_split = col_split,
    cluster_column_slices = FALSE,
    column_gap = unit(0, "mm"),
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    border = TRUE,
    border_gp = gpar(lwd = 0.5),
    
    left_annotation = rowAnnotation(
      group = cell_metadata$group,
      patient = cell_metadata$sample,
      col = list(group = GROUP_COLORS, patient = PATIENT_COLORS),
      annotation_name_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      annotation_name_side = "bottom",
      annotation_legend_param = list(
        group = list(
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        ),
        patient = list(
          grid_width = unit(2, "mm"),
          labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
          title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
        )
      )
    ),
    
    use_raster = FALSE,
  ) %>%
    draw(
      row_title = "cell type",
      row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      column_title = "chromosome",
      column_title_side = "bottom",
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ) 
}

(p <- plot_resexp_marrow())
ggsave_publication("S1d_resexp_marrow", plot = p,
                   type = "png", width = 18, height = 7.5)



# Tables ------------------------------------------------------------------

## S2 ----

dir_ls("data_raw/", recurse = TRUE, regex = "metrics_summary.csv$") %>% 
  map_dfr(read_csv, .id = "file") %>%
  rename_with(~str_glue("{.x} (%)"), !file & where(is.character)) %>% 
  mutate(across(!file & where(is.character), parse_number)) %>%
  extract(file, into = "sample", regex = "([\\w\\d_]*)_trans") %>% 
  left_join(
    read_csv("metadata/sample_groups.csv", comment = "#") %>% 
      select(bsf_id, sample, group) %>%
      mutate(sample = rename_patients(sample), group = rename_groups(group)),
    by = c(sample = "bsf_id")
  ) %>% 
  select(
    Sample = sample.y,
    `Estimated Number of Cells`,
    `Mean Reads per Cell`,
    `Median Genes per Cell`,
    `Number of Reads`,
    `Valid Barcodes (%)`,
    `Total Genes Detected`,
    `Median UMI Counts per Cell`,
    `Sequencing Saturation (%)`
  ) %>%
  filter(!is.na(Sample)) %>% 
  mutate(
    Sample = factor(Sample, levels = names(PATIENT_COLORS)),
    Group = GROUP_NAMES_LONG[str_sub(Sample, 1, 1)],
    .after = Sample
  ) %>% 
  arrange(Sample) %>% 
  save_table("S2_sequencing_statistics", "Sequencing Statistics")



## S3 ----

infercnv_data$regions_data %>% 
  filter(near(prob, 0.5)) %>%
  extract(
    cell_group_name,
    into = c("sample", "group"),
    regex = "malignant_(.*?)_([IV]+)"
  ) %>%
  filter(!is.na(sample)) %>% 
  transmute(
    sample = sample,
    type = "scRNA-seq",
    chr = str_sub(chr, 4) %>%
      factor(levels = snparray_data$chromosome_size$chr),
    start = as.integer(start),
    end = as.integer(end),
    delta_copy_number = as.numeric(state) - 3,
  ) %>%
  bind_rows(
    snparray_data$cnv_regions %>%
      select(!type:ploidy) %>%
      mutate(type = "SNP array")
  ) %>%
  mutate(
    sample = rename_patients(sample),
    group = GROUP_NAMES_LONG[str_sub(sample, 1, 1)],
    .after = sample
  ) %>%
  rename(
    Patient = sample,
    Group = group,
    Type = type,
    Chromosome = chr,
    Start = start,
    End = end,
    "Copy Number Difference" = delta_copy_number
  ) %>%
  arrange(Patient, Type, Chromosome, Start) %>%
  save_table("S3_cnv", "CNV")
