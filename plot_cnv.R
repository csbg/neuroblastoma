# Plot copy number variations.
#
# @DEPI infercnv_output

library(GenomeInfoDb)
library(infercnv)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyverse)
library(fs)
source("common_functions.R")
source("styling.R")

ht_opt(message = FALSE, show_parent_dend_line = FALSE)



# Load data ---------------------------------------------------------------

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

infercnv_data <- read_infercnv_data("data_generated/infercnv_output")
nb_data <- readRDS("data_generated/metadata.rds")

snparray_sample_names <- 
  read_csv("metadata/sample_groups.csv", comment = "#") %>% 
  select(snp_array_id, sample) %>% 
  filter(!is.na(snp_array_id)) %>% 
  deframe()
snparray_data <- read_snparray_data("data_raw/snp_array", snparray_sample_names)



# Plot residual expression ------------------------------------------------

#' Plot residual expression as calculated by infercnv.
#'
#' @param data infercnv data loaded by `read_infercnv_data()`.
#' @param metadata Cell metadata from `assemble_metadata.R`.
#' @param cells Plot either the "tumor" or "ref"erence cells.
#' @param cells_per_row_split If not `NULL`, only plot a limited number of cells
#'   per row plit, which are randomly selected.
#' @param seed Seed for random selection of cells.
#' @param row_title Heatmap row title; if `NULL`, use a default title.
#' @param filename Name of output file.
#'
#' @return A character vector giving the order of samples in the heatmap.
plot_residual_expression <- function(data, metadata,
                                     cells = c("tumor", "ref", "input"),
                                     cells_per_row_split = 50L, seed = 42,
                                     row_title = NULL, filename = NULL) {
  cells <- match.arg(cells)
  
  # colors for group and sample annotations
  group_colors <- c(  # colours.cafe 502
    I = "#f8f4e8",
    II = "#154d42",
    III = "#dcbc65",
    IV = "#b7336c"
  )
  sample_colors <-
    rainbow(nlevels(metadata$sample)) %>%
    set_names(levels(metadata$sample))
  
  # row titles
  if (is.null(row_title)) {
    row_title <- c(
      tumor = "sample",
      ref = "cell type",
      input = ""
    )[cells]  
  }
  
  # prepare cell metadata
  if (cells == "tumor") {
    # order of tumor samples
    tumor_samples <-
      metadata %>% 
      filter(group != "I") %>% 
      pull(sample) %>% 
      fct_drop() %>% 
      levels()
    
    cell_metadata <- 
      data$final_obj@tumor_subclusters$subclusters %>%
      enframe() %>% 
      filter(str_starts(name, "malignant")) %>% 
      extract(
        name,
        into = "row_split",
        regex = "malignant_(.*)_I"
      ) %>% 
      mutate(
        row_split = as_factor(row_split) %>% fct_relevel(tumor_samples)
      ) %>% 
      arrange(row_split)
  } else if (cells == "ref") {
    cell_metadata <- 
      data$final_obj@tumor_subclusters$subclusters %>%
      enframe(name = "row_split") %>% 
      filter(!str_starts(row_split, "malignant"))
  } else {
    # if cells == "input"
    cell_indices <- 
      data$final_obj@tumor_subclusters$subclusters %>% 
      map_dfr(~enframe(.[[1]], "cell", "cell_index"))
    cell_metadata <- 
      metadata %>% 
      select(cell, sample, group, row_split = cluster_50) %>% 
      left_join(cell_indices, by = "cell")
  }
  
  # post-processing for tumor and ref cells
  if (cells != "input") {
    cell_metadata <- 
      cell_metadata %>% 
      unnest_longer(value) %>%
      select(!value_id) %>%
      unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>% 
      left_join(metadata %>% select(cell, sample, group), by = "cell")
  }
  
  # cell_metadata is now a dataframe with five columns:
  # cell, sample, group (as usual),
  # cell_index (index in residual expression matrix), and
  # row_split (factor for splitting the heatmap rows).
  # return(cell_metadata)
  
  # optional: subset n cells per row_split
  if (!is.null(cells_per_row_split)) {
    set.seed(seed)
    cell_metadata <- 
      cell_metadata %>% 
      group_by(row_split) %>% 
      slice_sample(n = cells_per_row_split) %>% 
      ungroup()
  }
  
  # matrix with residual expression, rows ordered
  cnv_mat <-
    data$final_obj@expr.data %>% 
    magrittr::extract(, cell_metadata$cell_index) %>% 
    t()
  
  # annotation
  left_annotation <- rowAnnotation(
    sample = cell_metadata$sample,
    group = cell_metadata$group,
    col = list(sample = sample_colors, group = group_colors)
  )
  
  # draw heatmap
  p <- 
    Heatmap(
      cnv_mat,
      name = "residual\nexpression",
      col = circlize::colorRamp2(
        breaks = seq(0.85, 1.15, length.out = 7),
        colors = rev(brewer.pal(7, "RdBu"))
      ),
      show_row_names = FALSE,
      show_column_names = FALSE,
      
      row_split = cell_metadata$row_split,
      row_gap = unit(0, "mm"),
      row_title_rot = 0,
      
      cluster_columns = FALSE,
      column_split =
        data$final_obj@gene_order$chr %>%
        str_sub(4) %>%
        as_factor(),
      cluster_column_slices = FALSE,
      column_gap = unit(0, "mm"),
      
      border = TRUE,
      
      left_annotation = left_annotation,
      
      use_raster = FALSE,
    ) %>%
    draw(
      row_title = row_title,
      column_title = "chromosome"
    ) 
  
  ggsave_default(filename, plot = p)
  invisible(names(row_order(p)))
}

# tumor cells
sample_order <- plot_residual_expression(
  infercnv_data, nb_data, filename = "cnv/residual_expr_tumor"
)
sample_order

# reference cells
plot_residual_expression(
  infercnv_data, nb_data, cells = "ref", filename = "cnv/residual_expr_ref"
)

# small clusters
plot_residual_expression(
  infercnv_data,
  nb_data %>% filter(as.numeric(cluster_50) >= 15),
  cells = "input",
  row_title = "cluster",
  cells_per_row_split = NULL,
  filename = "cnv/residual_expr_clusters"
)



# Plot CNV regions --------------------------------------------------------

#' Plot CNV regions detected by the HMM algorithm.
#'
#' @param data infercnv data loaded by `read_infercnv_data()`.
#' @param sample_order Order of samples on the y-axis.
#' @param probability Show results for this cutoff probability, i.e., CNV
#'   regions whose posterior probability of being normal exceeds this value.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cnv_regions <- function(data, sample_order,
                             probability = 0.5, filename = NULL) {
  # get chromosome sizes from GenomeInfoDb
  chromosome_size <- 
    Seqinfo(genome = "GRCh38.p13") %>%
    seqlengths() %>% 
    enframe("chr", "end") %>% 
    filter(str_detect(chr, "^(\\d\\d?|X|Y)$")) %>%
    mutate(chr = as_factor(chr) %>% fct_inseq())
  
  # descriptions and colors for HMM states
  hmm_states <- tribble(
    ~label, ~state, ~color,
    "complete loss", "1", "#4393c3",
    "loss of one copy", "2", "#92c5de",
    "neutral", "3", "#f7f7f7",
    "addition of one copy", "4", "#fddbc7",
    "addition of two copies", "5", "#ef8a62",
    "addition of more than two copies", "6", "#b2182b"
  )
  hmm_labels <-
    hmm_states %>%
    select(!color) %>%
    deframe()
  hmm_colors <-
    hmm_states %>%
    select(!state) %>%
    deframe()
  
  # preprocess data
  plot_data <-
    data$regions_data %>% 
    filter(near(prob, probability)) %>% 
    extract(
      cell_group_name,
      into = c("sample", "group"),
      regex = "malignant_(.*?)_([IV]+)"
    ) %>% 
    arrange(group, sample) %>% 
    mutate(
      sample = as_factor(sample) %>% fct_relevel(sample_order) %>% fct_rev(),
      state = as_factor(state) %>% fct_inseq() %>% fct_recode(!!!hmm_labels),
      chr = chr %>% str_sub(4) %>% as_factor() %>% fct_inseq(),
      ymin = as.numeric(sample) - 0.5,
      ymax = ymin + 1
    )
  
  # make plot
  p <- 
    ggplot(chromosome_size, aes(xmax = end)) +
    geom_rect(aes(xmin = 0, ymin = -0.1, ymax = 0)) +
    geom_rect(
      data = plot_data,
      aes(xmin = start, ymin = ymin, ymax = ymax, fill = state)
    ) +
    scale_x_continuous(
      name = "chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "sample",
      breaks = seq_along(sample_order),
      labels = rev(sample_order),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = hmm_colors,
      guide = guide_legend(nrow = 1),
      drop = FALSE
    ) +
    coord_cartesian(ylim = c(0.5, length(sample_order) + 0.5)) +
    facet_grid(
      cols = vars(chr),
      space = "free_x",
      scales = "free_x",
      switch = "x"
    ) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(0, "mm"),
      panel.background = element_rect(color = NA, fill = cn_colors["neutral"]),
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_blank()
    )
  
  ggsave_default(filename, width = 297, height = 100)
  p
}

plot_cnv_regions(infercnv_data, sample_order)

walk(
  unique(infercnv_data$regions_data$prob),
  ~plot_cnv_regions(
    infercnv_data,
    sample_order,
    probability = .,
    filename = str_glue("cnv/regions_{.}")
  )
)



# Compare to SNP array data -----------------------------------------------

#' Compare CNVs detected in scRNA-seq data to CNVs from SNP array data.
#'
#' @param sc_data scRNA-seq data as loaded by `read_infercnv_data()`.
#' @param snp_data  SNP array data as loaded by `read_snparray_data()`.
#' @param selected_samples Samples to include. If `NULL`, include all samples
#'   with complete data.
#' @param probability Filtering probability for CNVs predicted from scRNA-seq.
#' @param logrr_prop Proportion of log R ratio data points to plot.
#' @param filename Name of output file.
#'
#' @return A ggplot object.
plot_cnv_data_comparison <- function(sc_data,
                                     snp_data,
                                     selected_samples = NULL,
                                     probability = 0.5,
                                     logrr_prop = 0.01,
                                     filename = NULL) {
  # labels and colors for copy numbers
  cn_metadata <-  # colorbrewer 11-class BrBG 
    tribble(
      ~delta_copy_number, ~label,         ~color,
      -2, "complete loss",                "#01665EFF",
      -1, "loss of one copy",             "#80CDC1FF",
       0, "neutral",                      "#F7F7F7FF",
       1, "gain of one copy",             "#DFC27DFF",
       2, "gain of two copies",           "#BF812DFF",
       3, "gain of more than two copies", "#8C510AFF",
       38, "amplification",               "#543005FF"
    ) %>%
    mutate(label = as_factor(label))
  cn_color <- circlize::colorRamp2(
    cn_metadata$delta_copy_number,
    cn_metadata$color
  )
  
  # prepare data for plotting
  plot_data_regions <-
    sc_data$regions_data %>% 
    filter(near(prob, probability)) %>% 
    extract(
      cell_group_name,
      into = c("sample", "group"),
      regex = "malignant_(.*?)_([IV]+)"
    ) %>% 
    transmute(
      sample = sample,
      type = "sc",
      chr = str_sub(chr, 4) %>% factor(levels = snp_data$chromosome_size$chr),
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
    mutate(fill = cn_color(delta_copy_number))
  # return(plot_data_regions)
  
  plot_data_logrr <-
    snp_data$logrr_data %>% 
    group_by(sample) %>% 
    slice_sample(prop = logrr_prop)
  # return(plot_data_logrr)
  
  # filter for selected or common samples
  if (is.null(selected_samples)) {
    selected_samples <- 
      plot_data_regions %>%
      group_by(type) %>% 
      summarise(samples = list(unique(sample))) %>% 
      deframe() %>% 
      {intersect(.$sc, .$snp)}
  }
  plot_data_regions <-
    plot_data_regions %>%
    filter(sample %in% selected_samples) %>% 
    mutate(sample = factor(sample, level = selected_samples))
  plot_data_logrr <-
    plot_data_logrr %>%
    filter(sample %in% selected_samples) %>% 
    mutate(sample = factor(sample, level = selected_samples))
  
  # make plot
  p <-
    ggplot(NULL, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
    geom_rect(
      data = snp_data$chromosome_size,
      aes(xmin = 0, ymin = -0.1, ymax = 0)
    ) +
    geom_rect(
      data = plot_data_regions,
      aes(fill = fill)
    ) +
    geom_line(
      data = plot_data_logrr,
      aes(x = start, y = smoothed + 1.5),
      inherit.aes = FALSE,
      size = .25,
      alpha = .5
    ) +
    geom_hline(yintercept = 1, size = 0.1) +
    scale_x_continuous(
      name = "chromosome",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = NULL,
      breaks = c(0.5, 1.5),
      labels = c("scRNA-seq", "SNP array"),
      expand = c(0, 0)
    ) +
    scale_fill_identity(
      name = NULL,
      guide = guide_legend(nrow = 1),
      breaks = cn_metadata$color,
      labels = cn_metadata$label
    ) +
    coord_cartesian(ylim = c(0, 2)) +
    facet_grid(
      cols = vars(chr),
      rows = vars(sample),
      space = "free_x",
      scales = "free_x",
      switch = "x"
    ) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "bottom",
      panel.spacing.x = unit(0, "mm"),
      panel.background = element_rect(color = NA, fill = cn_color(0)),
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_blank()
    )

  ggsave_default(filename, width = 420, height = 400)
  p
}

plot_cnv_data_comparison(infercnv_data, snparray_data,
                         # selected_samples = "2019_2495",
                         # logrr_prop = 0.001,
                         filename = "cnv/comparison")



# Publication figures -----------------------------------------------------

## Figure 1e ----

plot_resexp_tumor <- function(data,
                              metadata,
                              cells_per_sample = 50L) {
  # cell metadata
  cell_metadata <- 
    data$final_obj@tumor_subclusters$subclusters %>%
    enframe() %>%
    filter(str_starts(name, "malignant")) %>%
    unnest_longer(value) %>%
    select(value) %>%
    unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>%
    left_join(nb_data %>% select(cell, sample, group), by = "cell") %>%
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
  
  # matrix with residual expression, rows ordered
  cnv_mat <-
    data$final_obj@expr.data %>% 
    magrittr::extract(, cell_metadata$cell_index) %>% 
    t()
  
  # draw heatmap
  Heatmap(
    cnv_mat,
    name = "residual\nexpression",
    col = circlize::colorRamp2(
      breaks = seq(0.85, 1.15, length.out = 7),
      colors = rev(brewer.pal(7, "RdBu"))
    ),
    heatmap_legend_param = list(
      border = FALSE,
      grid_width = unit(2, "mm"),
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
    column_split =
      data$final_obj@gene_order$chr %>%
      str_sub(4) %>%
      as_factor(),
    cluster_column_slices = FALSE,
    column_gap = unit(0, "mm"),
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    border = TRUE,
    border_gp = gpar(lwd = 0.5),
    
    left_annotation = rowAnnotation(
      group = cell_metadata$group,
      col = list(group = GROUP_COLORS),
      annotation_name_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      annotation_name_side = "top",
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
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ) 
}

ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(1, "mm")
)
p <- plot_resexp_tumor(infercnv_data, nb_data, 50L)
ggsave_publication("1e_resexp_tumor", plot = p,
                   type = "png", width = 18, height = 8)


## Figure 1f ----

plot_cnv_data_comparison <- function(sc_data,
                                     snp_data,
                                     selected_samples = NULL,
                                     probability = 0.5,
                                     logrr_prop = 0.01,
                                     filename = NULL) {
  # labels and colors for copy numbers
  cn_metadata <-  # colorbrewer 11-class BrBG 
    tribble(
      ~delta_copy_number, ~label,        ~color,
      -2, "complete loss",               "#01665EFF",
      -1, "loss of one copy",            "#80CDC1FF",
      0, "neutral",                      "#F7F7F7FF",
      1, "gain of one copy",             "#DFC27DFF",
      2, "gain of two copies",           "#BF812DFF",
      3, "gain of more than two copies", "#8C510AFF",
      38, "amplification",               "#543005FF"
    ) %>%
    mutate(label = as_factor(label))
  cn_color <- circlize::colorRamp2(
    cn_metadata$delta_copy_number,
    cn_metadata$color
  )
  
  # prepare data for plotting
  plot_data_regions <-
    sc_data$regions_data %>% 
    filter(near(prob, probability)) %>% 
    extract(
      cell_group_name,
      into = c("sample", "group"),
      regex = "malignant_(.*?)_([IV]+)"
    ) %>% 
    transmute(
      sample = sample,
      type = "sc",
      chr = str_sub(chr, 4) %>% factor(levels = snp_data$chromosome_size$chr),
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
    mutate(fill = cn_color(delta_copy_number))
  # return(plot_data_regions)
  
  plot_data_logrr <-
    snp_data$logrr_data %>% 
    group_by(sample) %>% 
    slice_sample(prop = logrr_prop)
  # return(plot_data_logrr)
  
  # filter for selected or common samples
  if (is.null(selected_samples)) {
    selected_samples <- 
      plot_data_regions %>%
      group_by(type) %>% 
      summarise(samples = list(unique(sample))) %>% 
      deframe() %>% 
      {intersect(.$sc, .$snp)}
  }
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
      data = snp_data$chromosome_size,
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
      guide = guide_legend(nrow = 1),
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
      legend.margin = margin(0, 1, -2, 1, "mm"),
      legend.position = "top",
      panel.spacing.x = unit(-.5, "pt"),
      panel.spacing.y = unit(1, "mm"),
      panel.background = element_rect(color = NA, fill = cn_color(0)),
      strip.placement = "outside",
      strip.text.x = element_text(margin = margin()),
      strip.text.y.left = element_text(angle = 0, margin = margin())
    )
  
  p
}

plot_cnv_data_comparison(
  infercnv_data, snparray_data,
  selected_samples = c("2016_4503", "2005_1702", "2006_2684")
)
ggsave_publication("1f_cnv_comparison", width = 18, height = 6)



## Figure S1d ----

plot_resexp_marrow <- function(data,
                               metadata,
                               cells_per_type = 50L) {
  # cell metadata
  cell_metadata <- 
    data$final_obj@tumor_subclusters$subclusters %>%
    enframe(name = "row_split") %>% 
    filter(!str_starts(row_split, "malignant")) %>% 
    unnest_longer(value) %>%
    select(value) %>%
    unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>%
    left_join(
      nb_data %>%
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
  
  # matrix with residual expression, rows ordered
  cnv_mat <-
    data$final_obj@expr.data %>% 
    magrittr::extract(, cell_metadata$cell_index) %>% 
    t()
  
  # draw heatmap
  Heatmap(
    cnv_mat,
    name = "residual\nexpression",
    col = circlize::colorRamp2(
      breaks = seq(0.85, 1.15, length.out = 7),
      colors = rev(brewer.pal(7, "RdBu"))
    ),
    heatmap_legend_param = list(
      border = FALSE,
      grid_width = unit(2, "mm"),
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
    column_split =
      data$final_obj@gene_order$chr %>%
      str_sub(4) %>%
      as_factor(),
    cluster_column_slices = FALSE,
    column_gap = unit(0, "mm"),
    column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
    
    border = TRUE,
    border_gp = gpar(lwd = 0.5),
    
    left_annotation = rowAnnotation(
      group = cell_metadata$group,
      patient = cell_metadata$sample,
      col = list(group = GROUP_COLORS, patient = PATIENT_COLORS),
      annotation_name_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
      annotation_name_side = "top",
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
      column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT)
    ) 
}

p <- plot_resexp_marrow(infercnv_data, nb_data)
ggsave_publication("S1d_resexp_marrow", plot = p,
                   type = "png", width = 18, height = 8)
