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
      copy_number = X13
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

plot_cnv_data_comparison <- function(sc_data,
                                     snp_data,
                                     probability = 0.5,
                                     logrr_prop = 0.01,
                                     filename = NULL) {
  # labels and colors for copy numbers
  cn_metadata <-  # colorbrewer 11-class BrBG 
    tribble(
      ~copy_number, ~label, ~color,
      0L, "complete loss", "#01665e",
      1L, "loss of one copy", "#80cdc1",
      NA, "mosaic loss", "#c7eae5",
      2L, "neutral", "#f7f7f7",
      NA, "mosaic gain", "#f6e8c3",
      3L, "gain of one copy", "#dfc27d",
      4L, "gain of two copies", "#bf812d",
      5L, "gain of more than two copies", "#8c510a",
      NA, "amplification", "#543005"
    ) %>%
    mutate(label = as_factor(label))
  cn_colors <- 
    cn_metadata %>% 
    select(label, color) %>% 
    deframe()
  
  # prepare data for plotting
  plot_data_sc <-
    sc_data$regions_data %>% 
    filter(near(prob, probability)) %>% 
    extract(
      cell_group_name,
      into = c("sample", "group"),
      regex = "malignant_(.*?)_([IV]+)"
    ) %>% 
    transmute(
      sample = sample,
      chr = str_sub(chr, 4) %>% factor(levels = snp_data$chromosome_size$chr),
      start = as.integer(start),
      end = as.integer(end),
      copy_number = as.integer(state) - 1L
    ) %>% 
    left_join(cn_metadata, by = "copy_number")
  # return(plot_data_sc)
  
  plot_data_snp_regions <-
    snp_data$cnv_regions %>% 
    mutate(
      label = case_when(
        type == "amplification" ~ "amplification",
        type == "gain" & near(copy_number, 5) ~ "gain of more than two copies",
        type == "gain" & near(copy_number, 4) ~ "gain of two copies",
        type == "gain" & near(copy_number, 3) ~ "gain of one copy",
        type == "mosaic_gain" ~ "mosaic gain",
        type == "mosaic_loss" ~ "mosaic loss",
        type == "loss" & near(copy_number, 1) ~ "loss of one copy",
        type == "loss" & near(copy_number, 0) ~ "complete loss",
        TRUE ~ NA_character_
      )
    ) %>% 
    left_join(cn_metadata, by = "label")
  # return(plot_data_snp_regions)
  
  plot_data_snp_logrr <-
    snp_data$logrr_data %>% 
    group_by(sample) %>% 
    slice_sample(prop = logrr_prop)
  # return(plot_data_snp_logrr)
  
  # filter for common samples
  common_samples <- intersect(
    plot_data_sc$sample,
    plot_data_snp_regions$sample
  )
  plot_data_sc <-
    plot_data_sc %>%
    filter(sample %in% common_samples)
  plot_data_snp_regions <-
    plot_data_snp_regions %>%
    filter(sample %in% common_samples)
  plot_data_snp_logrr <-
    plot_data_snp_logrr %>%
    filter(sample %in% common_samples)
  
  # make plot
  p <-
    ggplot(NULL, aes(xmax = end)) +
    geom_rect(
      data = snp_data$chromosome_size,
      aes(xmin = 0, ymin = -0.1, ymax = 0)
    ) +
    geom_rect(
      data = plot_data_sc,
      aes(xmin = start, fill = label),
      ymin = 0,
      ymax = 1
    ) +
    geom_rect(
      data = plot_data_snp_regions,
      aes(xmin = start, fill = label),
      ymin = 1,
      ymax = 2
    ) +
    geom_line(
      data = plot_data_snp_logrr,
      aes(x = start, y = smoothed + 1.5),
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
    scale_fill_manual(
      name = NULL,
      values = cn_colors,
      guide = guide_legend(nrow = 1),
      drop = FALSE
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
      panel.background = element_rect(color = NA, fill = cn_colors["neutral"]),
      panel.border = element_rect(color = "black", fill = NA),
      strip.background = element_blank()
    )

  ggsave_default(filename, width = 420, height = 400)
  p
}

plot_cnv_data_comparison(infercnv_data, snparray_data,
                         filename = "cnv/comparison")
