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

infercnv_data <- read_infercnv_data("data_generated/infercnv_output")
nb_data <- readRDS("data_generated/metadata.rds")



# Plot residual expression ------------------------------------------------

#' Plot residual expression as calculated by infercnv.
#'
#' @param data infercnv data loaded by `read_infercnv_data()`.
#' @param metadata Cell metadata from `assemble_metadata.R`.
#' @param cells Plot either the "tumor" or "ref"erence cells.
#' @param cells_per_sample If not `NULL`, only plot a limited number of cells
#'   per sample, which are randomly selected.
#' @param seed Seed for random selection of cells.
#' @param filename Name of output file.
#'
#' @return A character vector giving the order of samples in the heatmap.
plot_residual_expression <- function(data, metadata, cells = c("tumor", "ref"),
                                     cells_per_sample = 50L, seed = 42,
                                     filename = NULL) {
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
  
  # data frame with tumor or reference cell indices
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
        into = "sample",
        regex = "malignant_(.*)_I"
      ) %>% 
      mutate(sample = as_factor(sample) %>% fct_relevel(tumor_samples)) %>% 
      arrange(sample)
  } else {
    cell_metadata <- 
      data$final_obj@tumor_subclusters$subclusters %>%
      enframe(name = "sample") %>% 
      filter(!str_starts(sample, "malignant"))
  }
  
  cell_metadata <- 
    cell_metadata %>% 
    unnest_longer(value) %>%
    select(!value_id) %>%
    unnest_longer(value, values_to = "cell_index", indices_to = "cell") %>% 
    left_join(metadata %>% select(cell, patient = sample, group), by = "cell")
  
  # optional: subset n cells per sample
  if (!is.null(cells_per_sample)) {
    set.seed(seed)
    cell_metadata <- 
      cell_metadata %>% 
      group_by(sample) %>% 
      slice_sample(n = cells_per_sample) %>% 
      ungroup()
  }
  
  # matrix with residual expression, rows ordered
  cnv_mat <-
    data$final_obj@expr.data %>% 
    magrittr::extract(, cell_metadata$cell_index) %>% 
    t()
  
  # annotation (group for tumor cells, sample for reference cells)
  if (cells == "tumor") {
    left_annotation <- rowAnnotation(
      group = cell_metadata$group,
      col = list(group = group_colors)
    )
  } else {
    patient_ids <- unique(cell_metadata$patient)
    left_annotation <- rowAnnotation(
      sample = cell_metadata$patient,
      group = cell_metadata$group,
      col = list(sample = sample_colors, group = group_colors)
    )
  }
  
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
      
      row_split = cell_metadata$sample,
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
      
      left_annotation = left_annotation
    ) %>%
    draw(
      row_title = if (cells == "tumor") "sample",
      column_title = "chromosome"
    ) 
  
  ggsave_default(filename, plot = p)
  invisible(names(row_order(p)))
}

# tumor cells
sample_order <- plot_residual_expression(
  infercnv_data, nb_data, filename = "cnv/residual_expr_tumor"
)

# reference cells
plot_residual_expression(
  infercnv_data, nb_data, cells = "ref", filename = "cnv/residual_expr_ref"
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
      panel.background = element_rect(color = NA, fill = "#f7f7f7"),
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
