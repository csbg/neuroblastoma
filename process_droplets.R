library(SingleCellExperiment)
library(scuttle)
library(DropletUtils)
library(tidyverse)
library(fs)
library(patchwork)
library(ggbeeswarm)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# folder that contains cellranger results
data_dir <- "data_raw/COUNT"

# path where the merged dataset is saved
out_file <- "data_generated/rna_merged_NEW.rds"

# data frame with three columns 'order', 'sample_name', and 'sample_file'
samples <-
  tibble(
    # sample_file = dir_ls(
    #   data_dir,
    #   recurse = TRUE,
    #   regexp = "raw_feature_bc_matrix.h5"
    # ),
    sample_file = dir_ls(
      data_dir,
      type = "directory",
      recurse = TRUE,
      regexp = "raw_feature_bc_matrix"
    ),
    bsf_id = str_match(sample_file, "([\\w\\d_]*)_trans")[, 2]
  ) %>%
  left_join(
    read_csv("metadata/sample_groups.csv", comment = "#"),
    by = "bsf_id"
  ) %>% 
  arrange(bsf_order) %>% 
  select(order = bsf_order, sample_file, sample_name = sample) %>% 
  drop_na()



# Functions ---------------------------------------------------------------

# for testing
# sample_file <- "data_raw/COUNT/MF220_NB_BM_Patient1_transcriptome/raw_feature_bc_matrix"
# sample_name <- "2018_1404"
# order <- 15
# FDR <- 0.001

process_droplets <- function(sample_file, sample_name, order, FDR = 0.001) {
  
  # load data
  info("Loading data in {sample_file}")
  sce <- read10xCounts(sample_file, sample.names = sample_name)
  sce$Barcode <- str_replace(sce$Barcode, "\\d+$", as.character(order))
  metadata(sce) <- list(
    sample_name = sample_name,
    sample_file = sample_file,
    bsf_order = order
  )
  rownames(sce) <- rowData(sce)$Symbol
  
  # barcode ranks
  info("  Calculating barcode ranks")
  barcode_ranks <- barcodeRanks(counts(sce))
  
  # empty droplets
  info("  Determining empty droplets")
  set.seed(42)
  empty_droplets <- emptyDrops(counts(sce), test.ambient = TRUE)
  
  info(
    "  Keeping {sum(empty_droplets$FDR <= FDR, na.rm = TRUE)} ",
    "of {nrow(empty_droplets)} droplets"
  )
  sce <- sce[, which(empty_droplets$FDR <= FDR)]
  
  # QC filtering
  info("  Performing QC filtering")
  is_mito <- rowData(sce)$Symbol %>% str_starts("MT-")
  qc_data <- perCellQCMetrics(sce, subsets = list(mito = is_mito))
  discard <- quickPerCellQC(qc_data, sub.fields = "subsets_mito_percent")$discard
  
  qc_data <- 
    qc_data %>% 
    as_tibble() %>% 
    select(
      library_size = sum,
      n_features = detected,
      percent_mito = subsets_mito_percent
    ) 
  
  info("  Keeping  {sum(!discard)} of {nrow(qc_data)} cells")
  colData(sce) <- cbind(colData(sce), qc_data)
  sce <- sce[, !discard]
  sce <- sce[rowSums(counts(sce)) > 0, ]
  
  # return values
  list(
    sce = sce,
    barcode_ranks = barcode_ranks,
    empty_droplets = empty_droplets,
    qc_data = qc_data %>% mutate(retain = !discard)
  )
}


plot_droplet_metrics <- function(data, FDR = 0.001, save_plot = FALSE) {
  vis_data <- 
    left_join(
      as_tibble(data$barcode_ranks, rownames = "drop_id"),
      as_tibble(data$empty_droplets, rownames = "drop_id"),
      by = "drop_id"
    ) %>% 
    mutate(retained = (FDR <= {{FDR}}) %>% replace_na(FALSE)) %>% 
    select(drop_id, rank, total, retained, p_adj = PValue, log_prob = LogProb) 
  
  default_color_scheme <- scale_color_manual(values = c("#dd1c77", "black"))
  default_theme <- 
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  # ranks
  p1 <-
    vis_data %>%
    filter(total > 0) %>%
    ggplot(aes(rank, total)) +
    geom_point(
      aes(color = retained),
      alpha = 0.25,
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = flatten_dbl(metadata(data$barcode_ranks)),
      linetype = c("solid", "dashed")
    ) +
    scale_x_log10("barcode rank") +
    scale_y_log10("total UMI count") +
    ggtitle("Droplet processing") +
    default_color_scheme +
    default_theme

  # p value histogram
  p2 <-
    vis_data %>%
    filter(total <= 100, total > 0) %>%
    ggplot(aes(p_adj)) +
    geom_histogram(bins = 20, fill = "gray80") +
    xlab("adjusted p value") +
    default_theme
  
  # counts vs log probability
  p3 <-
    vis_data %>%
    filter(total > 0) %>%
    ggplot(aes(total, -log_prob)) +
    geom_point(
      aes(color = retained),
      alpha = 0.25,
      size = 0.5,
      show.legend = FALSE
    ) +
    default_color_scheme +
    xlab("total UMI count") +
    ylab("-log probability") +
    default_theme
  
  # qc violin plots
  p4 <- 
    data$qc_data %>% 
    pivot_longer(!retain) %>% 
    ggplot(aes("1", value)) +
    geom_violin(color = "gray60", alpha = 0.2, width = 0.8) +
    geom_quasirandom(
      aes(color = retain),
      width = 0.4,
      bandwidth = 1,
      alpha = 0.25,
      size = 0.5,
      show.legend = FALSE
    ) +
    default_color_scheme +
    ylab("") +
    facet_wrap(vars(name), scales = "free_y") +
    ggtitle("Quality control") +
    default_theme +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  # %mt vs library size
  p5 <- 
    data$qc_data %>% 
    ggplot(aes(library_size, percent_mito)) +
    geom_point(
      aes(color = retain),
      alpha = 0.25,
      size = 0.5,
      show.legend = FALSE
    ) +
    default_color_scheme +
    default_theme
  
  # %mt vs feature count
  p6 <- 
    data$qc_data %>% 
    ggplot(aes(n_features, percent_mito)) +
    geom_point(
      aes(color = retain),
      alpha = 0.25,
      size = 0.5,
      show.legend = FALSE
    ) +
    default_color_scheme +
    default_theme
  
  layout <- "
    AADD
    BCEF"
  
  p <-
    wrap_plots(p1, p2, p3, p4, p5, p6) +
    plot_layout(design = layout) +
    plot_annotation(
      title = str_glue(
        "{data$sce@metadata$sample_name} ",
        "(file {data$sce@metadata$sample_file})"
      ),
      subtitle = str_glue(
        "{nrow(vis_data)} droplets → ",
        "{sum(vis_data$retained)} cells → ",
        "{ncol(data$sce)} cells pass QC ",
        "({nrow(data$sce)} features)"
      ),
      caption = str_glue(
        "red, discarded cells; ",
        "solid line, knee; ",
        "dashed line, inflection"
      )
    )
  if (save_plot)
    ggsave_default(
      str_glue(
        "qc/",
         "{data$sce@metadata$bsf_order}_",
         "{data$sce@metadata$sample_name}"
      ),
      height = 180
    )
  p
}



# Analysis ----------------------------------------------------------------

# for testing
# nb <- pmap(samples[15, ], process_droplets)
# nb[[1]]$sce

nb <- pmap(samples, process_droplets)

nb %>% walk(plot_droplet_metrics, save_plot = TRUE)



# Save results ------------------------------------------------------------

nb %>% 
  map(~.x$sce) %>% 
  saveRDS(out_file)


# xx <- readRDS("data_generated/rna_merged_NEW.rds")
# 
# xx %>% map_dbl(ncol) %>% sum()
# 
# xx
