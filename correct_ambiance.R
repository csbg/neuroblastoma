# Creates a cell_data_set with an additional assay `soupx_counts`, containing
# decontaminated counts as determined by SoupX.
#
# @DEPI raw data
# @DEPI rna_integrated_monocle.rds
# @DEPO rna_decontaminated.rds

library(SingleCellExperiment)
library(DropletUtils)
library(monocle3)
library(SoupX)
library(tidyverse)
library(fs)
library(ggpmisc)
library(patchwork)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# folder that contains cellranger results
data_dir <- "data_raw/COUNT"

# integrated dataset
in_file <- "data_generated/rna_integrated_monocle.rds"

# dataset with decontaminated counts
out_file <- "data_generated/rna_decontaminated.rds"

# data frame with three columns 'order', 'sample_name', and 'sample_file'
samples <-
  tibble(
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



# Load data ---------------------------------------------------------------

nb <- readRDS(in_file)
nb_metadata <- readRDS("data_generated/metadata.rds")



# Analysis ----------------------------------------------------------------

remove_soup <- function(bsf_order) {
  info("Removing soup for ID {bsf_order}")
  
  # subset filtered cells
  sce_cells <- nb[, str_ends(colnames(nb), str_glue("-{bsf_order}"))]
  sce_cells
  
  # load raw droplet data
  raw_file <-
    samples %>%
    filter(order == {{bsf_order}}) %>% 
    pull(sample_file)
  sce_droplets <- read10xCounts(raw_file)
  rownames(sce_droplets) <- make.unique(rowData(sce_droplets)$Symbol)

  # calculate contamination fraction
  soup_channel <-
    SoupChannel(tod = counts(sce_droplets), toc = counts(sce_cells)) %>%
    setClusters(
      nb_metadata %>%
        select(cell, cluster_50) %>%
        semi_join(tibble(cell = rownames(colData(sce_cells))), by = "cell") %>%
        deframe() %>%
        as.character()
    ) %>%
    setDR(
      nb_metadata %>%
        select(cell, umap_1_monocle, umap_2_monocle) %>%
        semi_join(tibble(cell = rownames(colData(sce_cells))), by = "cell") %>%
        column_to_rownames("cell"),
      reductName = "UMAP"
    ) %>% 
     autoEstCont()
  
  # return adjusted counts
  list(
    soup_channel = soup_channel,
    adjusted_counts = adjustCounts(soup_channel)
  )
}

soupx_results <- map(samples$order, remove_soup)

soupx_counts <-
  soupx_results %>%
  map(~.x$adjusted_counts) %>%
  reduce(cbind)

waldo::compare(rownames(soupx_counts), rownames(nb))
waldo::compare(colnames(soupx_counts), colnames(nb))

assays(nb)$soupx_counts <- soupx_counts



# Plots -------------------------------------------------------------------

plot_change_map <- function(sce, genes, filename = NULL) {
  coordinates <- reducedDim(sce, "UMAP")
  
  sce_filtered <- sce[genes, rownames(coordinates), drop = FALSE]
  
  old <-
    counts(sce_filtered) %>%
    colSums()
  new <-
    assay(sce_filtered, "soupx_counts") %>%
    colSums()
  rel_change <- (old - new) / old
  
  p <- 
    coordinates %>% 
    magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>% 
    as_tibble() %>% 
    bind_cols(rel_change = rel_change) %>% 
    arrange(!is.na(rel_change), rel_change) %>% 
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = rel_change), size = .1) +
    annotate(
      "text_npc", npcx = 0.5, npcy = 1,
      label = genes, fontface = "bold", hjust = 0.5
    ) +
    scale_color_distiller(
      "relative\nchange",
      palette = "YlOrRd",
      direction = 1,
      na.value = "gray80",
      limits = c(0, 1)
    ) +
    coord_fixed() +
    theme_classic() +
    theme(panel.grid = element_blank())
  
  ggsave_default(filename)
  p
}

# plot_change_map(nb, "IL32")


map(
  c("HBB", "IGKC", "LYZ", "GNLY", "GAL", "IL32"),
  plot_change_map,
  sce = nb
) %>% 
  wrap_plots(guides = "collect")
ggsave_default("ambiance/relative_change", width = 420, height = 297)


tibble(
  gene = rownames(nb),
  n_cells_orig = rowSums(counts(nb) > 0),
  n_cells_corr = rowSums(assay(nb, "soupx_counts") > 0)
) %>% 
  mutate(
    abs_change = n_cells_orig - n_cells_corr,
    rel_change = abs_change / n_cells_orig
  ) %>% 
  arrange(desc(abs_change)) %>%
  select(!rel_change) %>% 
  column_to_rownames("gene") %>%
  head(20) %>% 
  gridExtra::tableGrob() %>% 
  wrap_plots()
ggsave_default("ambiance/table")



# Save results ------------------------------------------------------------

nb %>% saveRDS(out_file)
