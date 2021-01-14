library(monocle3)
library(tidyverse)

cds <- new_cell_data_set(readRDS("test_data/cao_l2_expression.rds"),
                         cell_metadata = readRDS("test_data/cao_l2_colData.rds"),
                         gene_metadata = readRDS("test_data/cao_l2_rowData.rds"))
cds

cds <- preprocess_cds(cds, num_dim = 100, verbose = TRUE)
cds
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, verbose = TRUE)
cds <- reduce_dimension(cds, reduction_method = "tSNE", verbose = TRUE)
cds

plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cao_cell_type", group_label_size = 4)

plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

cds <- cluster_cells(cds, verbose = TRUE, resolution = 1e-5)
plot_cells(cds)
plot_cells(cds, group_cells_by = "partition")


marker_test_res <- top_markers(cds, group_cells_by = "partition", 
                               reference_cells = 1000, cores = 2)


top_specific_markers <-
  marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)
