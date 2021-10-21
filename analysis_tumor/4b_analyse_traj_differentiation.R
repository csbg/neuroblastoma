source("analysis_tumor/packages.R")



# Import data -------------------------------------------------------------

cds <- readRDS("analysis_tumor/data_generated/cds.rds")
DEG <- readRDS("analysis_tumor/data_generated/deg.rds")



# Step 1 ------------------------------------------------------------------

# assign pseudotime to cluster 2 principal points

set.seed(0)
cds_trajectory <- learn_graph(
  cds,
  use_partition = TRUE,
  verbose = TRUE,
  close_loop = TRUE,
  learn_graph_control = list(
    rann.k = 20, 
    minimal_branch_len = 10
  )
)

cds_trajectory <- order_cells(
  cds_trajectory,
  root_pr_nodes = c("Y_37", "Y_109", "Y_99", "Y_84"),
  reduction_method = "UMAP"
)

plot_cells(cds_trajectory, color_cells_by = "pseudotime")



# Step 1a (optional) ------------------------------------------------------

## Moran's I test ----

# find differentially expressed genes along trajectory
DEG_trajectory <- graph_test(
  cds_trajectory, 
  neighbor_graph = "principal_graph",
  reduction_method = "UMAP",
  k = 20,
  verbose = TRUE
)

# remove NAs & filter according to QC criteria
DEG_trajectory <- 
  DEG_trajectory %>% 
  na.omit() %>% 
  filter(q_value <= 0.05, morans_I > 0.1, q_value != 0) %>%
  as_tibble()


## Compare with detected clusterwise DEG ----

Overlap_DEG <-
  DEG_trajectory %>%
  filter(gene_short_name %in% DEG$gene_short_name)
  
Overlap_DEG_rev <-
  DEG %>%
  filter(gene_short_name %in% DEG_trajectory$gene_short_name)

Top_DEG <- list(
  "Overlap_Trajectory_DEG" =
    Overlap_DEG %>% 
    slice_max(n = 56, order_by = morans_I) %>%
    pull(gene_short_name) %>% 
    unique(), 
  "Overlap_DEG_Trajectory" = 
    Overlap_DEG_rev %>% 
    slice_max(n = 63, order_by = estimate) %>% 
    pull(gene_short_name) %>% 
    unique()
  )

# plot_cells(
#   cds_trajectory,
#   genes = Top_DEG[[2]],
#   cell_size = 0.5,
#   label_branch_points = FALSE,
#   label_roots = FALSE,
#   label_leaves = FALSE
# )

list(
  "DE along trajectory" = DEG_trajectory, 
  "Overlap_Trajectory_DEG" = Overlap_DEG, 
  "Overlap_DEG_Trajectory" = Overlap_DEG_rev
) %>% 
  saveRDS("analysis_tumor/data_generated/deg_trajectory.rds")



# Step 3: CytoTRACE -------------------------------------------------------

cytoMA <- CytoTRACE(
  assays(cds)$raw_counts %>% as.matrix(),
  batch = table(rownames(cds@colData), cds@colData$sample),
  enableFast = FALSE
)

# assign results to CDS trajectory metadata
cds_trajectory@colData$cytoTRACE <-
  cytoMA[1]$CytoTRACE %>% 
  as.numeric()



# Save results ------------------------------------------------------------

saveRDS(cds_trajectory, "analysis_tumor/data_generated/cds_trajectory.rds")
saveRDS(cytoMA, "analysis_tumor/data_generated/cytotrace.rds")
