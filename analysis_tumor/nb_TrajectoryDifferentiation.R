source("nb_Packages.R")

cds <- readRDS(paste0(path,"/nb_CDS.rds")) # change link
DEG <- readRDS(paste0(path,"/nb_DEG.rds")) # change link

################################################################################
# Step 1: learn Trajectory on monocle CDS & 
#           assign pseudotime to cluster 2 principal points

set.seed(0)
cds_trajectory <- learn_graph(
                    cds,
                    use_partition = T,
                    verbose = TRUE,
                    close_loop = TRUE,
                    learn_graph_control = 
                      list(rann.k = 20, 
                           minimal_branch_len = 10
                           )
                    )


cds_trajectory <- order_cells(
                    cds_trajectory,
                    root_pr_nodes = c("Y_37","Y_109","Y_99","Y_84"),
                    reduction_method = "UMAP"
                    )

plot_cells(cds_trajectory, 
           label_principal_points = T, 
           color_cells_by = "pseudotime")

################################################################################
##Optional: Moranś I test to find differential expressed genes along trajectory

DEG_trajectory <- graph_test(
                    cds_trajectory, 
                    neighbor_graph = "principal_graph",
                    reduction_method = "UMAP",
                    k = 20,
                    cores = 3,
                    verbose = TRUE
                    )

#remove NA's & filter according to QC criteria
DEG_trajectory <- 
  DEG_trajectory %>% 
  na.omit() %>% 
  filter(
    q_value <= 0.05 & morans_I > 0.1 & q_value != 0
    ) %>%
  as_tibble()
################################################################################
# optional
# Compare with detected clusterwise DEG
Overlap_DEG <-
  DEG_trajectory %>%
  filter(
    gene_short_name %in% DEG$gene_short_name
        )
  
Overlap_DEG_rev <-
  DEG %>%
  filter(
    gene_short_name %in% DEG_trajectory$gene_short_name
  )

Top_DEG <- list(
  "Overlap_Trajectory_DEG" =
    unique((Overlap_DEG %>% 
              slice_max(
                order_by = morans_I, 
                n = 56
                )
     )$gene_short_name), 
  "Overlap_DEG_Trajectory" = 
    unique((Overlap_DEG_rev %>% 
              slice_max(
                order_by = estimate, 
                n = 63
                )
     )$gene_short_name)
  )

#plot_cells(cds_trajectory,
#           genes = Top_DEG[[2]],
#           cell_size = 0.5,
#           label_branch_points = F,
#           label_roots = F,
#           label_leaves = F
#           )

saveRDS(list("DE along trajectory" = DEG_trajectory, 
             "Overlap_Trajectory_DEG" = Overlap_DEG, 
             "Overlap_DEG_Trajectory" = Overlap_DEG_rev
             ),
        file = paste0(path, "/nb_DEG_trajectory.rds")
        )

################################################################################
# Step 3:
# CytoTRACE for assigning differentiation potential
################################################################################

cytoMA <- 
  CytoTRACE(assays(cds)$raw_counts %>% as.matrix(),
            enableFast = F,
            batch = table(rownames(cds@colData),cds@colData$sample)
            )
# Assign results to CDS-Trajectory metadata
cds_trajectory@colData$cytoTRACE <- 
  cytoMA[1]$CytoTRACE %>% 
  as.numeric()

saveRDS(cds_trajectory,
        file = paste0(path,"/nb_CDS_trajectory.rds"))

# integrated cytoTRACE with multiple batch correction
iCytoTRACE()

pheno <- 
  cds@colData %>%
  as_tibble() %>%
  select(cluster) 

pheno <- 
  pheno$cluster %>% 
  as.character()

names(pheno) <- rownames(cds@colData)

#plotCytoTRACE(cytoMA,
#              phenotype = pheno, 
#              emb = reducedDim(cds, "UMAP")
#              )

#plotCytoGenes(cytoMA,numOfGenes = 50)

saveRDS(cytoMA,
        file = paste0(path, "/nb_CytoTRACE.rds")) # change link
