# Cell type assignment via SingleR for an integrated Seurat dataset
#
# Adds metadata on clustering to a Seurat object, as well as two columns
# SingleR_cells and SingleR_clusters, containing predicted cell types on
# the individual cell and cluster level, respectively.

library(magrittr)
library(Seurat)
library(SingleR)
library(celldex)


# small test datasets
seurat_infile <- "data_generated/5_datasets/nb_integrated.rds"
seurat_outfile <- "data_generated/5_datasets/nb_integrated_singler.rds"

# large complete datasets
# seurat_infile <- "data_generated/all_datasets/nb_integrated.rds
# seurat_outfile <- "data_generated/all_datasets/nb_integrated_singler.rds"



# Load data ---------------------------------------------------------------

nb <- readRDS(seurat_infile)

nb <- 
  nb %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = .5)

# for testing: choose only a few cells
# nb <- subset(nb, cells = sample(ncol(nb), 2000))

reference_cell_types <- HumanPrimaryCellAtlasData()



# Predict types of individual cells ---------------------------------------

predicted_cell_types <- SingleR(
  GetAssayData(nb$RNA),
  reference_cell_types,
  labels = reference_cell_types$label.main
)
dim(predicted_cell_types)


# Predict types of clusters -----------------------------------------------

predicted_cluster_types <- SingleR(
  GetAssayData(nb$RNA),
  reference_cell_types,
  labels = reference_cell_types$label.main,
  clusters = Idents(nb)
)
predicted_cluster_types


# Combine data ------------------------------------------------------------

nb <-
  nb %>%
  AddMetaData(predicted_cell_types$labels, "SingleR_cells") %>%
  AddMetaData(predicted_cluster_types[Idents(nb),]$labels, "SingleR_clusters")

# DimPlot(nb)
# DimPlot(nb, group.by = "SingleR_cells")
# DimPlot(nb, group.by = "SingleR_clusters")

saveRDS(nb, seurat_outfile)
