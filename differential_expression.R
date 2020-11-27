# WORK IN PROGRESS

library(tidyverse)
library(Seurat)

nb <- readRDS("rds/nb_integrated_5.rds")
nb@meta.data$sample <- str_match(nb@meta.data$sample, ".*/(.*)")[,2]
nb <- subset(nb, subset = sample != "R1_GNB_2014_0102")

nb@meta.data %>% head()

nb <- 
  nb %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = .5)

DimPlot(nb, group.by = "sample")
DimPlot(nb, split.by = "sample", ncol = 2, label = TRUE)
DimPlot(nb)


cluster11_conserved_markers <- FindConservedMarkers(nb, ident.1 = 11, grouping.var = "sample")
FindMarkers()

FeaturePlot(nb, features = c("CDHR3", "VPREB3", "CD79B", "CD24"))



nb_data <- inner_join(
  nb@meta.data %>% as_tibble(rownames = "cell"),
  nb$tsne@cell.embeddings %>% as_tibble(rownames = "cell"),
  by = "cell"
)

  
nb_data %>% 
  ggplot(aes(tSNE_1, tSNE_2)) +
  geom_point() +
  facet_wrap(vars(sample))
