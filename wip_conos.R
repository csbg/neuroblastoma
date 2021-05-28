# data integration using the conos package

library(conos)
library(pagoda2)
library(SeuratObject)
library(tidyverse)
source("common_functions.R")



# Preprocessing -----------------------------------------------------------

cell_types <-
  read_csv("data_wip/cellid_predictions.csv", col_types = "cfd") %>% 
  filter(log_q > 2) %>%
  select(!log_q) %>%
  deframe()

nb <- readRDS("data_generated/rna_qcpassed.rds")
nb_names <- 
  map_chr(
    nb,
    ~str_c(.@meta.data$group[1], .@meta.data$sample[1], sep = "_")
  ) %>%
  make.unique()
  
nb <- set_names(nb, nb_names)
nb <- nb[sort(names(nb))]

nb <- map(
  nb,
  ~basicP2proc(.$RNA@counts, n.cores = 6)
)



# Conos analysis ----------------------------------------------------------

con <- Conos$new(nb, n.cores = 6)

con$buildGraph(
  k = 30,
  k.self = 5,
  space = "CPCA",
  score.component.variance = TRUE,
  verbose = TRUE
)

saveRDS(con, "data_wip/rna_conos_1_cpca.rds")
# con <- readRDS("data_wip/rna_conos_1_cpca.rds")


# clustering
con$findCommunities(method = leiden.community, resolution = 3)
# con$findCommunities(method = igraph::walktrap.community, steps = 7)
# fc <- greedyModularityCut(con$clusters$walktrap$result, 40)

# dimensional reduction
con$embedGraph(method = "largeVis")
con$embedGraph(method = "UMAP", min.dist = 0.01,
               spread = 15, min.prob.lower = 1e-3)

saveRDS(con, "data_wip/rna_conos_2_umap.rds")
# con <- readRDS("data_wip/rna_conos_2_umap.rds")



# Plots -------------------------------------------------------------------

# individual clusters
con$plotPanel(clustering = "multilevel", use.local.clusters = TRUE)
ggsave_default("wip/conos_clusters_individual")

# common clusters
con$plotPanel()
ggsave_default("wip/conos_clusters_common")

plotClusterBarplots(con)
ggsave_default("wip/conos_clusters_barplot", width = 210, height = 297)

## largeVis
con$plotGraph(embedding = "largeVis", size = 0.1) +
  coord_fixed()
ggsave_default("wip/conos_largevis_clusters")

con$plotGraph(embedding = "largeVis", size = 0.1, gene = "PHOX2A") +
  coord_fixed()
ggsave_default("wip/conos_largevis_phox2a")

con$plotGraph(
  embedding = "largeVis",
  alpha = 0.1,
  groups = cell_types,
  plot.na = -1
) +
  coord_fixed()
ggsave_default("wip/conos_largevis_celltypes")

con$plotPanel(embedding = "largeVis", use.common.embedding = TRUE,
              size = 0.1, font.size = 2, title.size = 3)
ggsave_default("wip/conos_largevis_panel")

con$plotPanel(embedding = "largeVis", use.common.embedding = TRUE,
              size = 0.1, font.size = 2, title.size = 3, gene = "PHOX2A")
ggsave_default("wip/conos_largevis_panel_phox2a")



## UMAP
con$plotGraph(embedding = "UMAP", size = 0.1) +
  coord_fixed()
ggsave_default("wip/conos_umap_clusters")

con$plotGraph(embedding = "UMAP", size = 0.1, gene = "PHOX2A") +
  coord_fixed()
ggsave_default("wip/conos_umap_phox2a")

con$plotGraph(
  embedding = "UMAP",
  alpha = 0.1,
  groups = cell_types,
  plot.na = -1
) +
  coord_fixed()
ggsave_default("wip/conos_umap_celltypes")

con$plotPanel(embedding = "UMAP", use.common.embedding = TRUE,
              size = 0.1, font.size = 2, title.size = 3)
ggsave_default("wip/conos_umap_panel")

con$plotPanel(embedding = "UMAP", use.common.embedding = TRUE,
              size = 0.1, font.size = 2, title.size = 3, gene = "PHOX2A")
ggsave_default("wip/conos_umap_panel_phox2a")
