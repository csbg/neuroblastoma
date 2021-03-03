library(Seurat)
library(SeuratObject)
library(tidyverse)
library(fs)
library(patchwork)
library(viridis)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/assay_SCT_seurat.rds")
nb_data <- readRDS("data_generated/metadata.rds")

nb@meta.data <-
  nb_data %>%
  column_to_rownames("cell")



# Canonical cell type markers ---------------------------------------------

markers <- read_csv("metadata/cell_markers.csv", comment = "#")


# monocle clusters
Idents(nb) <- 
  fct_relevel(
    nb@meta.data$cluster_50,
    # T cells
    "3", "5", "18",
    
    # NK cells
    "4", "6",
    
    # B cells
    "2", "9", "12", "17", "19", "21",
    
    # monocyte
    "1", "15", "16", "22",
    
    # pDC
    "14",
    
    # erythroblast
    "13",
    
    # bone marrow
    "10", "11", 
    
    # other
    "7", "20",
    
    # NB
    "8"
  )

DotPlot(nb, features = rev(markers$gene)) +
  coord_flip() +
  scale_color_viridis(option = "inferno", direction = -1)
ggsave_default("monocle/markers/canonical_markers_clusters",
               height = 297, width = 210)
ggsave_default("monocle/markers/canonical_markers_clusters",
               type = "pdf", height = 297, width = 210, crop = FALSE)



# NB markers --------------------------------------------------------------

nb_markers <-
  markers %>%
  filter(cell_type == "neuroblastoma") %>%
  pull(gene)

# add monocle UMAP coordinates to the Seurat object
nb@reductions$umap_monocle <- CreateDimReducObject(
  embeddings = nb_data %>% 
    column_to_rownames("cell") %>% 
    select(UMAP_1 = umap_1_monocle, UMAP_2 = umap_2_monocle) %>% 
    as.matrix(),
  assay = "monocle_align_cds",
  global = TRUE,
  key = "UMAPMON_"
)

FeaturePlot(
  nb,
  nb_markers,
  reduction = "umap_monocle",
  min.cutoff = "q5",
  max.cutoff = "q95",
  coord.fixed = TRUE,
  ncol = 4,
  order = TRUE
) &
  scale_color_viridis(option = "cividis")
ggsave_default("monocle/markers/nb_markers", width = 420, height = 297)


map(
  levels(nb@meta.data$group),
  function(g) {
    message("Plotting dotplot for group ", g)
    group_cells <-
      nb@meta.data %>%
      as_tibble(rownames = "cell") %>%
      filter(group == g) %>%
      pull(cell)
    DotPlot(nb[, group_cells], features = nb_markers) +
      coord_flip() +
      scale_color_viridis(option = "inferno", direction = -1)
  }
) %>% 
  wrap_plots()
ggsave_default(str_glue("monocle/markers/nb_dotplot_groupwise"),
               height = 200, width = 400)
