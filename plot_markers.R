library(Seurat)
library(tidyverse)
library(fs)
library(patchwork)
library(viridis)

ggsave_default <- function(filename, width = 297, height = 210,
                           crop = TRUE, ...) {
  if (is.null(filename))
    return()
  
  filename <- str_glue("plots/{filename}.png")
  filename %>% path_dir() %>% dir_create()
  
  ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
         width = width, height = height, ...)
  
  if (crop) knitr::plot_crop(filename)
  
  invisible(filename)
}



# Load data ---------------------------------------------------------------

nb <- readRDS("data_generated/all_datasets_current/nb_assay_SCT.rds")
Idents(nb) <- "integrated_snn_res.0.5"

nb@meta.data <- 
  nb@meta.data %>% 
  rownames_to_column("cell") %>% 
  left_join(
    read_csv(
      "data_raw/metadata/sample_groups.csv",
      col_types = "ccf",
      comment = "#"
    ) %>%
      distinct(sample, group) %>% 
      mutate(group = fct_relevel(group, "I", "II", "III", "IV")),
    by = "sample"
  ) %>%
  column_to_rownames("cell")



# Mitochondrial genes -----------------------------------------------------

FeaturePlot(nb, "percent.mt", coord.fixed = TRUE) +
  scale_color_viridis(option = "viridis")
ggsave_default("qc_mtgene_umap", width = 200, height = 200)

ggplot(nb@meta.data, aes(integrated_snn_res.0.5, percent.mt)) +
  geom_violin(aes(fill = integrated_snn_res.0.5), show.legend = FALSE) +
  geom_jitter(alpha = .1) +
  xlab("Cluster") +
  ylab("% mitochondrial genes") +
  theme_classic()
ggsave_default("qc_mtgene_per_cluster", width = 200, height = 150)



# Canonical cell type markers ---------------------------------------------

markers <- read_csv("data_raw/metadata/cell_markers.csv", comment = "#")

Idents(nb) <- 
  fct_relevel(
    Idents(nb),
    "0", "14", "7", "21",               # T cells
    "2", "8",                           # NK cells
    "1", "6", "13", "20",               # B cells
    "3", "11", "19", "10", "12", "18",  # myeloid
    "15",                               # pDC
    "5", "9",                           # NB
    "17",                               # erythroblast
    "16",                               # CMP
    "4"                                 # other
  )

DotPlot(nb, features = rev(markers$gene)) +
  coord_flip() +
  scale_color_viridis(option = "inferno", direction = -1)
ggsave_default("markers/canonical_markers")



# NB markers --------------------------------------------------------------

nb_markers <-
  markers %>%
  filter(cell_type == "neuroblastoma") %>%
  pull(gene)

FeaturePlot(
  nb,
  nb_markers,
  min.cutoff = "q5",
  max.cutoff = "q95",
  coord.fixed = TRUE,
  ncol = 4,
  order = TRUE
) &
  scale_color_viridis(option = "cividis")
ggsave_default("markers/nb_markers", width = 420, height = 297)


# this only works on the cluster
walk(
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
    ggsave_default(str_glue("markers/nb_dotplot_{g}"),
                   height = 100, width = 200)
  }
)

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
ggsave_default(str_glue("markers/nb_dotplot"),
               height = 200, width = 400)