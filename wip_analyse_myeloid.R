library(monocle3)
library(scuttle)
library(tidyverse)
library(scico)
library(patchwork)
source("common_functions.R")
source("styling.R")



# Load data ---------------------------------------------------------------

nb_metadada <- readRDS("data_generated/metadata.rds")

myeloid_barcodes <- 
  nb_metadada %>% 
  filter(cellont_abbr == "M") %>% 
  pull(cell)

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  magrittr::extract(, myeloid_barcodes) %>% 
  logNormCounts(assay.type = "soupx_counts")

markers <- read_csv("metadata/myeloid_markers.csv")
# markers$gene %>% setdiff(rownames(cds_my))


# Analysis ----------------------------------------------------------------

set.seed(42)
cds_my <- 
  nb %>% 
  preprocess_cds(verbose = TRUE) %>% 
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE) %>% 
  align_cds(alignment_group = "sample", verbose = TRUE) %>% 
  reduce_dimension(
    reduction_method = "UMAP",
    preprocess_method = "Aligned",
    verbose = TRUE
  ) %>%
  cluster_cells(k = 20, random_seed = 42, verbose = TRUE)

my_metadata <-
  list(
    tibble(cell = rownames(colData(cds_my))),
    nb_metadada %>%
      select(cell, sample, group, cellont_cluster) %>% 
      mutate(group = rename_groups(group), sample = rename_patients(sample)),
    reducedDim(cds_my, "UMAP") %>%
      magrittr::set_colnames(c("UMAP1", "UMAP2")) %>% 
      as_tibble(rownames = "cell"),
    clusters(cds_my) %>% 
      enframe(name = "cell", value = "subcluster")
  ) %>% 
  reduce(left_join, by = "cell")



# Plots -------------------------------------------------------------------

my_metadata %>% 
  count(subcluster) %>% 
  mutate(n_rel = n / sum(n) * 100)



# UMAPs ----

cluster_labels <-
  my_metadata %>% 
  group_by(subcluster) %>% 
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

ggplot(my_metadata, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = subcluster), size = 0.1, show.legend = FALSE) +
  geom_text(data = cluster_labels, aes(label = subcluster)) +
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("myeloid/umap")


ggplot(my_metadata, aes(UMAP1, UMAP2)) +
  geom_hex(bins = 100) +
  scale_fill_scico(palette = "nuuk") +
  coord_fixed() +
  facet_wrap(vars(group)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("myeloid/clusters_groups")


my_metadata %>% 
  group_by(subcluster) %>% 
  count(group) %>% 
  mutate(n = n / sum(n) * 100) %>% 
  filter(as.numeric(subcluster) <= 11) %>% 
  ggplot(aes("", n, fill = group)) +
  geom_col() +
  scale_x_discrete(NULL) +
  scale_fill_manual(values = GROUP_COLORS) +
  coord_polar("y") +
  facet_wrap(vars(subcluster), nrow = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("myeloid/abundances_groups")


my_metadata %>% 
  group_by(subcluster) %>% 
  count(sample) %>% 
  mutate(n = n / sum(n) * 100) %>% 
  filter(as.numeric(subcluster) <= 11) %>% 
  ggplot(aes("", n, fill = sample)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = PATIENT_COLORS) +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  facet_wrap(vars(subcluster), nrow = 2) +
  theme_nb() +
  theme(
    legend.key.width = unit(2, "mm"),
    legend.key.height = unit(2, "mm"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("myeloid/abundances_patients", height = 50, width = 150)



## Markers ----

plot_cell_type <- function(t) {
  plot_dots(
    logcounts(cds_my),
    markers %>% 
      filter(cell_type == t) %>% 
      pull(gene),
    clusters(cds_my)
  ) +
    scale_y_discrete(t) +
    theme(legend.position = "none")
}



unique(markers$cell_type) %>% 
  map(plot_cell_type) %>% 
  wrap_plots(ncol = 1)
ggsave_default("myeloid/markers", height = 700)


walk(
  unique(markers$cell_type),
  function(t) {
    plot_cell_type(t)
    ggsave_default(str_glue("myeloid/markers_{t}"), width = 100)
  }
)
