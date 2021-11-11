library(Seurat)
library(monocle3)
library(scuttle)
library(muscat)
library(umap)
library(ProjecTILs)
library(tidyverse)
library(fs)
library(ComplexHeatmap)
library(scico)
source("common_functions.R")
source("styling.R")

ADRENAL_MEDULLA_COLORS <- c(
  "late SCPs" = "#921813",
  "SCPs" = "#be202e",
  "cycling SCPs" = "#be6867",
  "Bridge" = "#ef845f",
  "connecting Chromaffin cells" = "#33ac75",
  "Chromaffin cells" = "#006e3b",
  "late Chromaffin cells" = "#686c58",
  "cycling Neuroblasts" = "#aacedc",
  "Neuroblasts" = "#244a94",
  "late Neuroblasts" = "#303d63"
)
ADRENAL_MEDULLA_COLORS <- 
  ADRENAL_MEDULLA_COLORS %>% 
  colorspace::lighten(0.3) %>% 
  set_names(names(ADRENAL_MEDULLA_COLORS))

files_dong <- read_csv("data_wip/files_dong.csv", comment = "#")

  

# Prepare Dong data -------------------------------------------------------

# metadata_dong <-
#   read_csv("data_raw/GSE137804/GSE137804_tumor_dataset_annotation.csv") %>%
#   separate(cellname, into = c("sample", "cell"), sep = "_")
# 
# data_dong <- pmap(
#   files_dong %>% filter(high_risk),
#   function(file, tumor_id, bad_header, ...) {
#     if (bad_header) {
#       sce <- read_tsv(
#         path_join(c("data_raw", "GSE137804", file)),
#         skip = 1,
#         col_names =
#           read_tsv(path_join(c("data_raw", "GSE137804", file)), n_max = 0) %>%
#           colnames() %>%
#           prepend("Symbol")
#       )
#     } else {
#       sce <-
#         read_tsv(
#           path_join(c("data_raw", "GSE137804", file)),
#           col_types = cols(
#             Gene_ID = "c",
#             Symbol = "c",
#             .default = col_integer()
#           )
#         ) %>%
#         select(!Gene_ID)
#     }
# 
#     sce <-
#       sce %>%
#       mutate(Symbol = make.unique(Symbol)) %>%
#       column_to_rownames("Symbol") %>%
#       as.matrix() %>%
#       CreateSeuratObject(project = tumor_id)
# 
#     sce@meta.data <-
#       sce@meta.data %>%
#       as_tibble(rownames = "cell") %>%
#       rename(sample = orig.ident) %>%
#       left_join(metadata_dong, by = c("cell", "sample")) %>%
#       column_to_rownames("cell")
# 
#     sce %>%
#       subset(subset = celltype == "tumor")
#   }
# )
# 
# data_dong_names <- map_chr(data_dong, Project)
# data_dong %>%
#   set_names(data_dong_names) %>%
#   saveRDS("data_wip/tumor_data_dong.rds")



# Prepare reference atlas -------------------------------------------------

# adr <-
#   readRDS("~/Desktop/adrenal_medulla_Seurat.RDS") %>% 
#   RunUMAP(dims = 1:30, return.model = TRUE)
# adr
# 
# set.seed(1234)
# 
# varfeat <- adr$RNA@var.features
# refdata <- data.frame(t(adr$RNA@data[varfeat, ]))
# refdata <- refdata[, sort(colnames(refdata))]
# 
# ref.pca <- prcomp(refdata, rank. = 50, scale. = TRUE, center = TRUE, retx = TRUE)
# ref.umap <- umap(
#   ref.pca$x[, 1:10],
#   n_neighbors = 30,
#   min_dist = 0.3,
#   metric = "cosine",
#   random_state = 1234,
#   transform_state = 1234,
#   verbose = TRUE
# )
# 
# colnames(ref.umap$layout) <- c("UMAP_1", "UMAP_2")
# ref.umap
# 
# 
# 
# adr@reductions$umap@cell.embeddings <- ref.umap$layout
# adr@reductions$pca@cell.embeddings <- ref.pca$x
# adr@reductions$pca@feature.loadings <- ref.pca$rotation
# colnames(adr@reductions$pca@cell.embeddings) <-
#   gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl = TRUE)
# colnames(adr@reductions$pca@feature.loadings) <-
#   gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl = TRUE)
# adr@misc$pca_object <- ref.pca
# adr@misc$umap_object <- ref.umap
# adr@misc$projecTILs <- "adrenal_medulla_atlas"
# 
# adr %>% 
#   RenameAssays(RNA = "integrated") %>% 
#   saveRDS("data_wip/projectils/reference_atlas_adr.rds")



# ProjecTILs analysis -----------------------------------------------------

# evil hack to prevent conversion of human to mouse gene names
convert.orthologs.id <- function(x, ...) x
environment(convert.orthologs.id) <- asNamespace("ProjecTILs")
assignInNamespace("convert.orthologs", convert.orthologs.id, ns = "ProjecTILs")

ref <- load.reference.map("data_wip/projectils/reference_atlas_adr.rds")
ref$functional.cluster <- Idents(ref)

# this is the preferred approach, but alas, memory is scarce
# query_projected <- make.projection(
#   query = tumor_data,
#   ref = ref,
#   filter.cells = FALSE
# )

# (a) Dong
tumor_data <- readRDS("data_wip/tumor_data_dong.rds")

# process Dong datasets 8–10 individually due to memory issues
tumor_data <- tumor_data[8]   # T162
# tumor_data <- tumor_data[9]   # T200
# tumor_data <- tumor_data[10]  # T230
# gc()

# (b) Jansky
# tumor_data <-
#   readRDS("data_wip/tumor_data_jansky.rds") %>%
#   subset(subset = anno_new == "Tumor cells") %>% 
#   SplitObject("patientID")

# (c) own data
# tumor_data <-
#   readRDS("data_generated/rna_decontaminated.rds") %>%
#   counts() %>%
#   CreateSeuratObject(
#     project = "NB",
#     meta.data =
#       readRDS("data_generated/metadata.rds") %>%
#       column_to_rownames("cell")
#   ) %>%
#   subset(subset = cellont_abbr == "NB") %>%
#   SplitObject("sample")


# (a-c) common analysis function
iwalk(
  tumor_data,
  function(query_data, sample) {
    print(sample)
    gc()
    
    res <-
      make.projection(
        query = query_data,
        ref = ref,
        filter.cells = FALSE
      ) %>% 
      cellstate.predict(ref = ref)
    
    list(
      Embeddings(res, "umap") %>%
        as_tibble(rownames = "cell"),
      res@meta.data %>% 
        select(
          predicted_type = functional.cluster,
          confidence_score = functional.cluster.conf
        ) %>%
        as_tibble(rownames = "cell")
    ) %>% 
      reduce(left_join, by = "cell") %>% 
      write_csv(str_glue("data_wip/projectils/projectils_results_{sample}.csv"))
  }
)
q()


# Plots -------------------------------------------------------------------

projectils_data <-
  dir_ls("data_wip/projectils", regex = "results") %>% 
  map_dfr(read_csv, .id = "filename") %>% 
  extract(filename, into = "sample", regex = "results_(.+)\\.csv") %>% 
  mutate(
    sample = rename_patients(sample),
    group = if_else(str_ends(sample, ".M."), "Tm", str_sub(sample, 1, 1))
  )

Embeddings(ref, "umap") %>%
  as_tibble() %>% 
  mutate(cell_type = Idents(ref)) %>% 
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(
    aes(color = cell_type),
    size = 0.1
  ) +
  geom_density_2d(
    data = projectils_data %>% filter(confidence_score >= 0.5),
    contour_var = "ndensity",
    color = "black",
    na.rm = TRUE,
    size = 0.25
  ) +
  scale_color_manual(
    values = ADRENAL_MEDULLA_COLORS,
    guide = guide_legend(override.aes = list(size = 3))
  ) +
  coord_fixed() +
  facet_wrap(vars(group), nrow = 2) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("comparison/projection")


projectils_data %>%
  filter(confidence_score >= 0.5) %>%
  mutate(
    predicted_type = factor(
      predicted_type,
      levels = levels(ref$functional.cluster)
    )
  ) %>%
  group_by(sample, predicted_type) %>%
  summarise(n = n()) %>% 
  mutate(n_rel = n / sum(n) * 100, group = str_sub(sample, 1, 1)) %>% 
  ggplot(aes(sample, n_rel)) +
  geom_col(aes(fill = group)) +
  geom_hline(yintercept = 0) +
  xlab("sample") +
  ylab("relative abundance (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(GROUP_COLORS, "T" = "#433447")) +
  facet_wrap(vars(predicted_type), ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("comparison/predicted_cell_types", height = 210, width = 297)



# Correlation -------------------------------------------------------------

tumor_data_merged <-
  readRDS("data_wip/tumor_data.rds") %>% 
  {merge(.[[1]], .[-1])} %>% 
  as.SingleCellExperiment() %>% 
  logNormCounts()

colData(tumor_data_merged)$sample <- 
  colData(tumor_data_merged)$sample %>%
  rename_patients()

pb_tumor <-
  tumor_data_merged %>% 
  aggregateData(
    assay = "logcounts",
    fun = "mean",
    by = "sample"
  )

hvgs <-
  tumor_data_merged %>%
  scran::modelGeneVar() %>% 
  scran::getTopHVGs()

corr_mat <-
  assay(pb_tumor, 1) %>%
  magrittr::extract(hvgs, ) %>%
  cor(use = "pairwise.complete.obs")

distance <- as.dist(1 - corr_mat)

group_names <-
  colnames(corr_mat) %>%
  map_chr(str_sub, 1, 1)

mycn_status <- 
  files_dong %>%
  select(tumor_id, mycn_amplified) %>%
  deframe() %>% 
  append(
    PATIENT_ORDER %>% 
      str_detect("M") %>% 
      set_names(PATIENT_ORDER)
  ) %>% 
  magrittr::extract(colnames(corr_mat)) %>% 
  if_else("amplified", "normal")

p <- Heatmap(
  corr_mat,
  col = circlize::colorRamp2(
    seq(min(corr_mat), 1, length.out = 9),
    scico(9, palette = "davos", direction = -1),
  ),
  name = "correlation of\npseudobulk\nexpression",
  heatmap_legend_param = list(
    at = c(round(min(corr_mat), 2), 0.9, 1)
  ),
  
  clustering_distance_rows = distance,
  clustering_distance_columns = distance,
  
  width = unit(80, "mm"),
  height = unit(80, "mm"),
  
  show_column_dend = FALSE,
  
  left_annotation = rowAnnotation(
    group = group_names,
    mycn = mycn_status,
    col = list(
      group = c(GROUP_COLORS, "T" = "#433447"),
      mycn = c("normal" = "gray90", "amplified" = "#d35f5f")
    ),
    show_annotation_name = FALSE,
    show_legend = TRUE,
    annotation_legend_param = list(
      group = list(
        title = "group"
      ),
      mycn = list(
        title = "MYCN status"
      )
    )
  ),
)
p
ggsave_default("comparison/correlation", plot = p)
