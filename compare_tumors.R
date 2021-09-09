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

rename_dong <- partial(
  rename_str_or_fct,
  nm = tribble(
    ~old,       ~new,
    "T162",  "T162 (M)",
    "T200",  "T200 (M)",
    "T230",  "T230 (M)",
  )
)

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

files_dong <- tribble(
  ~tumor_id, ~file, ~bad_header, ~high_risk, ~mycn_amplified,
  "T10",  "GSM4088774_T10_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T19",  "GSM4088775_T19_gene_cell_exprs_table.xls",  FALSE, FALSE, NA,
  "T27",  "GSM4088776_T27_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T34",  "GSM4088777_T34_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T40",  "GSM4088778_T40_gene_cell_exprs_table.xls",  FALSE, FALSE, NA,
  "T44",  "GSM4088779_T44_gene_cell_exprs_table.xls",  FALSE, FALSE, NA,
  "T69",  "GSM4088780_T69_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T71",  "GSM4088781_T71_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T75",  "GSM4088782_T75_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T92",  "GSM4088783_T92_gene_cell_exprs_table.xls",  FALSE, TRUE, FALSE,
  "T162", "GSM4654669_T162_gene_cell_exprs_table.xls", TRUE, TRUE, TRUE,
  "T175", "GSM4654670_T175_gene_cell_exprs_table.xls", FALSE, FALSE, NA,
  "T188", "GSM4654671_T188_gene_cell_exprs_table.xls", FALSE, FALSE, NA,
  "T200", "GSM4654672_T200_gene_cell_exprs_table.xls", TRUE, TRUE, TRUE,
  "T214", "GSM4654673_T214_gene_cell_exprs_table.xls", TRUE, TRUE, FALSE,
  "T230", "GSM4654674_T230_gene_cell_exprs_table.xls", TRUE, TRUE, TRUE,
)



# Prepare tumor data ------------------------------------------------------

# nb_metadata <- readRDS("data_generated/metadata.rds")
# 
# nb_tumor <-
#   readRDS("data_generated/rna_decontaminated.rds") %>%
#   counts() %>% 
#   CreateSeuratObject(
#     project = "NB",
#     meta.data = nb_metadata %>% column_to_rownames("cell")
#   ) %>% 
#   subset(subset = cellont_abbr == "NB") %>% 
#   SplitObject("sample")
# 
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
#   append(nb_tumor) %>% 
#   saveRDS("data_wip/tumor_data.rds")



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
#   saveRDS("data_wip/reference_atlas_adr.rds")




# Projectils --------------------------------------------------------------

## Analysis ----

# evil hack to prevent conversion of human to mouse gene names
convert.orthologs.id <- function(x, ...) x
environment(convert.orthologs.id) <- asNamespace("ProjecTILs")
assignInNamespace("convert.orthologs", convert.orthologs.id, ns = "ProjecTILs")

tumor_data <- readRDS("data_wip/tumor_data.rds")
ref <- load.reference.map("data_wip/reference_atlas_adr.rds")
ref$functional.cluster <- Idents(ref)

# this is the preferred approach, but alas, memory is scarce
# query_projected <- make.projection(
#   query = tumor_data,
#   ref = ref,
#   filter.cells = FALSE
# )

# process datasets 8–10 individually due to memory issues
# tumor_data <- list(
#   T162 = tumor_data$T162 %>%
#     magrittr::extract(, sample(colnames(.), size = ncol(.) / 2))
# )
# tumor_data <- tumor_data[9]   # T200
# tumor_data <- tumor_data[10]  # T230
# gc()

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
      write_csv(str_glue("data_wip/projectils_results_{sample}.csv"))
  }
)


## Plots ----

projectils_data <-
  dir_ls("data_wip", regex = "projectils") %>% 
  map_dfr(read_csv, .id = "filename") %>% 
  extract(filename, into = "sample", regex = "results_(.+)\\.csv") %>% 
  mutate(
    sample =
      rename_patients(sample) %>%
      rename_dong() %>%
      as_factor() %>% 
      fct_relevel(c(PATIENT_ORDER, files_dong$tumor_id))
  )

Embeddings(ref, "umap") %>%
  as_tibble() %>% 
  mutate(cell_type = Idents(ref)) %>% 
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(
    aes(color = cell_type),
    size = 0.1
  ) +
  # geom_point(
  #   data = projectils_data %>% filter(sample %in% c("M1", "T200 (M)")),
  #   size = 0.1
  # ) +
  geom_density_2d(
    data = projectils_data %>%
      # filter(sample %in% c("M1", "T200 (M)")) %>% 
      {.},
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
  facet_wrap(vars(sample), nrow = 4) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("comparison/projection", width = 420, height = 297)


projectils_data %>%
  mutate(
    predicted_type = factor(
      predicted_type,
      levels = levels(ref$functional.cluster)
    )
  ) %>%
  ggplot(aes(predicted_type)) +
  geom_bar(aes(fill = predicted_type), show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  xlab(NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = ADRENAL_MEDULLA_COLORS) +
  coord_flip() +
  facet_wrap(vars(sample), scales = "free_x", nrow = 4) +
  theme_bw() +
  theme(
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("comparison/predicted_cell_types", width = 350)


projectils_data %>%
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
ggsave_default("comparison/predicted_cell_types2", height = 210, width = 297)



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



# Integration -------------------------------------------------------------

## Analysis ----

cds <-
  readRDS("data_wip/tumor_data.rds") %>%
  {merge(.[[1]], .[-1])} %>% 
  {new_cell_data_set(.$RNA@counts, cell_metadata = .@meta.data)}

cds <- 
  preprocess_cds(cds, verbose = TRUE) %>% 
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
plot_cells(cds) + coord_fixed()

set.seed(42)
cds_aligned <-
  align_cds(cds, alignment_group = "sample", verbose = TRUE) %>% 
  reduce_dimension(
    reduction_method = "UMAP",
    preprocess_method = "Aligned",
    verbose = TRUE
  ) %>% 
  cluster_cells(k = 20, random_seed = 42, verbose = TRUE)
  
plot_cells(cds_aligned) + coord_fixed()  

list(
  colData(cds) %>%
    as_tibble(rownames = "cell") %>%
    select(cell, sample),
  reducedDim(cds, "UMAP") %>%
    magrittr::set_colnames(c("umap_1_unaligned", "umap_2_unaligned")) %>% 
    as_tibble(rownames = "cell"),
  reducedDim(cds_aligned, "UMAP") %>%
    magrittr::set_colnames(c("umap_1_aligned", "umap_2_aligned")) %>% 
    as_tibble(rownames = "cell"),
  clusters(cds_aligned) %>%
    enframe(name = "cell", value = "cluster")
) %>%
  reduce(left_join, by = "cell") %>% 
  saveRDS("data_wip/metadata_tumor_aligned.rds")


## Plots ----

tumor_metadata <-
  readRDS("data_wip/metadata_tumor_aligned.rds") %>% 
  mutate(
    sample =
      rename_patients(sample) %>%
      rename_dong() %>%
      as_factor() %>% 
      fct_relevel(c(PATIENT_ORDER, files_dong$tumor_id)),
    group = str_sub(sample, 1, 1),
    mycn_high = group == "M" | sample %in% c("T162", "T200", "T230")
  ) %>% 
  slice_sample(prop = 1)
tumor_metadata

ggplot(tumor_metadata, aes(umap_1_unaligned, umap_2_unaligned)) +
  geom_point(aes(color = sample), size = 0.1) +
  scale_color_hue(guide = guide_legend(override.aes = list(size = 3))) +
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("comparison/umap_unaligned")

ggplot(tumor_metadata, aes(umap_1_aligned, umap_2_aligned)) +
  geom_point(aes(color = sample), size = 0.1) +
  coord_fixed() +
  scale_color_hue(guide = guide_legend(override.aes = list(size = 4))) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("comparison/umap_aligned")

tumor_metadata %>%
  # filter(sample %in% c("M1", "T230 (M)")) %>% 
  ggplot(aes(umap_1_aligned, umap_2_aligned)) +
  geom_point(
    data = tumor_metadata %>% select(!sample),
    size = 0.1,
    color = "gray80"
  ) +
  geom_point(aes(color = sample), size = 0.1) +
  coord_fixed() +
  scale_color_hue(guide = guide_legend(override.aes = list(size = 4))) +
  facet_wrap(vars(sample), nrow = 4) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
ggsave_default("comparison/umap_samples", width = 420, height = 297)
