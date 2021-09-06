library(Seurat)
library(monocle3)
library(scuttle)
library(umap)
library(ProjecTILs)
library(tidyverse)
library(fs)
source("common_functions.R")
source("styling.R")



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
# files_dong <- tribble(
#   ~tumor_id, ~file, ~bad_header, ~high_risk, ~do_work,
#   "T10",  "GSM4088774_T10_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T19",  "GSM4088775_T19_gene_cell_exprs_table.xls",  FALSE, FALSE, NA,
#   "T27",  "GSM4088776_T27_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T34",  "GSM4088777_T34_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T40",  "GSM4088778_T40_gene_cell_exprs_table.xls",  FALSE, FALSE, NA,
#   "T44",  "GSM4088779_T44_gene_cell_exprs_table.xls",  FALSE, FALSE, NA,
#   "T69",  "GSM4088780_T69_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T71",  "GSM4088781_T71_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T75",  "GSM4088782_T75_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T92",  "GSM4088783_T92_gene_cell_exprs_table.xls",  FALSE, TRUE, TRUE,
#   "T162", "GSM4654669_T162_gene_cell_exprs_table.xls", TRUE, TRUE, FALSE,
#   "T175", "GSM4654670_T175_gene_cell_exprs_table.xls", FALSE, FALSE, NA,
#   "T188", "GSM4654671_T188_gene_cell_exprs_table.xls", FALSE, FALSE, NA,
#   "T200", "GSM4654672_T200_gene_cell_exprs_table.xls", TRUE, TRUE, TRUE,
#   "T214", "GSM4654673_T214_gene_cell_exprs_table.xls", TRUE, TRUE, TRUE,
#   "T230", "GSM4654674_T230_gene_cell_exprs_table.xls", TRUE, TRUE, FALSE,
# )
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

# evil hack to prevent conversion of human to mouse gene names
convert.orthologs.id <- function(x, ...) x
environment(convert.orthologs.id) <- asNamespace("ProjecTILs")
assignInNamespace("convert.orthologs", convert.orthologs.id, ns = "ProjecTILs")

tumor_data <- readRDS("data_wip/tumor_data.rds")
ref <- load.reference.map("data_wip/reference_atlas_adr.rds")

# query_projected <- make.projection(
#   query = tumor_data,
#   ref = ref,
#   filter.cells = FALSE
# )
# 
# query_projected
# 
# projected_umap_data <- 
#   query_projected %>% 
#   map_dfr(
#     ~Embeddings(., "umap") %>% as_tibble(rownames = "cell"),
#     .id = "sample"
#   ) %>% 
#   mutate(sample = rename_patients(sample)) %>%
#   {.}
# 
# Embeddings(ref, "umap") %>%
#   as_tibble() %>% 
#   mutate(cell_type = Idents(ref)) %>% 
#   ggplot(aes(UMAP_1, UMAP_2)) +
#   geom_point(
#     aes(color = cell_type),
#     size = 0.1
#   ) +
#   geom_point(
#     data = projected_umap_data,
#     size = 0.1) +
#   scale_color_hue(guide = guide_legend(override.aes = list(size = 3))) +
#   coord_fixed() +
#   facet_wrap(vars(sample)) +
#   theme_bw() +
#   theme(panel.grid = element_blank())
# ggsave_default("comparison/projection")

projected_umap_data <- imap_dfr(
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
      Embeddings("umap") %>%
      as_tibble(rownames = "cell")
    res %>% write_csv(str_glue("data_wip/umap_{sample}.csv"))
    res
  },
  .id = "sample"
)

umap_data <- 
  dir_ls("data_wip", regex = "umap") %>% 
  map_dfr(read_csv, .id = "filename") %>% 
  extract(filename, into = "sample", regex = "umap_(.+)\\.csv") %>% 
  mutate(sample = rename_patients(sample))

Embeddings(ref, "umap") %>%
  as_tibble() %>% 
  mutate(cell_type = Idents(ref)) %>% 
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(
    aes(color = cell_type),
    size = 0.1
  ) +
  geom_point(
    data = umap_data,
    size = 0.1) +
  scale_color_hue(guide = guide_legend(override.aes = list(size = 3))) +
  coord_fixed() +
  facet_wrap(vars(sample)) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave_default("comparison/projection", width = 420, height = 297)

