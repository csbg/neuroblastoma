# Integrate scRNA-seq datasets via alignment in batchelor and export metadata.
#
# The intergrated dataset is exported to rna_integrated_batchelor.rds.
#
# Metadata is exported to metadata_monocle.csv with columns
# * cell – barcode as used by Seurat
# * umap_1 ↓
# * umap_2 – UMAP coordinates
# * cluster_[k] – cluster IDs for different numbers k of nearest neighbors
# * partition_[k] – partition IDs for different k
#
# @DEPI rna_merged.rds
# @DEPO rna_integrated_batchelor.rds
# @DEPO metadata_batchelor.csv

library(Seurat)
library(scuttle)
library(scran)
library(scater)
library(batchelor)
library(bluster)
library(tidyverse)
library(fs)
source("common_functions.R")



# Parameters --------------------------------------------------------------

# the merged dataset
merged_dataset <- "data_generated/rna_merged.rds"

# folder where results are saved
out_dir <- "data_generated"



# Load data ---------------------------------------------------------------

nb <-
  readRDS(merged_dataset) %>% 
  SplitObject("sample")



# Normalization & feature selection ---------------------------------------

get_hvgs <- function(seurat_object) {
  batch <-
    rownames(seurat_object@meta.data)[1] %>%
    str_extract("\\d+")
  info("Normalizing batch {batch}")
  
  sce <- 
    seurat_object %>% 
    as.SingleCellExperiment() %>% 
    logNormCounts()
  metadata(sce) <- list(batch = batch)
  
  gene_var <- modelGeneVar(sce)
  top_genes <- getTopHVGs(gene_var, prop = 0.1)
  
  vis_data <- tibble(
    gene = rownames(gene_var),
    mean = metadata(gene_var)$mean,
    var = metadata(gene_var)$var,
    trend = metadata(gene_var)$trend(mean),
    highly_variable = gene %in% top_genes
  )
  
  rowSubset(sce) <- top_genes
  list(sce = sce, gene_var = gene_var, vis_data = vis_data)
}

nb <- map(nb, get_hvgs)

nb %>%
  map_dfr(~.$vis_data, .id = "sample") %>% 
  arrange(highly_variable) %>% 
  ggplot(aes(mean, var)) +
  geom_point(aes(color = highly_variable), size = 0.5, alpha = .25) +
  geom_line(aes(y = trend), color = "#e7298a") +
  facet_wrap(vars(sample)) +
  xlab("mean normalized log-expression") +
  ylab("variance") +
  scale_color_manual(
    name = "highly\nvariable",
    values = c("#666666", "#66a61e"),
    guide = guide_legend(override.aes = list(size = 2, alpha = 1))
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave_default("mean_variance_relationship")

nb_gene_var <- map(nb, ~.$gene_var)
nb <- map(nb, ~.$sce)


# Dimensionality reduction ------------------------------------------------

run_pca <- function(sce) {
  info("Running PCA on batch {metadata(sce)$batch}")
  runPCA(sce, subset_row = rowSubset(sce))
}

nb <- map(nb, run_pca)



# Integration -------------------------------------------------------------

nb <- multiBatchNorm(nb)
combined_top_genes <- combineVar(nb_gene_var)$bio > 0

integrate_without_correction <- function(data, selected_genes) {
  info("Merging data ...")
  data <-
    nb %>% 
    map(
      function(sce) {
        rowData(sce) <- NULL
        sce
      }
    ) %>% 
    reduce(cbind)
  
  info("Running PCA ...")
  data <- runPCA(data, subset_row = selected_genes)
  
  info("Running UMAP ...")
  data <- runUMAP(data, dimred = "PCA", min_dist = 0.1)
  
  info("... done!")
  data
}

nb_uncorrected <- integrate_without_correction(nb, combined_top_genes)



integrate_with_correction <- function(data, selected_genes) {
  merge_order <- 
    map_dfr(
      nb,
      ~tibble(sample = colData(.)$sample[1], group = colData(.)$group[1])
    ) %>% 
    mutate(i = row_number()) %>% 
    group_by(group) %>% 
    summarise(i = list(i)) %>% 
    arrange(desc(group)) %>% 
    pull(i)
  
  set.seed(42)
  
  info("Performing fast MNN correction ...")
  data <- fastMNN(data, subset.row = selected_genes, merge.order = merge_order)
  
  info("Running UMAP ...")
  data <- runUMAP(data, dimred = "corrected", min_dist = 0.1,
                  metric = "cosine", nn_method = "annoy")
  
  info("... done!")
  data
}

nb_integrated <- integrate_with_correction(nb, combined_top_genes)

# metadata(nb_integrated)$merge.info$lost.var
plotUMAP(nb_integrated, colour_by = "batch")



# Clustering --------------------------------------------------------------

k_values <- c(20, 50)

snn_clusters <-
  k_values %>% 
  map(
    ~clusterRows(
      reducedDim(nb_integrated, "corrected"),
      NNGraphParam(cluster.fun = "louvain", k = .)
    )
  ) %>%
  set_names(str_c("cluster_", k_values))

plotUMAP(nb_integrated, colour_by = I(snn_clusters[[2]]))



# Save data ---------------------------------------------------------------

saveRDS(nb_integrated, path_join(c(out_dir, "rna_integrated_batchelor.rds")))


nb_metadata <-
  bind_cols(
    reducedDim(nb_uncorrected, "UMAP") %>%
      magrittr::set_colnames(c("umap_1_unaligned", "umap_2_unaligned")) %>%
      as_tibble(rownames = "cell"),
    reducedDim(nb_integrated, "UMAP") %>%
      magrittr::set_colnames(c("umap_1", "umap_2")) %>% 
      as_tibble(),
    bind_cols(snn_clusters)
  ) %>% 
  rename_with(~str_c(., "_batchelor"), .cols = !cell)

nb_metadata %>% write_csv(path_join(c(out_dir, "metadata_batchelor.csv")))
