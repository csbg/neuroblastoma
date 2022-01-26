# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds
# @DEPO dge_pb_results.rds

library(monocle3)
library(muscat)
library(scater)
library(tidyverse)
library(enrichR)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(viridis)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb <-
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  logNormCounts(assay.type = "soupx_counts")

nb_metadata <- readRDS("data_generated/metadata.rds")

# tumor infiltration rate
tif <-
  nb_metadata %>%
  group_by(sample) %>%
  summarise(tif = sum(cellont_abbr == "NB") / n())

nb@colData <-
  nb_metadata %>%
  mutate(
    sample =
      str_c(group, sample, sep = ".") %>% 
      as_factor() %>%
      fct_relevel(str_sort) %>% 
      fct_relabel(~str_extract(.x, "\\d+_\\d+")),
    Size_Factor = colData(nb)$Size_Factor
  ) %>% 
  left_join(tif, by = "sample") %>% 
  column_to_rownames("cell") %>% 
  as("DataFrame")
rowData(nb)[["gene_short_name"]] <- rownames(nb)

# clusters that contain more than 1% of total cells
used_clusters <- 
  nb_metadata %>% 
  count(cellont_cluster) %>% 
  mutate(n = n / sum(n)) %>% 
  filter(n > 0.01) %>% 
  pull(cellont_cluster)



# Analyse data ------------------------------------------------------------

# create pseudobulk data
nb <- prepSCE(                                           # on cluster level:
  nb[, colData(nb)$cellont_cluster %in% used_clusters],  # nb
  kid = "cellont_abbr",                                  # "cellont_cluster" 
  gid = "group",
  sid = "sample",
  drop = TRUE
)

pb <- aggregateData(
  nb,
  assay = "counts",
  fun = "sum",
  by = c("cluster_id", "sample_id")
)
pb

# assemble model matrix
exp_info <-
  metadata(pb)$experiment_info %>%
  left_join(tif, by = c(sample_id = "sample")) %>% 
  column_to_rownames("sample_id")


# either (a) default analysis ...
# model_mat <- model.matrix(~ group_id, data = exp_info)

# ...or (b) regress out TIF
model_mat <- model.matrix(~ group_id + tif, data = exp_info)

colnames(model_mat) <-
  colnames(model_mat) %>%
  str_replace("group_id(.+)", "\\1_vs_I")
model_mat


plot_model_matrix <- function(mat, filename = NULL) {
  p <- Heatmap(
    mat[levels(metadata(pb)$experiment_info$sample_id), ],
    name = "value",
    col = circlize::colorRamp2(
      breaks = seq(0, 1, length.out = 7),
      colors = inferno(7)
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
  )
  ggsave_default(filename, plot = p, width = 70, height = 100)
  p
}

plot_model_matrix(model_mat, filename = "dge/model_matrix")


# continue with common workflow: detect DE genes
dge <- map(
  2:ncol(model_mat),
  ~pbDS(
    pb,
    design = model_mat,
    coef = .x,
    min_cells = 1
  )
)

dge_results <-
  map_dfr(
    dge,
    ~resDS(nb, .x, frq = TRUE) %>%
      as_tibble()
  ) %>% 
  mutate(
    cluster_id =
      as_factor(cluster_id) %>%
      fct_expand(levels(colData(nb)$cluster_id)) %>%
      fct_relevel(levels(colData(nb)$cluster_id))
  ) %>% 
  rename(contrast = coef)



# Filter results ----------------------------------------------------------

#' Filter DGE results and add a column "direction".
#'
#' @param data A dataframe as returned by `muscat::resDS()`.
#' @param contrast_frq Named list of character vectors indicating which groups
#'   should be included for frequency filtering in each contrast.
#' @param max_p Maximum p value.
#' @param max_p_adj Maximum adjusted p value. Note: The FDR is calculated by
#'   muscat/edgeR via p.adjust(method = "BH").
#' @param min_abs_log_fc Minimum absolute log fold change.
#' @param min_freq Minimum gene expression frequency (fracions of cells that
#'   expresses a given gene.)
#' @param remove_ribosomal If true, remove ribosomal proteins.
#'
#' @return A data frame.
filter_dge_results <- function(data,
                               contrast_frq,
                               max_p = Inf,
                               max_p_adj = 0.05,
                               min_abs_log_fc = 1,
                               min_freq = 0.1,
                               remove_ribosomal = TRUE) {
  # lookup table of gene frequencies in groups
  gene_frq <-
    data %>% 
    select(gene, cluster_id, ends_with("frq")) %>%
    distinct() %>% 
    pivot_longer(
      ends_with("frq"),
      names_to = "frq_col",
      names_pattern = "(.+)\\.",
      values_to = "frq"
    )
  
  contrast_frq <- map_dfr(
    contrast_frq,
    ~gene_frq %>% 
      filter(frq_col %in% .x) %>% 
      group_by(gene, cluster_id) %>% 
      summarise(frq = max(frq)),
    .id = "contrast"
  )
  
  res <-
    data %>%
    left_join(
      contrast_frq,
      by = c("gene", "cluster_id", "contrast")
    ) %>% 
    filter(
      frq >= min_freq,
      p_val <= max_p,
      p_adj.loc <= max_p_adj,
      abs(logFC) >= min_abs_log_fc
    ) %>% 
    mutate(direction = if_else(logFC > 0, "up", "down")) %>% 
    select(!matches(".frq"))
  
  if (remove_ribosomal) {
    ribo_proteins <-
      msigdbr(species = "Homo sapiens", category = "C2") %>% 
      filter(gs_name == "KEGG_RIBOSOME") %>% 
      pull(human_gene_symbol)
    res <- 
      res %>%
      filter(!gene %in% ribo_proteins)
  }
  
  res
}

contrast_frq <- list(
  II_vs_I = c("I", "II"),
  III_vs_I = c("I", "III"),
  IV_vs_I = c("I", "IV"),
  tif = c("I", "II", "III", "IV")
)

# for visualization and EnrichR analysis, filter with default settings
dge_results_filtered <- filter_dge_results(dge_results, contrast_frq)

# for GSEA, only remove low-frequent and ribosomal genes
dge_results_filtered_gsea <- filter_dge_results(
  dge_results,
  contrast_frq,
  max_p_adj = Inf,
  min_abs_log_fc = 0
)



# Enrichment analysis -----------------------------------------------------

#' Perform enrichment analysis via enrichr for upregulated genes in a given
#' contrast and cluster.
#'
#' @param data A DGE dataset.
#' @param contrast The selected contrast.
#' @param cluster The selected cluster.
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#' @param direction Use genes that are up- or downregulated, respectively.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding a column "db".
enrich_genes <- function(data,
                         contrast,
                         cluster,
                         dbs,
                         direction = c("up", "down")) {
  info("Cluster {cluster}, contrast {contrast}, {direction}regulated genes")
  
  direction <- match.arg(direction)
  
  top_genes <- 
    data %>% 
    filter(
      cluster_id == {{cluster}},
      contrast == {{contrast}},
      direction == {{direction}}
    ) %>% 
    pull(gene)
  
  if (length(top_genes) > 0) {
    info("genes for enrichr: {str_c(top_genes, collapse = ', ')}")
    enrichr(top_genes, dbs) %>% 
      keep(~nrow(.) > 0) %>% 
      bind_rows(.id = "db") 
  }
}

#' Perform enrichment anlalysis for all contrasts and clusters in a dataset.
#'
#' @param data A DGE dataset.
#' @param dbs Character vector containing valid enrichr databases
#'   as returned by `enrichR::listEnrichrDbs()`.
#'
#' @return A dataframe, which combines the dataframes returned by
#'   `enrichR::enrichr()` (empty results are removed) by adding three columns
#'   "contrast", "cluster", and "db".
enrich_all_genes <- function(data, dbs) {
  queries <- 
    data %>%
    distinct(contrast, cluster = cluster_id) %>% 
    mutate(direction = list(c("up", "down"))) %>% 
    unnest_longer(direction)
  
  enrichr_data <-
    queries %>%
    pmap(enrich_genes, data = data, dbs = dbs)
  
  queries %>%
    mutate(data = enrichr_data) %>%
    rowwise() %>%
    filter(!is.null(data)) %>%
    filter(nrow(data) > 0) %>%
    unnest(data) %>%
    mutate(contrast = as_factor(contrast) %>% fct_relevel(str_sort)) %>% 
    separate(
      Overlap,
      into = c("overlap_size", "geneset_size"),
      convert = TRUE
    )
}

enrichr_dbs <- c(
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "WikiPathways_2019_Human",
  "MSigDB_Hallmark_2020",
  "TRRUST_Transcription_Factors_2019"
)

enrichr_results <- enrich_all_genes(dge_results_filtered, enrichr_dbs)



# GSEA --------------------------------------------------------------------

#' Download enrichr databases in a format that can be used by fgsea.
#'
#' @param dbs Databases to download.
#'
#' @return A named list with names deriving from values in `dbs`. Each element
#'   is a named list. Names correspond to enrichr terms, values are character
#'   vectors that comprise all genes associated with the respective term.
get_enrichr_genesets <- function(dbs) {
  dbs %>% 
    map(
      function(db) {
        info("Downloading {db}")
        url <- paste0(
          "https://maayanlab.cloud/Enrichr/geneSetLibrary",
          "?mode=text&libraryName=",
          db
        )
        read_lines(url)
      }
    ) %>% 
    set_names(dbs) %>% 
    map(
      function(db) {
        m <- str_match(db, "(.+?)\\t\\t(.+)")
        terms <- m[, 2]
        genes <- m[, 3] %>% str_split("\\t")
        genes %>% 
          map(stringi::stri_remove_empty) %>% 
          set_names(terms)
      }
    )
}


#' Perform gene set enrichment analysis.
#'
#' @param data DGE data as returned by `filter_dge_results()`.
#' @param gene_sets Gene set as returned by `get_enrichr_genesets()`.
#'
#' @return A dataframe, comprising columns "db", "contrast", "cluster", as well
#'   as all columns in the result of `fgseaMultilevel()`.
perform_gsea <- function(data, gene_sets) {
  data %>% 
    distinct(contrast, cluster_id) %>% 
    mutate(db = list(names(gene_sets))) %>% 
    unnest_longer(db) %>% 
    pmap_dfr(
      function(contrast, cluster_id, db) {
        ranked_genes <-
          data %>% 
          filter(contrast == {{contrast}}, cluster_id == {{cluster_id}}) %>%
          select(gene, logFC) %>%
          deframe()
        ranked_genes <- ranked_genes[!is.na(ranked_genes)]
        
        info("GSEA of contrast {contrast}, cluster {cluster_id}, ",
             "db {db} ({length(ranked_genes)} genes)")
        
        fgseaMultilevel(
          gene_sets[[db]],
          ranked_genes,
          eps = 0,
          nPermSimple = 10000
        ) %>%
          as_tibble() %>%
          mutate(
            db = {{db}},
            contrast = {{contrast}},
            cluster = {{cluster_id}},
            .before = 1
          )
      }
    )
}

# already defined in secion 'enrichment analysis',
# but we repeat it here for convenience
enrichr_dbs <- c(
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "WikiPathways_2019_Human",
  "MSigDB_Hallmark_2020",
  "TRRUST_Transcription_Factors_2019"
)

enrichr_genesets <- get_enrichr_genesets(enrichr_dbs)

# remove mouse genes from TTRUST database
enrichr_genesets$TRRUST_Transcription_Factors_2019 <-
  enrichr_genesets$TRRUST_Transcription_Factors_2019 %>% 
  magrittr::extract(imap_lgl(., ~str_detect(.y, "human"))) %>% 
  set_names(str_extract, "\\w+")

gsea_results <- perform_gsea(dge_results_filtered_gsea, enrichr_genesets)



# DGE in tumor cluster ----------------------------------------------------

# preprocessing as above
nb_tumor <-
  nb[
    ,
    colData(nb)$group_id %in% c("II", "IV")
    & colData(nb)$cluster_id == "NB"
  ] %>% 
  prepSCE()

pb_tumor <- aggregateData(nb_tumor, assay = "counts")

model_mat_tumor <- model.matrix(
  ~ group_id + tif,
  data =
    metadata(pb_tumor)$experiment_info %>%
    mutate(group_id = group_id %>% fct_drop() %>% fct_rev) %>% 
    left_join(tif, by = c(sample_id = "sample")) %>% 
    column_to_rownames("sample_id")
)

colnames(model_mat_tumor) <-
  colnames(model_mat_tumor) %>%
  str_replace("group_id(.+)", "\\1_vs_IV")
model_mat_tumor

# perform DGE
dge_results_tumor <-
  pbDS(
    pb_tumor,
    design = model_mat_tumor,
    coef = 2,  # only test II_vs_IV
    min_cells = 1
  ) %>% 
  {resDS(nb_tumor, ., frq = TRUE)} %>%
  as_tibble() %>% 
  rename(contrast = coef)

# filter results
dge_results_filtered_tumor <- filter_dge_results(
  dge_results_tumor,
  list(II_vs_IV = c("II", "IV"))
)



# Save results ------------------------------------------------------------

list(
  cds = nb,
  metadata = nb_metadata,
  used_clusters = used_clusters,
  results = dge_results,
  results_filtered = dge_results_filtered,
  results_filtered_gsea = dge_results_filtered_gsea,
  enrichr = enrichr_results,
  gsea = gsea_results,
  cds_tumor = nb_tumor,
  results_tumor = dge_results_tumor,
  results_tumor_filtered = dge_results_filtered_tumor
) %>% 
  saveRDS("data_generated/dge_pb_results.rds")
