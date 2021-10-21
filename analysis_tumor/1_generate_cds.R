source("analysis_tumor/packages.R")


features <- read_tsv(
 "data_raw/COUNT/MF220_NB_BM_Patient1_transcriptome/filtered_feature_bc_matrix/features.tsv.gz",
  col_names = c("ENSEMBL", "gene_short_name", "X1")
)

cds <- readRDS("data_generated/rna_decontaminated.rds")

# add cell metadata
cds@colData <-
  read_rds("data_generated/metadata.rds") %>%
  mutate(Size_Factor = colData(cds)$Size_Factor) %>%     
  column_to_rownames("cell") %>% 
  as("DataFrame")

# use soupX counts as count matrix
assayNames(cds) <- c("raw_counts", "counts")  

# add additional gene metadata (ENSEMBL ID)
rowData(cds) <- features

# select NB tumor cells
cds <- cds[, colData(cds)$cellont_abbr == "NB"]

# export data
saveRDS(cds, "analysis_tumor/data_generated/cds_raw.rds")

