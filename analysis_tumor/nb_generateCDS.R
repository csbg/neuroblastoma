# Options:
source("nb_Packages.R")           # Load packages

# SAS change file link!
# import & assign gene-metadata to variable so that ENSMBL-ids are integrated 
#   (required for CellCycle determination)
# comes with every raw dataset 
# (e.g., COUNT/MF220_NB_BM_Patient1_transcriptome/filtered_feature_bc_matrix/features.tsv.gz; 
# but the file should be the same for all datasets)

features <- read_tsv(
 "/home/ubuntu/mnt/agfortelny/PROJECTS/Neuroblastoma/data_raw/geo_submission/S5_features.tsv.gz",
  col_names = c("ENSEMBL", "gene_short_name", "X1")
) #change link!

# SAS change file link!
#CellDataSet import after passed QC (from GFS)
cds <- read_rds(
  paste0(path_DATA, "/rna_decontaminated.rds")
) #change link!

# SAS change file link!
# add cell metadata to CDS from WES GFS
cds@colData <- read_rds(
  paste0(path_DATA,"/metadata.rds")
) %>%
  mutate(Size_Factor = colData(cds)$Size_Factor) %>%     
  column_to_rownames("cell") %>% 
  as("DataFrame") #change link!

#use soupX counts as count matrice
assayNames(cds) <- c("raw_counts", "counts")  

#add additional gene metadata (ENSEMBL ID)
rowData(cds) <- features

################################################################################
#subset to Neuron-cluster ->  select NB tumor cells
################################################################################
sub <- row.names(
  subset(pData(cds), cluster_20 == 5 | cluster_20 == 32)
)

cds <- cds[, sub]

#export data
saveRDS(cds,file = paste0(path,"/nb_CDS_raw.rds")) #change link
