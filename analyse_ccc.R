# @DEPI rna_decontaminated.rds
# @DEPI metadata.rds
# @DEPO ccc_cellchat_object.rds
# @DEPO ccc_signaling_data.rds

library(SingleCellExperiment)
library(CellChat)
library(tidyverse)



# Load data ---------------------------------------------------------------

nb_metadata <- 
  read_rds("data_generated/metadata.rds") %>% 
  select(
    cell, library_size, n_features, sample,
    group, cellont_name, cellont_abbr, 
  ) %>% 
  filter(cellont_abbr != "other", group != "I") %>% 
  mutate(
    cellont_abbr =
      cellont_abbr %>% 
      fct_drop() %>% 
      fct_relevel("NB", "pDC", "M", "B", "T", "NK", "SC", "E")
  ) %>% 
  column_to_rownames("cell")

expr_data <- 
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  assay("soupx_counts") %>% 
  normalizeData(scale.factor = 1e6) %>% 
  magrittr::extract(, rownames(nb_metadata))



# Analyze data ------------------------------------------------------------

cellchat <- createCellChat(
  expr_data,
  meta = nb_metadata,
  group.by = "cellont_abbr",
)
cellchat@DB <- CellChatDB.human

# only keep genes in CellChatDB
cellchat <- subsetData(cellchat)

# find overexpressed ligands/receptors
cellchat <- identifyOverExpressedGenes(cellchat)

# find overexpressed interactions
cellchat <- identifyOverExpressedInteractions(cellchat)

# calculate communication probability between cell groups
cellchat <- computeCommunProb(cellchat, seed.use = 1)

# filter communications occurring with less than 10 cells in a cell group
cellchat <- filterCommunication(cellchat)

# calculate communication probability on the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# calculate aggregated communication network
cellchat <- aggregateNet(cellchat)

# determine network centrality pathways
cellchat <- netAnalysis_computeCentrality(cellchat)

# single LR dotplot data
df_net <-
  cellchat %>% 
  subsetCommunication(thresh = NA) %>% 
  as_tibble() %>% 
  group_by(interaction_name_2) %>% 
  mutate(prob.norm = prob / max(prob)) %>% 
  ungroup() %>% 
  mutate(
    source = factor(source, c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")),
    target = factor(target, c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")),
  )



# Save data ---------------------------------------------------------------

cellchat %>% saveRDS("data_generated/ccc_cellchat_object.rds")
df_net %>% saveRDS("data_generated/ccc_signaling_data.rds")
