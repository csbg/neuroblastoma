library(SingleCellExperiment)
library(CellChat)
library(tidyverse)



# Load data ---------------------------------------------------------------

expr_data <- 
  readRDS("data_generated/rna_decontaminated.rds") %>% 
  assay("soupx_counts") %>% 
  normalizeData(scale.factor = 1e6)


meta_NB <-
  readRDS("analysis_tumor/data_generated/cds.rds") %>% 
  clusters() %>% 
  enframe("cell", "cluster")
  
  
meta <-
  read_rds("data_generated/metadata.rds") %>% 
  left_join(meta_NB, by = "cell") %>%
  mutate(
    ID =
      cluster %>%
      fct_relabel(~paste0("NB-", .)) %>%
      coalesce(cellont_abbr)
  ) %>% 
  select(
    library_size, n_features, sample,
    group, cellont_name, cellont_abbr,
    cell, ID
  ) %>%
  relocate(cell)

# subset count matrix & metadata
meta_overall <- 
  meta %>% 
  filter(cellont_abbr != "other", group != "I") %>% 
  mutate(
    ID = fct_drop(ID),
    cellont_abbr =
      cellont_abbr %>% 
      fct_drop() %>% 
      fct_relevel("NB", "pDC", "M", "B", "T", "NK", "SC", "E")
  ) %>% 
  column_to_rownames("cell")

expr_data_overall <- expr_data[, rownames(meta_overall)]

# check
table(meta_overall$cellont_abbr)



# CCC inference -----------------------------------------------------------

cellchat <- createCellChat(
  expr_data_overall,
  meta = meta_overall,
  group.by = "cellont_abbr",
)

# add CellChatDB to cellchat object
cellchat@DB <- CellChatDB.human

# subset to relevant genes (contained in DB) to save computational cost
cellchat <- subsetData(cellchat)

# pre-processing: identify overexpressed signaling genes associated with each ID
cellchat <- identifyOverExpressedGenes(cellchat)

# pre-processing: identify overexpressed LR-interactions
cellchat <- identifyOverExpressedInteractions(cellchat)

# save detected overexpressed LR interactions
overexpressedLR <- identifyOverExpressedInteractions(
  cellchat,
  return.object = FALSE
)

overexpressedLR %>% saveRDS("data_generated/ccc_overexpressed_LR.rds")


cellchat <- computeCommunProb(
  cellchat, 
  seed.use = 1
)

# filter communications occurring with less than 10 cells in a cell group
cellchat <- filterCommunication(cellchat)

# calculate CCC at signalling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# calculate aggregated CCC network
cellchat <- aggregateNet(cellchat)

# determine network centrality pathways
cellchat <- netAnalysis_computeCentrality(cellchat)

cellchat %>% saveRDS("data_generated/ccc_cellchat_object.rds")



# Signaling data ----------------------------------------------------------

# extract enriched single LR pairs in all pathways in microenvironment
pairLR <- extractEnrichedLR(
  cellchat, 
  signaling = cellchat@netP$pathways
)

pairLR %>% saveRDS("data_generated/ccc_enriched_LR_interactions.rds")
 

# single LR dotplot data
df_net <-
  subsetCommunication(cellchat, thresh = NA) %>% 
  as_tibble() %>% 
  group_by(interaction_name_2) %>% 
  mutate(prob.norm = prob / max(prob)) %>% 
  ungroup() %>% 
  mutate(
    source = factor(source, c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")),
    target = factor(target, c("NB", "pDC", "M", "B", "T", "NK", "SC", "E")),
  )

df_net %>% saveRDS("data_generated/ccc_signaling_data.rds")
