source("analysis_tumor/packages.R")
source("styling.R")



# Load data ---------------------------------------------------------------

DEG <- readRDS("analysis_tumor/data_generated/deg.rds")
GSEA <- readRDS("analysis_tumor/data_generated/gsea.rds")
DEG_trajectory <- readRDS("analysis_tumor/data_generated/deg_trajectory.rds")
cds <- readRDS("analysis_tumor/data_generated/cds.rds")



# Create tables -----------------------------------------------------------

cds@colData %>% 
  as_tibble(rownames = "cell") %>%
  transmute(
    Cell = cell,
    Patient = rename_patients(sample),
    # Patient ID = sample,
    Group = rename_groups(group),
    Cluster = cluster,
    Phase = phase
  ) %>%
  save_table("Sx_overview", "Overview")

DEG %>% 
  select(
    Gene = gene_short_name,
    Cluster = cluster,
    Estimate = estimate,
    Q_Value = q_value
  ) %>% 
  arrange(desc(Estimate)) %>% 
  save_table("Sx_cluster_de", "Cluster DE")


GSEA %>%
  select(!c(Old.P.value, Old.Adjusted.P.value)) %>%
  filter(
    !db %in% c(
      "Human_Gene_Atlas",
      "Chromosome_Location",
      "Jensen_TISSUES",
      "Jensen_DISEASES",
      "ChEA_2016",
      "ENCODE_TF_ChIP-seq_2015"
    )
  ) %>% 
  save_table("Sx_GSEA", "GSEA")

DEG_trajectory[["DE along trajectory"]] %>% 
  save_table("Sx_trajectory_de", "Trajectory DE")

DEG_trajectory[["Overlap_Trajectory_DEG"]] %>% 
  save_table("Sx_overlap_trajectory_deg", "Overlap trajectory")

DEG_trajectory[["Overlap_DEG_Trajectory"]] %>% 
  save_table("Sx_overlap_cluster_de", "Overlap cluster DE")
