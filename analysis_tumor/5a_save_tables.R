#Export of DATA to excel etc.
source("nb_Packages.R")
source("styling.R")

DEG <- readRDS(
  paste0(path,"/nb_DEG.rds")
  )


GSEA <- readRDS(
  paste0(path, "/nb_GSEA.rds")
  )

DEG_trajectory <- readRDS(
  paste0(path, "/nb_DEG_trajectory.rds")
)

cds <- 
  readRDS(
  paste0(path, "/nb_CDS.rds")
  )

Overview <- 
  cds@colData %>% 
  as_tibble() %>%
  select("sample","group","cluster", "phase") %>%
  transmute("ID" = rownames(cds@colData),
            "Patient" = rename_patients(sample),
            "Patient ID" = sample,
            "Group" = rename_groups(group),
            "Cluster" = cluster,
            "Phase" = phase)

WriteXLS(x =
          list(
            Overview,
            DEG %>% 
              rename("Estimate" = "estimate",
                     "Q_Value" = "q_value"
                     ),
            GSEA %>%
              select(!c(Old.P.value, Old.Adjusted.P.value)) %>%
              filter(db != "Human_Gene_Atlas",
                     db != "Chromosome_Location",
                     db != "Jensen_TISSUES",
                     db != "Jensen_DISEASES",
                     db != "ChEA_2016",
                     db != "ENCODE_TF_ChIP-seq_2015"
                     ),
            DEG_trajectory[[1]],
            DEG_trajectory[[2]],
            DEG_trajectory[[3]]
          ),
         SheetNames = c("Overview",
                        "1-Cluster DE",
                        "2-GSEA",
                        "3-Trajectory DE",
                        "4-Overlap trajectory",
                        "5-Overlap cluster DE"
                        ),
         ExcelFileName = paste0(path,"/nb_SupplementTables.xls"), 
         AdjWidth = TRUE, 
         BoldHeaderRow = TRUE, 
         FreezeRow = 1,
         AutoFilter = T)

