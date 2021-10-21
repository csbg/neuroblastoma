source("nb_Packages.R")
source("nb_PlotFunctions.R")
source("styling.R")

paths <- read_csv("nb_Paths.csv", col_names = c("1","2","3","4","5","6"))
path <- paths$`1` #ubuntu data
path_plot <- paths$`2` #ubuntu plots

################################################################################
# Data Import:
################################################################################
# monocle Data import:

# Unaligned Dataset
unaligned_cds <- readRDS(
  paste0(path,"/nb_CDS_unaligned.rds")
  )

# rename metadata of CDS according to styling.R
colData(unaligned_cds) <- cbind(
  colData(unaligned_cds), 
  tibble(
  Group = rename_groups(colData(unaligned_cds)$group),
  Patient = rename_patients(colData(unaligned_cds)$sample)
  )
)
# Import of 4 Cluster CDS
cds <- 
  readRDS(
    paste0(path,"/nb_CDS.rds")
  ) 
  
colData(cds) <- cbind(
  colData(cds), 
  tibble(
    Group = rename_groups(colData(cds)$group),
    Patient = rename_patients(colData(cds)$sample)
  )
)

#NB subcluster "markers"

#subc_markers <- readRDS(
#  paste0(path,"/NB_cluster_markers.rds")
#) 



# GSEA Data Import
EnrData <- readRDS(
  file = paste0(path,"/nb_GSEA.rds")
  )

GSEA_data <- 
  EnrData %>% 
  as_tibble()

# Trajectory data import
trajectory_cds <- readRDS(
  file = paste0(path,"/nb_CDS_trajectory.rds")
  )

# Patient-Cluster Enrichment data import
Enr_Patients <- readRDS(
  file = paste0(path,"/nb_Patient_Enrichment_Cluster.rds")
  )

# Literature marker enrichment data import
Enr_Literature <- readRDS(
  file = paste0(path,"/nb_DEG_Enrichment_Jansky.rds")
)

# Differentially Expressed genes Cl4
DEG <- readRDS(
  file = paste0(path,"/nb_DEG.rds")
)

################################################################################
#  Plots:
################################################################################
################################################################################
# Upper Row:
# Settings
Side_quadratic_panel <- 3.6
Side_long_panel_x <- 6
Side_long_panel_y <- 7.2
Side_Heatmap_x <- 3
Side_Violin_x <- 3.8

################################################################################

# unaligned UMAP + group coloring, 1 patient highlighted
plot_UMAP_patient(
  unaligned_cds
  ) 

ggsave_publication("unaligned_Patients",
                   height = Side_quadratic_panel,
                   width = Side_quadratic_panel,
                   legends = F
                   )

ggsave_publication("guides_Patients",
                   height = Side_quadratic_panel,
                   width = Side_quadratic_panel*1.5,
                   legends = T
                   )

################################################################################
# 4x Aligned UMAP + cluster/Patients for Group II/IV coloring
plot_UMAP_group_patient(
  cds,
  "II"
)

ggsave_publication(
  "CL4_II_UMAP_Patient",
  height = Side_quadratic_panel,
  width = Side_quadratic_panel,
  legends = F
  )

plot_UMAP_group_patient(
  cds,
  "III"
)

ggsave_publication(
  "CL4_III_UMAP_Patient",
  height = Side_quadratic_panel,
  width = Side_quadratic_panel,
  legends = F
)


plot_UMAP_group_patient(
  cds,
  "IV"
)

ggsave_publication(
  "CL4_IV_UMAP_Patient",
  height = Side_quadratic_panel,
  width = Side_quadratic_panel,
  legends = F
)

################################################################################
# NB Cluster coloring 
plot_UMAP_cluster(
  cds
)

ggsave_publication(
  "CL4_cluster_UMAP",
  height = Side_quadratic_panel,
  width = Side_quadratic_panel,
  legends = F
  )


################################################################################
# Bottom Row
################################################################################

################################################################################
# UMAP aligned trajectory + colored by cytoTRACE

plot_trajectory(
  trajectory_cds
  )

ggsave_publication("NB_trajectory",
                   legends = F,
                   height = Side_quadratic_panel,
                   width = Side_quadratic_panel
                   )
# guides
ggsave_publication("Trajectory_Guides", 
                   legends = T,
                   height = Side_quadratic_panel,
                   width = Side_quadratic_panel
                   )


################################################################################
# Enrichment patients in clusters

plot_enrich_ClusterPatient(
  data = Enr_Patients %>% 
    mutate(Patient = 
             factor(x=Patient, 
                    levels = c("M1","M2","M3","M4",
                               "A1","A2",
                               "S1","S2","S3","S4","S5"
                               )
                    )
    ),
  circle_significant = F
)

ggsave_publication(
  filename = "nb_Patients_Clusters_Enrich_Cl4",
  height = Side_quadratic_panel, 
  width = 1.6,
  legends = F
)

ggsave_publication(
  filename = "Patients_Clusters_Enrich_guides",
  height = Side_quadratic_panel, 
  width = 2.7,
  legends = T
)

################################################################################
# Enrichment Literature Markers
plot_enrich_Literature(
  data = Enr_Literature,
  circle_significant = F
)

ggsave_publication(
  filename = "nb_Literature_Markers_Enrich_Cl4",
  height = Side_quadratic_panel, 
  width = Side_quadratic_panel*0.96,
  legends = F
)

ggsave_publication(
  filename = "Literature_Markers_Enrich_guides",
  height = Side_quadratic_panel, 
  width = Side_quadratic_panel*1.3,
  legends = T
)

################################################################################
# DEG genes Heatmap

p <- plot_DEG_heatmap(
  data_DEG = DEG,
  data_CDS = cds,
  legend = T
  )

ggsave_publication(
  plot = p,
  filename = "nb_Heatmap_DEG",
  height = Side_long_panel_y, 
  width = Side_Heatmap_x,
  legend = F
  )

ggsave_publication(
  plot = p,
  filename = "Heatmap_Guides",
  height = Side_long_panel_y, 
  width = Side_Heatmap_x,
  legend = T
)

################################################################################
# Violinplot DE
            
# Variant 2: ggplot-violin from monocle object 

top_markers <- list(
  Cl1 = tibble(gene=c("STAT3","JUN")),
  Cl2 = tibble(gene=c("E2F1","PCNA")),
  Cl3 = tibble(gene=c("SOX11","KIF5B"))
)

ggsave_publication(
  filename = "nb_violin_CL1",
  plot = 
    plot_violin2(
      cds,
      top_markers,
      1
    ),
  type ="pdf",
  height = Side_long_panel_y/3, 
  width = Side_Violin_x,
  dpi=1000,
  legend = F
)

ggsave_publication(
  filename = "nb_violin_CL2",
  plot = 
    plot_violin2(
      cds,
      top_markers,
      2
    ),
  type ="pdf",
  height = Side_long_panel_y/3, 
  width = Side_Violin_x,
  dpi=1000,
  legend = F
)

ggsave_publication(
  filename = "nb_violin_CL3",
  plot = 
    plot_violin2(
      cds,
      top_markers,
      3
    ),
  type ="pdf",
  height = Side_long_panel_y/3, 
  width = Side_Violin_x,
  dpi=1000,
  legend = F
)

################################################################################
# GSEA Dotplot

plot_gsea(
  data = EnrData[["Selected"]],
  min_OR = 1,
  circle_significant = F
)

ggsave_publication(
  "NB_GSEA", 
  width = Side_long_panel_x, 
  height = Side_long_panel_y,
  legends = F
  )
log10(1.5)
#guides
ggsave_publication(
  "GSEA_guides", 
  width = Side_long_panel_x*1.25, 
  height = Side_long_panel_y,
  legends = T
  )

################################################################################

# Supplementary Figure S2

################################################################################

# 3 Top + 1 both = width = 4

# unaligned UMAP colored by aligned cluster
plot_UMAP_alignedcluster(unaligned_cds)

ggsave_publication(
  "S2_UMAP_aligned_CL", 
  width = 4, 
  height = 4,
  legends = F
)

ggsave_publication(
  "S2_UMAP_aligned_CL_guides", 
  width = 4, 
  height = 4,
  legends = T
)

# CellCycle bargraph 

plot_CC_patient(cds, 
                labs = F
                )

ggsave_publication(
  "S2_CellCycle_Patients", 
  width = 4, 
  height = 4,
  legends = F
)

ggsave_publication(
  "S2_CellCycle_Patients_guides", 
  width = 4, 
  height = 8,
  legends = T
)


# CytoTRACE Score violin plot

plot_Violin_cyto(cds,
                 labs = T
                 )

ggsave_publication(
  "S2_Violin_CytoTRACE", 
  width = 4, 
  height = 4,
  legends = F
)

# Violinplot library_size

plot_Violin_cluster(cds)

ggsave_publication(
  "S2_Violin_library", 
  width = 4, 
  height = 4,
  legends = F
)

#Violinplot n_features

plot_Violin_cluster(cds)

ggsave_publication(
  "S2_Violin_features", 
  width = 4, 
  height = 4,
  legends = F
)

# Violinplot percent_mito

plot_Violin_cluster(cds)

ggsave_publication(
  "S2_Violin_mitopercent", 
  width = 4, 
  height = 4,
  legends = F
)

# Depracated
################################################################################
# Variant 1: seurat object wrapper ggplot expression function 
p1 <- list( Cl1 = plot_violin(
  cds,
  subc_markers[[1]] %>% 
    filter(
      gene == "BDNF" |
        gene == "STAT3" 
    )
),
Cl2 = plot_violin(
  cds,
  subc_markers[[2]] %>% 
    filter(
      gene == "CDK1" |
        gene == "E2F1" |
        gene == "UBE2T"|
        gene == "PCNA"
    )
),
Cl3 =  plot_violin(
  cds,
  subc_markers[[3]] %>% 
    filter(
      gene == "CDC20"|
        gene == "NEFL" |
        gene == "PRKCA"|
        gene == "SYT4"
    )
),
Cl4 = plot_violin(
  cds,
  subc_markers[[4]] %>% 
    filter(
      gene == "TP53" |
        gene == "MCM7" |
        gene == "MDM4" | 
        gene == "CAND1"
    )
)
)
