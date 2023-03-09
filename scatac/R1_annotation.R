# Load required packages
pacman::p_load(Signac, Seurat, 
               scater, tidyverse, ggplot2, 
               patchwork, hrbrthemes, ggbeeswarm, 
               enrichR, UCell, RColorBrewer)
source('/media/AGFORTELNY/people/rohit/projects/neuroblastoma/styling_atac.R')

# Load data (Signac object)
nb  <- 
  readRDS(
    file = '/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized.rds')

nb <- NormalizeData(
  object = nb,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(nb$nCount_RNA))

DefaultAssay(nb) <- 'RNA'
load("/media/AGFORTELNY/people/rohit/projects/markers.RData")

# Add marker information to Signac object
nb <- AddModuleScore_UCell(nb, features = markers, 
                           assay = 'RNA', slot = "counts")
signature.markers <- paste0(names(markers), "_UCell")

# Plot signature Violin plot

VlnPlot(nb, features = signature.markers, log = TRUE, pt.size = 0, ncol = 1) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb() &
  theme(legend.position = 'none') 

colfunc <- colorRampPalette(c("#e66101", "#5e3c99"))
colfunc(10)

# Plot signature Feature plot

FeaturePlot(
  nb, reduction = "umap_scopen", pt.size = 0.25,
  features = signature.markers[1:7], ncol = 3, 
  cols = c("#e66101","grey90", "#2c7bb6"), order = T) &
  theme_nb() &
  theme(legend.key.size = unit(0.3, 'cm'), 
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'), 
        legend.title = element_text(size=14), 
        legend.text = element_text(size=10)) 

########################## CELL MARKERS ########################################
################################################################################

# https://panglaodb.se/
panglaodb <- read.table('~/projects/PanglaoDB_markers_27_Mar_2020.tsv', 
                        sep = '\t', header = T)
# DOI: 10.1016/j.cell.2019.04.040
baryawno <- read.table('~/projects/Baryawno_media-3.csv',
                       sep = ';', header = TRUE, comment.char = '#')
# http://biocc.hrbmu.edu.cn/CellMarker/
cell_markers <- read.table("~/projects/CellMarker.csv", 
                           sep = ",", header = T)

# B cell
panglaodb_bcell <- panglaodb %>% 
  filter(cell.type == "B cells") %>% 
  pull(official.gene.symbol)
baryawno_bcell <- baryawno %>% 
  pull(Mature.B)
bcell <- unique(c(panglaodb_bcell, baryawno_bcell))

# Erythroid 
panglaodb_erythroid <- panglaodb %>% 
  filter(cell.type == 'Erythroid-like and erythroid precursor cells') %>% 
  pull(official.gene.symbol)
baryawno_erythroid <- baryawno %>% 
  pull(Erythroid)
erythroid <- unique(c(panglaodb_erythroid, baryawno_erythroid ))

# T cell
baryawno_tcell <- baryawno %>% 
  pull(Treg)
tcell <- unique(c(baryawno_tcell, 
                  c("CD2",  "CD3D", "CD3E", "CD3G", "CD7" )))
# Monocytes
panglaodb_monocytes <- panglaodb %>% 
  filter(cell.type == 'Monocytes') %>% 
  pull(official.gene.symbol)
baryawno_monocytes <- baryawno %>% 
  pull(Monocytes)
monocytes <- unique(panglaodb_monocytes, baryawno_monocytes)

# 'NK cells'
panglaodb_nk <- panglaodb %>% 
  filter(cell.type == 'NK cells') %>% 
  pull(official.gene.symbol)
baryawno_nk <- baryawno %>% 
  pull(NK)
nk <- unique(c(panglaodb_nk, baryawno_nk))

# Plasmacytoid dendritic cells
panglaodb_nk <- panglaodb %>% 
  filter(cell.type == 'Plasmacytoid dendritic cells') %>% 
  pull(official.gene.symbol)
baryawno_pdc <- baryawno %>% pull(pDC)

# NK T cells
panglaodb_nkt <- panglaodb %>% 
  filter(cell.type == 'Natural killer T cells') %>% 
  pull(official.gene.symbol)
baryawno_nkt <- baryawno %>% 
  pull(NKT)
nkt <- unique(c(panglaodb_nkt, baryawno_nkt))

# Plasmacytoid dendritic cells (pDCs)
baryawno_pdc <- baryawno %>% 
  pull(pDC)

# pDCs
baryawno_pdcs <- baryawno %>% 
  pull(pDC)
markers[['pDC']] <- baryawno_pdc[baryawno_pdc != ""]

# Basophils
panglaodb_basophils <- panglaodb %>% 
  filter(cell.type == 'Basophils') %>% 
  pull(official.gene.symbol)

# Neuroblastoma
featuresWolf <- c(
  "PHOX2B", "PHOX2A", "NCAM1", "B4GALNT1", "MYCN", "GATA2", 
  "DBH", "CHGB", "CHGA", "SYP","TH")
featuresFormicola <- c(
  'MYCN', 'HOXC6', 'SOX4', 'FOXP1',  'GFRA3', 'PTPRH', 
  'ADCY1', 'PRKACB', 'AKR1C1')
nbmarker <- unique(c(featuresFormicola, featuresWolf))

# memory B cell
baryawno_memBcell <- baryawno %>% 
  pull(memBcell)
panglaodb_memBcell <- panglaodb %>% 
  dplyr::filter(cell.type == 'B cells memory') %>% 
  dplyr::pull(official.gene.symbol)
memBcell <- unique(c(baryawno_memBcell, panglaodb_memBcell ))

# stem cells
stem1 <- c('OLIG1', 'OLIG2', 'SOX2', 'SOX8', 'ASCL1', 'POU3F3', 'HES6', 
           'POU3F2', 'SOX21', 'HEY2', 'SOX5', 'RFX4', 'KLF15', 'CITED1', 
           'LHX2', 'VAX2', 'MYCL1', 'SALL2', 'SOX1')
stem2 <- c('DLK1', 'CD114', 'G-CSFR', 'BMI', 'CD44', 'CD133', 
           'CD117', 'CD24', 'FZD6', 'ALDH1', 'LGR5', 'GPR49', 
           'TLX', 'ABCG2', 'NESTIN', 'JARID1B', 'SPDYA', 'TRPM7', 'L1-CAM')

stem3 <- c('OLIG1', 'OLIG2', 'SOX2', 'SOX8', 'ASCL1', 'POU3F3', 'HES6', 
           'POU3F2', 'SOX21', 'HEY2', 'SOX5', 'RFX4', 'KLF15', 'CITED1', 
           'LHX2', 'VAX2', 'MYCL1', 'SALL2', 'SOX1', 'DLK1', 'CD114', 
           'G-CSFR', 'BMI', 'CD44', 'CD133', 'CD117', 'CD24', 'FZD6', 
           'ALDH1', 'LGR5', 'GPR49', 'TLX', 'ABCG2', 'NESTIN', 'JARID1B', 
           'SPDYA', 'TRPM7', 'L1-CAM', 'NOTCH1', 'PlGF2', 'TRKB', 'LNGFR', 
           'CD133', 'KIT', 'NOTCH1', 'GPRC5C', 'PIGF2', 'TRKB',  'LNGFR')
stem <- unique(c(stem1, stem2, stem3))

# Hematopoietic stem cells
panglaodb_hsc <- panglaodb %>% 
  dplyr::filter(cell.type == 'Hematopoietic stem cells') %>% 
  dplyr::pull(official.gene.symbol)

#
bmsc <- cell_markers %>% 
  filter(Cell.Type == 'Bone marrow stem cell') %>% 
  pull(Cell.Marker) %>% noquote() %>% 
  gsub(',', ' ', .)  %>% 
  stringr::str_split(., "\\s+") %>% unlist()

markers <- list()
markers$tcell   <- c('CD7', 'CD3G', 'CD3E', 'CD3D')
markers$nktcell <- c('IL2RB', 'KLRB1', 'KLRD1', 'GATA3')
markers$bcell   <- c('VPREB1', 'TNFRSF13C', 'TNFRSF13B', 'TLR9', 
                     'CD24', 'BLNK',  'CD79B', 'CD79A', 'CD72',  
                     'TLC1A', 'POU2AF1', 'PAX5', 'MS4A1', 'IL21R', 'IGLL5',
                     'IGLL1', 'FCRLA', 'FAM129C', 'EBF1', 'DOK3', 'BCL11A' ) 
markers$monocytes <- c('APOBEC3A', 'CCL3', 'CCR2', 'CD163', 'CD33', 'CD36', 'CLEC7A',
                       'CMKLR1', 'CSF1R','CSF3R', 'FCN1', 'IFI30', 'ICAM1', 'IL1B', 
                       'ITGAX', 'LYZ', 'MEFV', 'PILRA', 'S100A9', 'S100A8' )
markers$erythroids <- c('AHSP', 'CA1', 'ERMAP', 'GYPA', 'GYPB', 'HBA1', 'HBB', 'HBD', 
                        'HBM', 'HEMGN', 'KLF1', 'PKLR')
markers$nblast <-     c('ADCY1', 'CHGA', 'CHGB', 'DBH', 'GATA2', 'NCAM1', 'PHOX2A', 'PHOX2B', 
                        'SOX4')
markers$stem <-       c('SOX8', 'POU3F2', 'KLF15')
save(markers, file = "/media/AGFORTELNY/people/rohit/projects/markers.RData")
