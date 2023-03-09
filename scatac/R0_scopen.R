
# Load packages and set seed ----------------------------------------------

pacman::p_load(Signac, Seurat, future, remote,
               GenomeInfoDb, UCell,
               EnsDb.Hsapiens.v86, 
               ggplot2, patchwork, 
               hrbrthemes)
set.seed(1234)
plan("multiprocess", workers = 16)
plan()

# Pre-processing workflow -------------------------------------------------
source('/media/AGFORTELNY/people/rohit/projects/neuroblastoma/styling_atac.R')
pathTo10X <- 
  '/media/AGFORTELNY/PROJECTS/Neuroblastoma/data_raw/bsf/COUNT/AGGR_ALL_ATAC'
counts <- Read10X_h5(file.path(pathTo10X, 'filtered_peak_bc_matrix.h5'))

metadata <- read.csv(
  file = file.path(pathTo10X,'singlecell.csv'),
  header = TRUE,
  row.names = 1)

neuro_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = file.path(pathTo10X,'fragments.tsv.gz'),
  min.cells = 0)

nblast <- CreateSeuratObject(
  counts = neuro_assay,
  assay = 'peaks',
  project = 'NB',
  names.field = 2, 
  names.delim = "-",
  meta.data = metadata)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(nblast) <- annotations

# compute nucleosome signal score per cell
nblast <- NucleosomeSignal(object = nblast)

nblast$nucleosome_group <- 
  ifelse(nblast$nucleosome_signal > 3, 'NS > 3', 'NS < 3')

nucleosome_score <- 
    FragmentHistogram(
    object = nblast, 
    group.by = 'nucleosome_group', 
    region = 'chr3-1-10000000') &  
    #scale_fill_brewer(palette="Set1") &
    scale_fill_manual(values = c('black', 'grey')) &
    theme(legend.position = "none") &
    labs(title = 'Nucleosome Signal (NS) score') &
    theme_nb()

# compute TSS enrichment score per cell
nblast <- 
  TSSEnrichment(
    nblast, fast = FALSE)

nblast$high.tss <- 
  ifelse(
    nblast$TSS.enrichment > 1, 'High', 'Low')

tss_score <- 
  TSSPlot(
  nblast, 
  group.by = 'high.tss') & 
  #scale_color_brewer(palette="Set1") & 
  scale_color_manual(values = c('black', 'grey')) &
  theme_nb()

# add blacklist ratio and fraction of reads in peaks
nblast$pct_reads_in_peaks <- 
  nblast$peak_region_fragments / nblast$passed_filters * 100
nblast$blacklist_ratio <- 
  nblast$blacklist_region_fragments / nblast$peak_region_fragments

add_to_metadata <- function(nb_metadata){
  # add patient group information
  groups <- list(c('2', '3', '6', '15', '16'),
                 c('1', '4', '5', '11'),
                 c('7','12'),
                 c('8', '9', '10','13', '14'))
  # rename patients 
  nb_metadata <-
    nb_metadata %>% 
    mutate(
      new.ident = recode_factor(
        orig.ident,  
        `2` = '2014_0102', `3` = '2020_1288', `6` = '2018_4252',
        `15` = '2016_1853', `16` = '2016_2950', `1` = '16_4503',
        `4` = '2018_1404', `5` = '2019_5754', `11` = '2019_5022',
        `7` = '2016_3924', `12` = '2005_1702', `8` = '2018_6056',
        `9` = '2020_1667', `10` = '2006_2684', `13` = '2018_1625',
        `14` = '2019_2495'),
      sample = recode_factor(
        new.ident,  
        `2014_0102` = 'C1', `2016_1853`	= 'C2', `2016_2950`	= 'C3',
        `2018_4252`	= 'C4', `2020_1288`	= 'C5', `16_4503`	= 'M1',
        `2018_1404`	= 'M2', `2019_5022`	= 'M3', `2019_5754`	= 'M4',
        `2005_1702`	= 'A1', `2016_3924`	= 'A2', `2006_2684`	= 'S1',
        `2018_1625`	= 'S2', `2018_6056`	= 'S3',`2019_2495`	= 'S4',
        `2020_1667`	= 'S5'),
      group = dplyr::case_when(orig.ident %in% groups[[1]] ~"I",
                               orig.ident %in% groups[[2]] ~"II", 
                               orig.ident %in% groups[[3]] ~"III",
                               orig.ident %in% groups[[4]] ~"IV"))
  return(nb_metadata)
}

nblast@meta.data <- add_to_metadata(nblast@meta.data)

# nblast$group <- nblast@meta.data$group

features <- c('pct_reads_in_peaks','peak_region_fragments',
              'TSS.enrichment', 'nucleosome_signal',
              'blacklist_ratio')

titles <- c('Percent reads in peaks', 'Fragments in peak region',
            'TSS Enrichment', 'Nucleosome signal score', 
            'Encode blacklisted ratio')

qc_plots_pre_filter <-
  VlnPlot(
  object = nblast,
  features = features,
  group.by = 'new.ident', 
  split.by =  'group',
  pt.size = 0,
  ncol = 2 ) & 
  theme_nb() &
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90)) 

qc_plots_pre_filter
ggsave_publication(filename = "qc_plots_pre_filter", 
                   width = 8*2.54, height = 11*2.54)

nblast <- subset(
  x = nblast,
  subset = peak_region_fragments > 500 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 25 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 3 &
    TSS.enrichment > 1
)
dim(nblast)

sample_colors <- 
  c("black", "black", "black", "black", "black",
    "darkolivegreen", "darkolivegreen", "darkolivegreen", "darkolivegreen",
    "firebrick2", "firebrick2",
    "blue3", "blue3", "blue3", "blue3", "blue3")

plot_VlnPlot <- VlnPlot( object = nblast,
           features = features , pt.size = 0, 
           group.by = 'sample', 
           split.by =  'group',
           ncol = 4, combine = FALSE) 

for (i in 1:length(x = plot_VlnPlot)) {
  plot_VlnPlot[[i]] <- 
    plot_VlnPlot[[i]] + 
    theme_nb() + 
    labs(title = titles[i]) +
    theme(legend.title = element_blank(),
          #legend.position = "none",
          #axis.text.x = element_text(angle = 90)
          ) +
    coord_flip() +
    theme(plot.title = element_text(size=6, face="bold", family="sans"))
}

qc_plots <- 
  ( (plot_VlnPlot[[1]] + theme(legend.position = 'none')) | 
      (plot_VlnPlot[[2]] + theme(legend.position = 'none')) | 
      (plot_VlnPlot[[3]] + theme(legend.position = 'none')) | 
      (plot_VlnPlot[[4]]+ theme(legend.position = 'none')) |
      (plot_VlnPlot[[5]])) / 
  ((umap_atac + guides(colour = guide_legend(override.aes = list(size=1.5))) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=0.5))) | 
     (umap_rna + guides(colour = guide_legend(override.aes = list(size=1.5)))+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(colour = "black", size=0.5)))) +
  plot_annotation(tag_levels = 'a') +
  theme(plot.tag = element_text(face="bold", size = 10))

qc_plots

ggsave_publication(filename = "qc_plots_patchwork2", 
                   width = 8.5*2.54, height = 6*2.54)

# plot_spacer()
qc_plots_reviewers <-
  ((plot_VlnPlot[[1]] + theme(legend.position = 'none')) | 
     (plot_VlnPlot[[2]]+ theme(legend.position = 'none')) | 
     (plot_VlnPlot[[3]] + theme(legend.position = 'none')) | 
     (plot_VlnPlot[[4]])) / 
     ((plot_VlnPlot[[5]])| 
     (tss_score & theme(legend.position="none",
       plot.title = element_text(size=6, face="bold",family="sans"))) |
      (nucleosome_score & theme(
        legend.position="none",
        plot.title = element_text(size=6, face="bold",family="sans")))) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face="bold", size = 10))

qc_plots_reviewers
ggsave_publication(filename = "qc_plots_reviewers3", 
                   width = 8.5*2.54, height = 4.5*2.54)
# Run scopen with reticulate ----------------------------------------------

library(reticulate)
conda_list()
use_condaenv(condaenv = '/home/people/rohit/miniconda3/envs/mypython3', required = TRUE)
sc <- import("scopen")
matDR <- t(as.matrix(sc$Main$scopen_dr(nblast@assays$peaks@counts)))

colnames(matDR) <- paste0("scOpen_", 1:ncol(matDR))
rownames(matDR) <- colnames(nblast@assays$peaks@counts)
nblast@reductions[['scOpen']] <- 
  CreateDimReducObject(embeddings = matDR, 
                       assay = DefaultAssay(nblast))
DepthCor(nblast, reduction = "scOpen", n = 20)


nblast <- 
  RunUMAP(object = nblast, 
          reduction = 'scOpen', 
          dims = 1:30, 
          reduction.name = "umap_scopen")
nblast <- 
  FindNeighbors(object = nblast, 
                reduction = 'scOpen', 
                dims = 1:30)
nblast <- 
  FindClusters(object = nblast, 
               nblastverbose = FALSE, 
               algorithm = 3)
DimPlot(object = nblast, 
        label = TRUE, 
        reduction = "umap_scopen") &
  hrbrthemes::theme_ipsum()

gene.activities <- GeneActivity(nblast)

# add the gene activity matrix to the Seurat object as a new assay and normalize it

nblast[['RNA']] <- 
  CreateAssayObject(counts = gene.activities)

nblast <- NormalizeData(
  object = nblast,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(nblast$nCount_RNA)
)

nbmarkers <- c('MYCN', 'DDX1','ARID1A' , 'HOXD10', 'EVX2', 'NBAS')

DefaultAssay(nblast) <- 'RNA'
FeaturePlot(
  object = nblast,
  features = nbmarkers,
  pt.size = 0.3,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey80', 'red'),
  reduction = "umap_scopen"
)

# saveRDS(nblast, 
#   file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized.rds")
nblast <- 
  readRDS(
    file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized.rds")

# saveRDS(nblast,
#   file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized_motifs_added.rds")

nblast <- 
  readRDS(
    file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized_motifs_added.rds")


# Add markers info with Ucell ---------------------------------------------

# annotate clusters
library(UCell)
load("/media/AGFORTELNY/people/rohit/projects/markers.RData")


nblast <- AddModuleScore_UCell(nblast, features = baryawno_list, 
                               assay = 'RNA', slot = "counts")
baryawno_list.markers <- paste0(names(baryawno_list), "_UCell")
VlnPlot(nblast, features = baryawno_list.markers[11], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") #&
#scale_y_continuous(limits = c(1,1.015)) &
#theme_nb()
pdf(file = 'plots/violin/ucell_baryawno.pdf',   # The directory you want to save the file in
    width = 8, 
    height = 4)
for (i in 1:length(baryawno_list.markers)){
  print(VlnPlot(nblast, features = baryawno_list.markers[i], 
                log = TRUE, pt.size = 0) &
          stat_summary(fun = median, geom='point', 
                       size = 2, colour = "black"))
  
}
dev.off()

nblast <- AddModuleScore_UCell(nblast, features = markers, 
                               assay = 'RNA', slot = "counts")
signature.markers <- paste0(names(markers), "_UCell")

# Signature markers

(t_ucell <-
    VlnPlot(nblast, features = signature.markers[1], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    #scale_y_continuous(limits = c(1,1.015)) &
    theme_nb())
(nkt_ucell <-
  VlnPlot(nblast, features = signature.markers[2], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb())
(b_ucell <- VlnPlot(nblast, features = signature.markers[3], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb())
(monocytes_ucell <- 
  VlnPlot(nblast, features = signature.markers[4], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb())
(erythroid_ucell <-
  VlnPlot(nblast, features = signature.markers[5], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb())
(nblast_ucell <-
  VlnPlot(nblast, features = signature.markers[6], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb())
(stem_ucell<-
  VlnPlot(nblast, features = signature.markers[7], log = TRUE, pt.size = 0) &
  stat_summary(fun = median, geom='point', size = 2, colour = "black") &
  theme_nb())
(nk_ucell<-
    VlnPlot(nblast, features = signature.markers[9], log = TRUE, pt.size = 0) &
    stat_summary(fun = median, geom='point', size = 2, colour = "black") &
    theme_nb())

(nblast_ucell/ monocytes_ucell/ b_ucell/ erythroid_ucell/ nkt_ucell/ t_ucell ) +
  plot_layout(guides = "collect")
ggsave_publication('ucell_marker_violin', 
                   width = 19, height= 28)


# Tcell_signature = c("CD2","CD3E","CD3D")
# Myeloid_signature = c("SPI1","FCER1G","CSF1R")

Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[1:6], 
  ncol = 3, 
  cols = c("grey90", "magenta"),
  order = T) &
  scale_colour_gradientn(
    colours = rev(brewer.pal(n = 5, 
                             name = "RdYlGn"))) 

Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[1:9], 
  ncol = 3, 
  #cols = rev(Blue2Gray8Steps),
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
  scale_fill_manual(values=ModifiedSpectralScheme11Steps) +
(sm1 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[1],
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) & 
    labs(title = "T cell markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(colours = rev(brewer.pal(
    #   n = 5,
    #   name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5)) )
  
(sm2 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[2], 
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "NKT cell markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) & 
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm3 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[3],
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "B cell markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm4 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[4],
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "Myeloid markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm5 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[5], 
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "Erythroid markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm6 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[6],
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "NB markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 

(sm8 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[8],
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "NK markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 
(sm9 <- Seurat::FeaturePlot(
  nblast,     
  reduction = "umap_scopen", 
  features = signature.markers[9],
  cols = c("grey90", "#ef8a62", "#2166ac"),
  order = T) &
    labs(title = "pDC markers") &
    theme_nb() &
    #theme(legend.position="none") &
    # scale_colour_gradientn(
    #   colours = rev(brewer.pal(n = 5, 
    #                            name = "RdYlGn"))) &
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=5), #change legend title font size
          legend.text = element_text(size=5))) 

(sm1 | sm8   | sm3) /(sm4 |sm5 |sm6 ) +
  plot_annotation(tag_levels = 'A')
ggsave_publication('umap_signatures2', 
                   width = 18, height= 9)

ggsave_publication('umap_signatures', 
                   width = 18, height= 9)
# "#04598C" "#FEFEB2"
markers$tcell
tcelldb <- dotplotdb(nblast, peak_annotation, tcell, "T cell")
bcelldb <- dotplotdb(nblast, peak_annotation, bcell, "B cell")

ggplot(data = bcelldb) +
  aes(x= Cluster, y = GENE) +
  geom_point(aes(size = Pct,  color= Avg)) +
  theme_ipsum() + ylab('') + xlab('Clusters') 
+
  lemon::facet_rep_grid(Name ~ ., scales = "free", space="free", 
                        repeat.tick.labels = 'bottom')
anno_path <- 
  '/media/AGFORTELNY/people/rohit/projects/neuroblastoma/idents3.csv'

# Change Seurat identities
annotate <- read.table(anno_path,
                       sep = ',',
                       header = FALSE,
                       strip.white = TRUE )
annot <- as.character(annotate$V2)
names(annot) <-annotate$V1
nblast <- Seurat::RenameIdents(nblast, annot)
nblast$activ.ident <- Idents(nblast)
nblast$seurat_clusters

SeuratObject <- readRDS(
  file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized.rds")

SeuratObject <- Seurat::RenameIdents(SeuratObject, 
                                     annot)

umap <- 
  DimPlot(
  object = nblast, 
  label = TRUE, 
  label.size = 3, 
  pt.size = 0.2,
  #cols = CELL_TYPE_COLORS,
  #group.by = 'new.ident',
  reduction = "umap_scopen") &
  #scale_colour_brewer(palette = "Set1") &
  #hrbrthemes::theme_ipsum() &
  theme_nb() &
  scale_color_manual(
    name = "cell type",
    values = CELL_TYPE_COLORS,
    labels =
      str_glue(
        "{names(CELL_TYPE_ABBREVIATIONS)} ({CELL_TYPE_ABBREVIATIONS})"
      )) 
  #theme(legend.position = "none")
  #guides(color = guide_legend(override.aes = list(size=10), ncol=1) ) &
  #theme(legend.text=element_text(size=10))
  
ggsave_publication('umap_wes_theme_nb_annot', 
                     width = 9, height= 6)

DimPlot(
  object = nblast, 
  label = TRUE, 
  group.by = 'seurat_clusters',
  reduction = "umap_scopen") &
  theme_nb() &
  labs(title = "Seurat clusters") &
  umap &
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank())
ggsave_publication('seurat_umap', 
                   width = 8, height= 7)
  
  # scale_colour_brewer(palette = "Set1") &
  # hrbrthemes::theme_ipsum()

cowplot::plot_grid(umap,  cov_plot1 / cov_plot2 ) 
(umap |cov_plot1 | cov_plot2 / (p3_tnfa /p3_e2f))

umap + plot_spacer() + 
  (cov_plot1 / cov_plot2) + 
  plot_spacer() + (p3_tnfa /p3_e2f) + plot_spacer()

ggpubr::ggarrange(umap, cov_plot1, theme_void(), theme_void(),
                  theme_void(), cov_plot2, theme_void(),theme_void(),
                  p3_tnfa , theme_void(), theme_void(),theme_void(),
                  p3_e2f,
          #labels = c("A", "B", "C"),
          ncol = 4, nrow = 4)
multiplot(umap,
          cov_plot1,cov_plot2,
          p3_tnfa,p3_e2f,ncol=4) 
multiplot(all_promoters, all_others, ncol=2)
ggpubr::ggarrange(all_promoters, all_others)

all_promoters
ggsave_publication('all_promoters_motif', 
                   width = 10, height= 15)



# Run Footprinting Analysis -----------------------------------------------

pacman::p_load(JASPAR2020, TFBSTools)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

DefaultAssay(nblast) <- 'peaks'
# add motif information
nblast <- AddMotifs(
  object = nblast,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# List of added motifs
motif.list <- 
  nblast@assays$peaks@motifs@motif.names
sorted.motifs <-  
  motif.list[order(unlist(motif.list),
                   decreasing=FALSE)]

# gather the footprinting information for sets of motifs
nblast <- Footprint(
  object = nblast,
  motif.name = c("GATA2", "CEBPA", "EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)
moi <- c('NFE2', 'SPIB', "PHOX2B", 
         "PHOX2A", 'IRF8',  'IRF1', 
         'IRF3', 'ETS1', 'ELF5', 
         'REST', "MYCN", "JUNB")

moi <- c('SPIB',"PHOX2A")
moi_not_found <- c('Bach1::Mafk', 'Runx1', 'NFE2L2')

nblast <- Footprint(
  object = nblast,
  motif.name = moi,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
nblast <- Footprint(
  object = nblast,
  motif.name = "JUNB",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

fp1 <- PlotFootprint(nblast, features = 'JUNB')
fp2 <- PlotFootprint(nblast, features = 'NFE2')
fp3 <- PlotFootprint(nblast, features = 'ELF5')
fp4 <- PlotFootprint(nblast, features = 'SPIB')
fp5 <- PlotFootprint(nblast, features = 'PHOX2B')
fp6 <- PlotFootprint(nblast, features = 'PHOX2A')
(fp1 | fp2 |fp3 ) /(fp4 |fp5 |fp6) &
  plot_layout(guides = 'collect') &
  theme_nb()


prow <- cowplot::plot_grid(
  fp1 & theme(legend.position="none"),
  fp2 & theme(legend.position="none"),
  fp3 ,
  fp4 & theme(legend.position="none"),
  fp5 & theme(legend.position="none"),
  fp6 ,
  align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F"),
  hjust = -1,
  nrow = 2,
  ncol = 3
)
prow

fp_1 <- PlotFootprint(nblast, features = 'JUNB') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells") 
fp_2 <- PlotFootprint(nblast, features = 'NFE2') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_3 <- PlotFootprint(nblast, features = 'ELF5') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_4 <- PlotFootprint(nblast, features = 'SPIB') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_5 <- PlotFootprint(nblast, features = 'PHOX2B') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")
fp_6 <- PlotFootprint(nblast, features = 'PHOX2A') & scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells")

prow2 <- cowplot::plot_grid(
  fp_1 & theme(legend.position="none") ,
  fp_2 & theme(legend.position="none") ,
  fp_3 ,
  fp_4 & theme(legend.position="none") ,
  fp_5 & theme(legend.position="none") ,
  fp_6 ,
  align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F"),
  hjust = -1,
  nrow = 2,
  ncol = 3
)
prow2
ggsave_publication('cell_type_footprinting', 
                   width = 30, height= 15)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  fp1 & theme(legend.box.margin = margin(0, 0, 0, 12))
)

cowplot::plot_grid(prow, legend, rel_widths = c(3, .4))

fp1 & 
  guides(fill="none") &
  theme(legend.position="none")
PlotFootprint(nblast, 
features = c('JUNB', 'NFE2', 'ELF5', 'SPIB', 'PHOX2B','PHOX2A')) &
  scale_color_manual(values = CELL_TYPE_COLORS) &
  labs(color = "Cells") &
  plot_layout(guides = 'collect') &
  theme_nb()


# plot the footprint data for each group of cells
p2 <- PlotFootprint(nblast, features = c("GATA2", "CEBPA", "EBF1"))
p2 + patchwork::plot_layout(ncol = 1)

p1 <- PlotFootprint(nblast, features =  moi[1:3])
p1 + patchwork::plot_layout(ncol = 3)

p3 <- PlotFootprint(nblast, features =  moi[1:3])
p3 + patchwork::plot_layout(ncol = 1)

p4 <- PlotFootprint(nblast, features =  moi[8:10])
p4 + patchwork::plot_layout(ncol = 1)

p5 <- PlotFootprint(nblast, features =  moi[11])
p5 + patchwork::plot_layout(ncol = 1)

names(table(nblast@active.ident))

# get ids of different cell types

create_bigwig_cell_ids <- function(object, cell_type) {
  sub_obj <- subset(x = nblast, subset = active.ident == 'NB')
  split_obj <- SplitObject(sub_obj, split.by = 'orig.ident' ) # "new.ident"
  split_obj_names <- names(split_obj)
  for (i in 1:length(split_obj)) {
    name <- split_obj_names[[i]]
    path <- 
      paste0('/media/AGFORTELNY/people/rohit/projects/bigwig/', 
                   tolower(cell_type))
    dir.create(file.path(path))
    filename <- file.path(path, paste0(name, '.csv'))
    cell_ids <- rownames(split_obj[[i]]@meta.data)
    write.csv(cell_ids, file=filename, row.names=F)
  }
}
my_list <- setdiff(names(table(nblast@active.ident)), "Memory_B_cell")
for (cell in my_list) {
  create_bigwig_cell_ids(nblast, cell)
}

create_bigwig_cell_ids(nblast, "NB")
# sub_obj <- subset(x = sub_obj, subset = orig.ident != '2' )
# sub_obj <- subset(x = sub_obj, subset = orig.ident != '3' )
# sub_obj <- subset(x = sub_obj, subset = orig.ident != '15' )
# sub_obj <- subset(x = sub_obj, subset = orig.ident != '16' )

sub_obj <- subset(x = nblast, subset = active.ident == 'NB')
sub_obj <- subset(x = sub_obj, subset = orig.ident != '16' )
split_obj <- SplitObject(sub_obj, split.by ='new.ident') # "new.ident" 'orig.ident' 
split_obj_names <- names(split_obj)

# Motif analysis ----------------------------------------------------------
EnrichedMotifs <- function(SeuratObject){
  enr_mot <- list()
  for (cell_type in levels(Idents(SeuratObject))) {
    message(cell_type)
    da_peaks <- FindMarkers(
      object = SeuratObject,
      ident.1 = cell_type,
      only.pos = TRUE,
      test.use = 'LR',
      min.pct = 0.05,
      latent.vars = 'nCount_peaks'
    )
    # get top differentially accessible peaks
    top.da.peak <- 
      da_peaks[da_peaks$p_val < 0.005, ] %>%
      mutate(odds_ratio = pct.1/pct.2) %>% 
      filter(quantile(odds_ratio, 0.75) < odds_ratio) %>% 
      row.names(.)
    # test enrichment
    enriched.motifs <- 
      FindMotifs(
      object = nblast,
      features = top.da.peak
    )
    # filter enriched motifs
    enriched.motifs <-
      enriched.motifs %>%
      mutate(cell_type = cell_type) # %>%
    # filter(qvalue < 0.05 & abs(log_odds_ratio) >0.75) 
    enr_mot[[cell_type]] <- enriched.motifs
  }
  return(enr_mot)
}

mot_enr_list <- EnrichedMotifs(nblast)

test_mot_m <- mot_enr_list$M %>% 
  mutate(log_odds_ratio = log(percent.observed/percent.background),
         qvalue = qvalue::qvalue(mot_enr_list$M$pvalue)$qvalues) %>% 
  filter(pvalue < 0.05 &
           qvalue < 0.05 & 
           abs(fold.enrichment) >2) %>% 
  arrange(qvalue) %>% head(100)


mot_enr_list_new <- list()
for (i in names(mot_enr_list)){
  print(i)
  my_df <- mot_enr_list[[i]]
  temp_df <-
    my_df %>% 
    mutate(
      log_odds_ratio = log(percent.observed/percent.background),
      qvalue = qvalue::qvalue(my_df$pvalue)$qvalues) %>% 
    filter(
      pvalue < 0.05 & qvalue < 0.05 ) %>% 
    arrange(qvalue) %>% head(15)
  mot_enr_list_new[[i]] <- temp_df
}


library(circlize)
col_fun = colorRamp2(c(0, 2, 4), c("khaki", "white", "red"))
col_fun(seq(0, 4))

do.call(rbind, mot_enr_list_new) %>% 
  mutate(cell_type = plyr::revalue(cell_type, c("NKT" = "NK"))) %>% 
  #ggplot(.,aes(x = cell_type, y = motif.name)) + 
  #geom_point(aes(size= fold.enrichment, 
  #               color = log_odds_ratio))
  dplyr::select(motif.name, cell_type, fold.enrichment) %>% 
  pivot_wider(names_from = "cell_type",
              values_from = "fold.enrichment",
              values_fill = 0) %>% 
  column_to_rownames(var = "motif.name") %>% 
  as.matrix() %>% t() %>% 
  ComplexHeatmap::Heatmap(col = col_fun, name = "fold change") 
  #pheatmap::pheatmap(.)

# top.da.peak_nb <- rownames(da_peaks_nb[da_peaks_nb$p_val < 0.005, ])
# Save NB cells -----------------------------------------------------------

nb_cluster <- subset(x = nblast, subset = active.ident == "NB")
saveRDS(nblast,
   file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_nb.rds")
nb_cluster <- readRDS(
  file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_nb.rds")

# # add patient group information
# groups <- list(c('2', '3', '6', '15', '16'),
#                c('1', '4', '5', '11'),
#                c('7','12'),
#                c('8', '9', '10','13', '14'))
# 
# nb_cluster@meta.data$group <- nb_cluster@meta.data %>% 
#   dplyr::mutate(
#     group = dplyr::case_when(orig.ident %in% groups[[1]] ~"C",
#                              orig.ident %in% groups[[2]] ~"M", 
#                              orig.ident %in% groups[[3]] ~"A",
#                              orig.ident %in% groups[[4]] ~"S")) %>% 
#   pull(group)
# 
# nb_cluster$group <- nb_cluster@meta.data$group

nb_cluster <- RunTFIDF(nb_cluster, method = 4)
nb_cluster <- FindTopFeatures(nb_cluster, min.cutoff = 'q5')
nb_cluster <- RunSVD(object = nb_cluster)

DepthCor(nb_cluster)

nb_cluster <- RunUMAP(
  object = nb_cluster,
  reduction = 'lsi',
  umap.method = 'uwot', 
  dims = 2:30)

nb_cluster <- FindNeighbors(
  object = nb_cluster,
  reduction = 'lsi',
  nn.method = 'rann',
  dims = 2:30)

nb_cluster <- FindClusters(
  object = nb_cluster,
  algorithm = 3,
  resolution = 1,
  verbose = FALSE)

DimPlot(
  object = nb_cluster, 
  #label = TRUE, 
  group.by = 'orig.ident',
  reduction = "umap") &
  hrbrthemes::theme_ipsum()


# Coverage plots ----------------------------------------------------------
pacman::p_load(Signac, Seurat, tidyverse, patchwork)

nblast <- 
  readRDS(
    file = "/media/AGFORTELNY/people/rohit/projects/nblast_scopen_gene_activity_normalized_motifs_added.rds")

source('/media/AGFORTELNY/people/rohit/projects/neuroblastoma/styling_atac.R')

# Coverage plots ----------------------------------------------------------
patients_col <- c(
  "I-2014_0102" = "navy" ,  "I-2018_4252" = "navy",  "I-2020_1288" = "navy",  
  "I-MF244_GNB" = "navy",   "I-MF244_GNM" = "navy", "II-16_4503"   = "royalblue",
  "II-2018_1404" = "royalblue", "II-2019_5022" = "royalblue", "II-2019_5754" = "royalblue", 
  "III-2005_1702" = "red",  "III-2016_3924" = "red", "IV-2006_2684"  = "seagreen",
  "IV-2018_1625" = "seagreen", "IV-2018_6056" = "seagreen", "IV-2019_2495" = "seagreen", 
  "IV-2020_1667" = "seagreen")

# add patient group information
groups <- list(c('2', '3', '6', '15', '16'),
               c('1', '4', '5', '11'),
               c('7','12'),
               c('8', '9', '10','13', '14'))

nblast@meta.data$group <- nblast@meta.data %>% 
  dplyr::mutate(
    group = dplyr::case_when(orig.ident %in% groups[[1]] ~"C",
                             orig.ident %in% groups[[2]] ~"M", 
                             orig.ident %in% groups[[3]] ~"A",
                             orig.ident %in% groups[[4]] ~"S")) %>% 
  pull(group)
nblast$group

group_col <- c('C' = 'navy',
               'M' = 'royalblue',
               'A' = 'red',
               'S' = 'seagreen')

myeloid <- subset(x = nblast, subset = activ.ident == "Monocytes")

peak <-'chr4-102521897-102522522'
ranges.show <- StringToGRanges(peak)
ranges.show$color <- "coral"
  
cov_plot_group <- CoveragePlot(
  object = myeloid,
  region.highlight = ranges.show,
  region = peak,
  #features = "NFKB1",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=16))

cov_plot_group <- CoveragePlot(
  object = nblast,
  region.highlight = ranges.show,
  region = peak,
  #features = "NFKB1",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=16))

peak1 <-'chr4-102521897-102522522'
ranges.show1 <- Signac::StringToGRanges(peak1)
ranges.show1$color <- "coral"
  
cov_plot1 <- Signac::CoveragePlot(
  object = myeloid,
  region.highlight = ranges.show1,
  region = peak1,
  #features = "NFKB1",
  annotation = TRUE,
  extend.upstream = 1000,
  extend.downstream  = 1000,
  group.by = 'group', 
  peaks = FALSE
) &
  scale_fill_manual( values = group_col ) &
  theme(text = element_text(size=10),
        legend.position = 'none') &
  guides(fill="none") &
  theme_nb()
cov_plot1

cov_plot1

peak2 <-'chr5-132603930-132604445' 
ranges.show2 <- StringToGRanges(peak2)
ranges.show2$color <- "royalblue"
  cov_plot2 <- CoveragePlot(
    object = myeloid,
    region.highlight = ranges.show2,
    region = peak2,
    #features = "RAD50",
    annotation = TRUE,
    extend.upstream = 1000,
    extend.downstream  = 1000,
    group.by = 'group', 
    peaks = FALSE
  ) &
    scale_fill_manual( values = group_col ) &
    theme(text = element_text(size=10),
          legend.position = 'none') &
    guides(fill="none") &
    theme_nb()
  cov_plot2
  
  cov_plot1 |cov_plot2 + plot_layout(guides = "collect")
  
  
  
  chr12-66212986-66213443
  chr16-79379258-79379415
  il15 <- chr4-141645655-141645941
  GACAT3 <- chr2-16148619-16230261
  kdm2b <- chr12-121570550-121570768
  cflar <- chr2-201135337-201135742
  ano6 <- chr12-45232871-45233355
  
  peak3 <-'chr12-121570550-121570768' 
  ranges.show3 <- Signac::StringToGRanges(peak3)
  ranges.show3$color <- "coral"
    
  cov_plot3 <- Signac::CoveragePlot(
    object = myeloid,
    region.highlight = ranges.show3,
    region = peak3,
    #features = "RAD50",
    annotation = TRUE,
    extend.upstream = 1000,
    extend.downstream  = 1000,
    group.by = 'group', 
    peaks = FALSE
  ) &
    scale_fill_manual( values = group_col ) &
    theme(text = element_text(size=10),
          legend.position = 'none') &
    guides(fill="none") &
    theme_nb()
  cov_plot3
  
  cov_plot1|cov_plot3
  ggsave_publication('coverage_plot_myeloid', 
                     width = 3.8*2.54, height= 2.5*2.54)  
  
  
  cowplot::plot_grid(cov_plot1 , 
                     cov_plot2, 
                     ncol = 1, 
                     nrow = 2)
  ggsave_publication('coverage_plot_nb2', 
                     width = 12, height= 15)
  (cov_plot1 & theme_nb() )|
    (cov_plot2 & theme_nb() ) 
  ggsave_publication('coverage_plot_nb2', 
                     width = 20, height= 8)
  
  (cov_plot1 & theme(legend.position = "none")) |(cov_plot3  & theme(legend.position = "none"))
  ggsave_publication('coverage_plot_nb_kdm2b', 
                     width = 10, height= 7)
  
  
  cov_plot_sample <- CoveragePlot(
    object = nblast,
    region.highlight = ranges.show,
    region = peak,
    #features = "NFKB1",
    annotation = TRUE,
    extend.upstream = 5000,
    extend.downstream  = 5000,
    group.by = 'sample_id', 
    peaks = FALSE
  ) &
    scale_fill_manual( values = patients_col ) &
    theme(text = element_text(size=40))
  
  
  .draw_coverage_plot(nblast, "chr17-60609440-60609846", 'group', group_col)
  # Save Coverage plots -----------------------------------------------------
  
  save_coverage_plots <- function(dataframe, path ){
    .draw_coverage_plot <- function(object, peak, group_by, color_scheme) {
      peak <- peak
      ranges.show <- StringToGRanges(peak)
      ranges.show$color <- "coral"
        cov_plt <- 
          CoveragePlot(
            object = object,
            region.highlight = ranges.show,
            region = peak,
            annotation = TRUE,
            extend.upstream = 10000,
            extend.downstream  = 10000,
            group.by = group_by , #'group' 
            peaks = FALSE
          ) &  
          scale_fill_manual( values = color_scheme ) &
          theme(text = element_text(size=28)) 
        return(cov_plt)
    }
    genes_up  <- 
      dataframe %>% 
      filter(!is.na(pathway)) %>% 
      group_by(group, pathway) %>% 
      slice_max(order_by = logFC_ATAC, n = 10) 
    genes_up_list <- genes_up$peak
    genes_down <- 
      dataframe %>% 
      filter(!is.na(pathway)) %>%
      group_by(group, pathway) %>% 
      slice_min(order_by = logFC_ATAC, n = 10) 
    genes_down_list <- genes_down$peak
    
    pat = path
    
    for (i in 1:length(genes_down_list)){
      peak_dn <- genes_down_list[i]
      print(peak_dn)
      filename <- 
        paste0(genes_down[i, ]$pathway, '_' , genes_down[i, ]$gene ,"_dn.pdf")
      peak_dn_aggr <- .draw_coverage_plot(nblast, peak_dn, 'group', group_col)
      peak_dn_sample <- .draw_coverage_plot(nblast, peak_dn, 'sample_id', patients_col)
      
      pdf(file.path(pat, filename), width = 40, height = 20)
      print(peak_dn_aggr)
      print(peak_dn_sample)
      dev.off()
    }
  }
  save_coverage_plots(monocytes_atac_df,
                      "/media/AGFORTELNY/people/rohit/projects/neuroblastoma/browser_tracks")
  # filename_aggr <- 
  #   paste0(genes_up[i, ]$pathway, '_' , genes_up[i, ]$gene ,"_aggr.pdf")
  # filename_sample <- 
  #   paste0(genes_up[i, ]$pathway, '_' ,genes_up[i, ]$gene ,"_sample.pdf")
  
  no_gene <-  c("chr12-92603948-92604523", 
                "chr12-92551861-92552230",
                "chr2-152200475-152200993",
                "chr1-66464406-66465145",
                "chr12-92559040-92559218",
                "chr17-34203317-34203629")
  nfkb1 <- 'chr4-102521897-102522522'
  manba <- "chr4-102627693-102628315"
  itpr2 <- "chr12-26334745-26335048"
  hipk2 <- "chr7-139748889-139749340"
  cdc37 <- "chr19-10386069-10386429"
  oxsr1 <- "chr3-38185669-38186114"
  
  # tile_plot <- TilePlot(
  #   object = nblast,
  #   group.by = 'group',
  #   region = 'chr4-102521897-102522522',
  # )
  # tile_plot
  # 
  # expr_plot <- ExpressionPlot(
  #   object = nblast,
  #   features = "NFKB1",
  #   group.by = 'group',
  #   assay = "RNA"
  # )
  # expr_plot
  
  # dotplotdb <- function(seuratObj, peakAnnotation, geneList, name) {
  #   require(matrixStats)
  #   require(stringr)
  #   require(tidyr)
  #   require(dplyr)
  #   bar_cluster <- split(row.names(seuratObj@meta.data), 
  #                        seuratObj@meta.data$seurat_clusters)
  #   
  #   peakAnnotation$peak2 <- gsub("_", "-", peakAnnotation$peak)
  #   
  #   poi <- peakAnnotation %>% 
  #     dplyr::filter(gene %in% geneList) %>% 
  #     dplyr::filter(peak2 %in% row.names(seuratObj@assays$peaks@data)) %>%
  #     #dplyr::filter(peak2 %in% clusterPeaks) %>% 
  #     dplyr::filter(abs(as.numeric(distance)) < 5000)
  #   
  #   m <- nblast@assays$peaks@data[poi$peak2,]
  #   
  #   rmean <- sapply(bar_cluster, function(bars){
  #     Matrix::rowMeans(m[,bars]!= 0) })
  #   row.names(rmean) <- with(poi, paste(gene, peak_type))
  #   
  #   rpercentage <- sapply(bar_cluster, function(bars){
  #     100 * Matrix::rowSums(m[,bars] != 0)/length(bars) })
  #   row.names(rpercentage) <- with(poi, paste(gene, peak_type))
  #   
  #   avg <- rmean %>% 
  #     data.frame() %>% 
  #     mutate(GENE = str_extract(row.names(.), '[[:alnum:]]+')) %>%
  #     group_by(GENE) %>% 
  #     summarise(across(everything(), mean)) %>% 
  #     data.frame() %>% 
  #     `rownames<-`(.[,'GENE']) %>% select(-GENE) %>% 
  #     mutate(mean_all = rowMeans(.), stdev=rowSds(as.matrix(.))) %>% 
  #     #filter(stdev > 0.01 & mean_all >0.05) %>% 
  #     tibble::rownames_to_column("GENE") %>% 
  #     rename_with(., ~ gsub("X", "", .x, fixed = TRUE)) %>% 
  #     tidyr::pivot_longer(cols = -c(GENE, mean_all, stdev), 
  #                         names_to = "Cluster", values_to = "Avg") 
  #   
  #   pct <- rpercentage %>% 
  #     data.frame() %>% 
  #     mutate(GENE = str_extract(row.names(.), '[[:alnum:]]+')) %>% 
  #     group_by(GENE) %>% 
  #     summarise(across(everything(), mean)) %>% 
  #     data.frame() %>% 
  #     `rownames<-`(.[,'GENE']) %>% select(-GENE) %>% 
  #     mutate(mean_all = rowMeans(.), stdev=rowSds(as.matrix(.))) %>% 
  #     #filter(stdev > 0.1 & mean_all >1) %>% 
  #     tibble::rownames_to_column("GENE") %>% 
  #     rename_with(., ~ gsub("X", "", .x, fixed = TRUE)) %>% 
  #     tidyr::pivot_longer(cols = -c(GENE, mean_all, stdev), 
  #                         names_to = "Cluster", values_to = "Pct")  
  #   
  #   merdb <- merge(avg, pct, by = c('GENE', "Cluster")) %>% 
  #     dplyr::mutate(Name = name) 
  #   return(merdb)
  # }

# scratch -----------------------------------------------------------------

# umap_label <- 
#   FetchData(nblast, 
#             vars = c("ident", "UMAP_1", "UMAP_2", "orig.ident")) %>% 
#   data.table::setDT(.) ->xDT
# # xDT is a data.table with the following columns: UMAP_1, UMAP_2, Clusters, rn (the cell barcode)
# hex.obj <- hexbin::hexbin(x=xDT$UMAP_1, y=xDT$UMAP_2, xbins = 100, IDs=TRUE)
# pDT <- cbind(xDT, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
# pDT <- pDT[,.N, by=c("hex.x", "hex.y", "ident")]
# pDT[, sum := sum(N), by=c("hex.x", "hex.y")]
# pDT[, frac := N / sum]
# pDT <- pDT[order(frac, decreasing = TRUE)][,head(.SD, 1), by=c("hex.x", "hex.y")]
# 
# ggplot(pDT, aes(x=hex.x, y=hex.y, color=ident)) + 
#   geom_point() + 
#   labs(x = "UMAP_1", y = "UMAP_2", 
#        color = "") + 
#   scale_color_manual(values = CELL_TYPE_COLORS) +
#   theme_nb() 
# ggsave_publication('umap_wes_theme_nb_light', 
#                    width = 10, height= 8)
# # 'TNFA' = 'royalblue', 'MYC' = 'gold',
# # 'Inflammation' = 'chartreuse')
# big_data %>% 
#   filter(cell_type == "Cluster_Monocytes") %>% 
#   View()
# 
# Seurat::Idents(object=nblast) <- nblast$new.ident
# 
# # Add patient groups information
# 
# nblast@meta.data$group <- nblast@meta.data %>% 
#   dplyr::mutate(
#     group = dplyr::case_when(orig.ident %in% groups[[1]] ~"C",
#                              orig.ident %in% groups[[2]] ~"M", 
#                              orig.ident %in% groups[[3]] ~"A",
#                              orig.ident %in% groups[[4]] ~"S")) %>% 
#   pull(group)
# nblast$group <- nblast@meta.data$group
# 
# # Add patient information with group information
# nblast@meta.data$sample_id <-
#   nblast@meta.data %>% 
#   dplyr::mutate(
#     sample_id = dplyr::case_when(orig.ident %in% groups[[1]] ~ paste0("I-",new.ident),
#                              orig.ident %in% groups[[2]] ~paste0("II-",new.ident), 
#                              orig.ident %in% groups[[3]] ~paste0("III-",new.ident),
#                              orig.ident %in% groups[[4]] ~paste0("IV-",new.ident))) %>% 
#   pull(sample_id)
# nblast$sample_id <- nblast@meta.data$sample_id
# 
# # Rename Monocytes to Myeloids cells
# nblast <- RenameIdents(object = nblast, `Monocytes` = "Myeloid")
# 
# cellMarkers <- list()
# cellMarkers$t_cell <- c('CCR4', 'CCR6')
# cellMarkers$b_cell <- c('CD10', 'CD19', 'CD20','CD138')
# cellMarkers$monocytes <- c('CD14', 'CD16', 'CD31')
# cellMarkers$erythrocyte <- c('ALDH')
# cellMarkers$myeloid <- c('CD10', 'CD106', 'CD115', 'CD117', 'CD11b', 
#                          'CD11c', 'CD122', 'CD123', 'CD13', 'CD130', 
#                          'CD14', 'CD141', 'CD15', 'CD16', 'CD163', 
#                          'CD165', 'CD169', 'CD177', 'CD178', 'CD183', 
#                          'CD192', 'CD193', 'CD194', 'CD195', 'CD198', 
#                          'CD200R', 'CD203c', 'CD205', 'CD206', 'CD212', 
#                          'CD217', 'CD218a', 'CD23', 'CD24', 'CD244', 
#                          'CD25', 'CD282', 'CD284', 'CD289', 'CD294', 
#                          'CD3', 'CD31', 'CD32', 'CD33', 'CD34', 'CD38', 
#                          'CD4', 'CD40', 'CD43', 'CD45', 'CD45RA', 'CD45RO', 
#                          'CD48', 'CD49d', 'CD49f', 'CD51', 'CD62L', 'CD66b', 
#                          'CD68', 'CD69', 'CD78', 'CD80', 'CD81', 'CD84', 
#                          'CD86', 'CD9', 'CXCR3')
# cellMarkers$dendridic <- c('CD1a','CD1c' ,'CD11c', 'SIRPA')
# cellMarkers$stromal <- c('CD90', 'CD105', 'CD73', 'CD34', 'CD44', 'CD45', 'CD271')
# cellMarkers$progenitor <- c('CD34', 'CD133')
# cellMarkers$macrophage <-c('CD206', 'CD163')
# cellMarkers$stem <- c('CD105', 'CD90', 'CD73', 'CD34')
# 
# plot_VlnPlot <- function(feature, title){
#   VlnPlot( object = nblast,
#            features = feature , pt.size = 0, 
#            split.by = 'orig.ident', ncol = 2,
#           combine = FALSE
#   )  & 
#     #coord_flip() & 
#     theme_nb() &
#     labs(title = title) &
#     theme(legend.title = element_blank(),
#           legend.position = "none") &
#     scale_fill_manual(
#       values = sample_colors)
# }
# qc_after1 <- plot_VlnPlot('pct_reads_in_peaks', 'Percent reads in peaks')
# qc_after2 <- plot_VlnPlot('peak_region_fragments', 'Fragments in peak region')
# qc_after3 <- plot_VlnPlot('TSS.enrichment',  'TSS Enrichment')
# qc_after4 <- plot_VlnPlot('blacklist_ratio',  'Encode blacklisted ratio')
# qc_after5 <- plot_VlnPlot('nucleosome_signal','Nuleosome signal score')


#scale_y_continuous(limits = c(1,1.015)) &
# theme(legend.position = 'none',
#       text = element_text(size=18, face = "bold"),
#       axis.text.x = element_text(angle=0,  hjust = 1, size=12, color="black" ),
#       axis.text.y = element_text( size=12, color="black" ),
#       axis.title    = element_text( size=15, face = "bold" ),
#       legend.text = element_text( size=40, color="black"),
#       legend.title = element_text( size=15, color="black"),
#       axis.line.x = element_line(size=1.5, color="black"),
#       axis.line.y = element_line(size=1.5, color="black"))

# plot_list <- plot_VlnPlot
# plot_list[[6]] <- umap_combined + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.5))
# plot_list[[7]] <- umap_atac + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.5))
# plot_list[[8]] <- umap_rna + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(colour = "black", size=0.5))
# 
# cowplot::plot_grid(plotlist = plot_list,
#                      ncol = 4, labels = "AUTO"
#                      ) 
# plot_VlnPlot2 <- 
#   lapply(plot_VlnPlot, `+`, labs(title = titles) +
#          scale_fill_manual(
#            values = sample_colors) +
#          theme_nb() + 
#          theme(legend.title = element_blank(),
#                legend.position = "none") +
#          scale_fill_manual(
#            values = sample_colors))
### We next perform dimension reduction using scOpen. 

### Note here we have to use reticulate package to allow us call python function in R.
# qc_plots_reviewers <-
#   ((plot_VlnPlot[[1]] + theme(legend.position = 'none') | plot_VlnPlot[[2]]) / 
#     (plot_VlnPlot[[3]] + theme(legend.position = 'none') | plot_VlnPlot[[4]]) / 
#     (plot_VlnPlot[[5]] + theme(legend.position = 'none')| plot_spacer())/
#        (tss_score & theme(
#          legend.position="none",
#          plot.title = element_text(size=6, face="bold",family="sans"))|
#           nucleosome_score & theme(
#             legend.position="none",
#             plot.title = element_text(size=6, face="bold",family="sans")))) +
#   plot_annotation(tag_levels = 'a') & 
#   theme(plot.tag = element_text(face="bold", size = 10))
# nblast <- RenameIdents(object = nblast, 
#                        "T-cell" = "T_cell", 
#                        "B-cell" = "B_cell",
#                        "Memory_B-cell" = "Memory_B_cell")
# nblast <- RenameIdents(object = nblast, 
#                        `NKT` = "NK")
# nblast <- RenameIdents(object = nblast, 
#                        `Monocytes` = "M")
# 
# nblast <- RenameIdents(object = nblast, 
#                        `T_cell` = "T",
#                        `Monocytes` = "M",
#                        `B_cell` = "B",
#                        `Memory_B_cell` = "MB",
#                        `Erythroblasts` = "E")
# Patient composition to each cell type  ----------------------------------

df_melt <- 
  nblast@meta.data %>% 
  dplyr::select(new.ident, activ.ident  ) %>% 
  dplyr::group_by( new.ident, activ.ident) %>% 
  summarise(no = n()) %>% 
  spread (new.ident, no) %>% View()
  melt( ., id.vars="activ.ident", value.name="Cells", variable.name="Patient" ) 

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 16
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

ggplot(df_melt) +
  geom_col(aes(x=as.factor(activ.ident), y=Cells, fill = Patient)) +
  scale_fill_manual(values = mycolors) +
  labs(title = "Sample composition of each cell type",
       x = "Cell type", y = 'Cells', fill = "Patients") +
  hrbrthemes::theme_ipsum() 

# PatchWork ---------------------------------------------------------------

tmp_cov_plot_1 <- 
  cov_plot1 & 
  theme_nb() & 
  theme(legend.position = "none")
tmp_cov_plot_2 <- 
  cov_plot2 & 
  theme_nb()
cov_plot1 | cov_plot2  
tmp_umap <- umap + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank())
ggsave_publication('tmp_umap', 
                   width = 9, height= 6)


layout <- "
AA####
AA####
BC####
DE####
"
umap  + scatter + 
  IRF1 +theme(legend.position = "none")  + 
  NFKB1 + theme(legend.position = "none") + 
  fos +theme(legend.position = "none") + 
  bach1 + theme(legend.position = "none") + 
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')
# umap +  promoter_dotplot  + distal_dotplot + 

ggsave_publication('FigureLayout_patchwork_trial', 
                   width = 10, height= 20)

scatter <- 
p3_e2f / p3_tnfa + 
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom') 
tmp_umap + p3_e2f + p3_igr + 
  tmp_cov_plot_1 + tmp_cov_plot_2 +  
  IRF1 + IRF2 + 
  fos + bach1 + 
  distal_dotplot + promoter_dotplot +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

ggsave_publication('FigureLayout_patchwork', 
                   width = 20, height= 29)

(tmp_umap |(p3_e2f/p3_igr))/
  (tmp_cov_plot_1 | tmp_cov_plot_2 | plot_spacer() |plot_spacer() ) / 
  ((IRF1 & theme(legend.position = "none") | IRF2 |plot_spacer() |plot_spacer()) /
     (fos & theme(legend.position = "none")| bach1 | plot_spacer()|plot_spacer())) +
  #plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

(umap /(p3_e2f / p3_tnfa/(IRF1 |NFKB1)/(fos|bach1))) | ((promoter_dotplot | distal_dotplot))

ggsave_publication('FigureLayout_patchwork_trial', 
                   width = 20, height= 26)