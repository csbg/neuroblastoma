pacman::p_load(Seurat, Signac, SingleCellExperiment, 
               tidyverse, nebula, scater)

Cluster.neunonal <- readRDS(
  #file = "/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/neuronal_cluster.rds",
  file = "/media/AGFORTELNY/people/rohit/projects/neuronal_cluster.rds")

# add patient group information
groups <- list(c('2', '3', '6', '15', '16'),
               c('1', '4', '5', '11'),
               c('7','12'),
               c('8', '9', '10','13', '14'))

Cluster.neunonal@meta.data$group <- Cluster.neunonal@meta.data %>% 
  dplyr::mutate(
    group = dplyr::case_when(orig.ident %in% groups[[1]] ~"1",
                             orig.ident %in% groups[[2]] ~"2", 
                             orig.ident %in% groups[[3]] ~"3",
                             orig.ident %in% groups[[4]] ~"4")) %>% 
  pull(group)

Cluster.neunonal$group <- Cluster.neunonal@meta.data$group

Cluster.neunonal <- RunTFIDF(Cluster.neunonal, method = 4)
Cluster.neunonal <- FindTopFeatures(Cluster.neunonal, min.cutoff = 'q5')
Cluster.neunonal <- RunSVD(object = Cluster.neunonal)

DepthCor(Cluster.neunonal)

# Non-linear dimension reduction and clustering ---------------------------

Cluster.neunonal <- RunUMAP(
  object = Cluster.neunonal,
  reduction = 'lsi',
  umap.method = 'uwot', 
  dims = 2:30)

Cluster.neunonal <- FindNeighbors(
  object = Cluster.neunonal,
  reduction = 'lsi',
  nn.method = 'rann',
  dims = 2:30)

Cluster.neunonal <- FindClusters(
  object = Cluster.neunonal,
  algorithm = 3,
  resolution = 1,
  verbose = FALSE)


umap_seurat_ident <- DimPlot(Cluster.neunonal)

Cluster.neunonal <- RenameIdents(
  object = Cluster.neunonal,
  '0' = 'clust1', '1' = 'clust1', '2' = 'clust1',  '3' = 'clust1',
  '4' = 'clust2', '5' = 'clust2', '6' = 'clust2'
)

umap_assign_ident <- DimPlot(Cluster.neunonal)

sce <- SeuratWrappers::as.cell_data_set(Cluster.neunonal)

#sce <- createSingleCellObject(Cluster.neunonal)

sce <- logNormCounts(sce) #, assay.type = "ln_counts")

model_mat <- model.matrix(~group , data = colData(sce))
model_mat <- model_mat[, colSums(model_mat != 0) > 0] 

sce_grouped <- group_cell(counts(sce), 
                          id = colData(sce)$orig.ident, 
                          pred = model_mat, 
                          offset = NULL)

result <- nebula(sce_grouped$count, 
                 sce_grouped$id, 
                 pred=sce_grouped$pred)

.get_chipseeker_annot <- function(peaks) {
  require(Signac)
  require(ChIPseeker)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  granges <- StringToGRanges(peaks)
  peakAnno <- 
    annotatePeak(granges, 
                 tssRegion=c(-3000, 3000),
                 TxDb=txdb, 
                 annoDb="org.Hs.eg.db")
  
  chipseeker_annot <- 
    peakAnno@anno %>% 
    Repitools::annoGR2DF() %>% 
    dplyr::mutate(gene = paste0(chr,'-' ,start,'-', end)) 
  return(chipseeker_annot)
}


result$summary <- result$summary %>%
  as_tibble() %>% 
  mutate(
    algorithm = result$algorithm,
    convergence = result$convergence,
    overdispersion_subject = result$overdispersion$Subject,
    overdispersion_cell = result$overdispersion$cell
  )
closest_genes <- .get_chipseeker_annot(result$summary$gene)

closest_genes <- ClosestFeature(Cluster.neunonal, 
                                regions = result$summary$gene) %>% 
  mutate(gene = query_region)

results_with_closest_gene <- merge(result$summary, 
                                   closest_genes, by = 'gene')
results_with_closest_gene <- 
  results_with_closest_gene %>% 
  dplyr::mutate(q_group3 = qvalue::qvalue(results_with_closest_gene$p_group3)$qvalues) %>%
  dplyr::mutate(q_group4 = qvalue::qvalue(results_with_closest_gene$p_group4)$qvalues) %>% 
  dplyr::mutate(logFC_group3 = logFC_group3/log(2),
                logFC_group4 = logFC_group4/log(2)) 
myfile = '/media/AGFORTELNY/PROJECTS/Neuroblastoma/analysis/rohit/nbatac/data/nebula/nb.csv'
write.table(results_with_closest_gene, file = myfile)

results_with_closest_gene <- read.table(myfile, header= TRUE, sep = ' ')
results_with_closest_gene %>% 
  filter(p_group3 < 0.05) %>% dim()

# old_result <- readRDS(file = '/media/AGFORTELNY/people/rohit/projects/nebula_neuroblastoma.rds')

# Enrichment analysis

results <- 
  read.table('/media/AGFORTELNY/PROJECTS/Neuroblastoma/analysis/rohit/nbatac/data/nebula/nb.csv')

grp3 <-
  results %>% 
  dplyr::filter(q_group3 < 0.5) %>% 
  dplyr::mutate(logFC = logFC_group3) %>% 
  dplyr::select(logFC, SYMBOL) %>% 
  dplyr::arrange(desc(logFC))

grp4 <-
  results %>% 
  dplyr::filter(q_group4 < 0.5) %>% 
  dplyr::mutate(logFC = logFC_group4) %>% 
  dplyr::select(logFC, SYMBOL) %>% 
  dplyr::arrange(desc(logFC))

grp3_genes <- grp3$logFC
names(grp3_genes) <- grp3$SYMBOL
eg <-  
  clusterProfiler::bitr(
    names(grp3_genes), 
    fromType="SYMBOL", 
    toType="ENTREZID", 
    OrgDb="org.Hs.eg.db")
head(eg$ENTREZID)

genes_entrez <- grp3 %>% 
  left_join(eg, by = "SYMBOL") %>% 
  dplyr::select(logFC, ENTREZID) %>% 
  arrange(desc(logFC))

geneList <- genes_entrez$logFC
names(geneList) <- genes_entrez$ENTREZID

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

ego <- enrichGO(gene          = names(grp3_genes),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "ALL",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)

ego2 <- gseDO(grp3_genes,
              keyType       = 'SYMBOL',
              ont           = "ALL",
              pvalueCutoff  = 1,
              qvalueCutoff  = 1)
enp <- enrichPathway(gene=eg$ENTREZID, 
                     pvalueCutoff = 0.5, readable=TRUE)

gsep <- gsePathway(geneList, 
                pvalueCutoff = 1,
                pAdjustMethod = "BH", 
                verbose = FALSE)
#
barplot(ego, showCategory=20) 
dotplot(ego, showCategory=20) + ggtitle("dotplot for GO")
#
barplot(enp)
#
ridgeplot(gsep)
cnetplot(gsep, foldChange=geneList)

gseaplot2(gsep, geneSetID = 1:10, subplots = 1)
gsearank(gsep, 1, title = gsep[1, "Description"])



#plotGOgraph(gsep)






