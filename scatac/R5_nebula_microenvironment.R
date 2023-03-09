
# Create big_data data frame ----------------------------------------------

# Mounted on local machine 
path = '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj'

# Singularity mounted on CAME 
path = '/media/AGFORTELNY/people/rohit/projects/clusterObj'

fnames = list.files(path, pattern="csv")

nebula_results_df <- function(filename) {
  cell_type <- str_replace(
    basename(filename), ".csv", "")
  results_with_closest_gene <- read.table(filename, header= TRUE, sep = ' ')
  results_with_closest_gene <- results_with_closest_gene %>% 
    
    dplyr::mutate(q_groupII = qvalue::qvalue(results_with_closest_gene$p_groupII)$qvalues) %>% 
    dplyr::mutate(q_groupIII = qvalue::qvalue(results_with_closest_gene$p_groupIII)$qvalues) %>%
    dplyr::mutate(q_groupIV = qvalue::qvalue(results_with_closest_gene$p_groupIV)$qvalues) %>% 
    dplyr::mutate(logFC_groupII = logFC_groupII/log(2),
           logFC_groupIII = logFC_groupIII/log(2),
           logFC_groupIV = logFC_groupIV/log(2)) %>% 
    dplyr::filter(!-log10(q_groupII) > 10) %>% 
    dplyr::filter(!-log10(q_groupIII) > 10) %>%
    dplyr::filter(!-log10(q_groupIV) > 10) %>%
    dplyr::mutate(significanceII = case_when(
      (logFC_groupII < -1 & q_groupII < 0.05 )~ 'Down',
      (logFC_groupII > 1 & q_groupII < 0.05 )~ 'Up', TRUE ~ 'No change')) %>%
    dplyr::mutate(significanceIII = case_when(
      (logFC_groupIII < -1 & q_groupIII < 0.05 )~ 'Down',
      (logFC_groupIII > 1 & q_groupIII < 0.05 )~ 'Up', TRUE ~ 'No change')) %>%
    dplyr::mutate(significanceIV = case_when(
      (logFC_groupIV < -1 & q_groupIV < 0.05 )~ 'Down',
      (logFC_groupIV > 1 & q_groupIV < 0.05 )~ 'Up', TRUE ~ 'No change')) %>%
    dplyr::mutate(delabel = case_when(significanceIII != 'No change' ~ gene_name,
                               significanceII != 'No change' ~ gene_name,
                               significanceIV != 'No change' ~ gene_name)) %>% 
    dplyr::mutate(cell_type = paste0('Cluster_',cell_type))
  return(results_with_closest_gene)
}

file_list = list()
for (i in file.path(path, fnames)) {
  cluster_name = paste0('Cluster_',
                        sub(pattern = "(.*)\\..*$", 
                     replacement = "\\1", 
                     basename(i)))
  file_list[[cluster_name]]  = nebula_results_df(i)
}
big_data = do.call(rbind, file_list)

# Section_1 ---------------------------------------------------------------

# test 
big_data %>% 
  filter(cell_type == "Cluster_Monocytes") %>% 
  filter(gene_name == "BRMS1L") %>% View()

all_da_peaks_granges <- 
  big_data %>% 
  pull(gene) %>% 
  unique() %>% 
  Signac::StringToGRanges() 

peak_bed <- data.frame(
  seqnames = seqnames(all_da_peaks_granges),
  start = start(all_da_peaks_granges),
  end = end(all_da_peaks_granges),
  strand = strand(all_da_peaks_granges))

peak_bed <- peak_bed %>% 
  dplyr::select(seqnames, start, end )
write.table(peak_bed, 
            '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/bigwig/peaks.bed', 
             quote = FALSE, sep="\t", row.names = FALSE)

# homer_motif_analysis ----------------------------------------------------

run_gsea <- function(dataframe) {
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2)
  require(EnsDb.Hsapiens.v75)
  require(org.Hs.eg.db)
  require(DOSE)
  organism = org.Hs.eg.db
  original_gene_list <- dataframe$logFC
  # name the vector
  names(original_gene_list) <- dataframe$gene_name
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="MF", 
               keyType = "GENENAME", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.5, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  return(gse)
}
run_msigdb <- function(dataframe) {
  # require(clusterProfiler)
  # require(enrichplot)
  # require(ggplot2)
  # require(EnsDb.Hsapiens.v75)
  # require(org.Hs.eg.db)
  require(msigdbr)
  require(DOSE)
  
  original_gene_list <- dataframe$logFC
  # name the vector
  names(original_gene_list) <- dataframe$gene_name
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol) %>% 
    dplyr::mutate(gs_name = str_replace(gs_name, "HALLMARK_", ""))  %>% 
    dplyr::mutate(gs_name = str_replace_all(gs_name, "_", " ")) 
  em2 <- GSEA(gene_list, TERM2GENE = H_t2g, pvalueCutoff = 1)
  return(em2)
}
run_gprofiler <- function(namedList){
  require(gprofiler2)
  gostres <- gost(query = namedList, 
                  organism = "hsapiens", ordered_query = TRUE, 
                  multi_query = TRUE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
}

table(big_data$cell_type)
# peak, logFC, gene, annotation, distanceToTSS, p_value, q_value, comparison
bigdf_chipseeker <- list()




for (i in levels(as.factor(big_data$cell_type))){
  print(i)
}



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

.get_promoter_peaks <- function(granges) {
  require(ChIPseeker)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  peakAnno <- 
    annotatePeak(granges, 
                 tssRegion=c(-3000, 3000),
                 TxDb=txdb, 
                 annoDb="org.Hs.eg.db")
  message('Promoter peaks for MYCN vs Control')
  promoter_peaks <- 
    peakAnno@anno %>% 
    Repitools::annoGR2DF() %>% 
    dplyr::filter(str_detect(annotation, "Promoter|intron 1 of"))  %>% 
    dplyr::mutate(peak = paste0(chr,'-' ,start,'-', end)) %>% 
    dplyr::pull(peak) 
  return(promoter_peaks)
}
.filter_df_for_promoters <- function(dataframe, promoters_lst) {
    dataframe %>%  
    dplyr::filter(gene %in% promoters_lst)
}
.save_bed_files <- function(dataframe, path, string) {
  write.table(dataframe, 
              file = file.path(path, paste0(string, ".bed")), 
              row.names = FALSE, sep = '\t', quote = FALSE)
}
.create_granges_from_peak_df <- function(df, colname){
  df %>% 
    dplyr::pull(colname) %>% 
    Signac::StringToGRanges() 
}
.group_filtered_df <- function(df, celltype, group,  logfc) {
  df %>% 
    dplyr::filter(cell_type == celltype,
                  q_groupII <= 0.05, 
                  case_when((group == "MYCN" & logfc == 'Up') ~ logFC_groupII > 0,
                            (group == "MYCN" & logfc == 'Down') ~ logFC_groupII < 0,
                            (group == "ATRX" & logfc == 'Up') ~ logFC_groupIII > 0,
                            (group == "ATRX" & logfc == 'Down')~ logFC_groupIII < 0,
                            (group == "Sporadic" & logfc == 'Up') ~ logFC_groupIV > 0,
                            (group == "Sporadic" & logfc == 'Down') ~ logFC_groupIV < 0)
                  )
}  



cell_type_go_data <- function(big_data_frame, celltype){
  #
  gprofiler_list <-  list()
  message('Running MYCN vs Control...')
  message('Upregulated ')
  mycn_ctrl_up <- .group_filtered_df(big_data,'Cluster_Monocytes',"MYCN", "Up")
  message('Get promoter peaks..')
  mycn_ctrl_up_promoters <- 
    .get_promoter_peaks(
    .create_granges_from_peak_df(mycn_ctrl_up,'gene'))
  message('Get promoter filtered Upregulated')
  mycn_up <- 
    .filter_df_for_promoters(mycn_ctrl_up, 
                             mycn_ctrl_up_promoters) %>% 
    dplyr::select(gene_name, logFC_groupII) %>% 
    dplyr::mutate(logFC = logFC_groupII) %>% 
    arrange(logFC)
  gprofiler_list[['MYCN Open']] <- mycn_up$gene_name
  message('Downregulated ')
  mycn_ctrl_down <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupII <= 0.05,
                  logFC_groupII < 0)  
  #
  mycn_ctrl_down <- .group_filtered_df(big_data,'Cluster_Monocytes',"MYCN", "Down")
  message('Get promoter peaks..')
  mycn_ctrl_down_promoters <- 
    .get_promoter_peaks(
      .create_granges_from_peak_df(mycn_ctrl_down,'gene'))
  message('Get promoter filtered Downregulated')
  mycn_down <- 
    .filter_df_for_promoters(mycn_ctrl_down, 
                             mycn_ctrl_down_promoters) %>% 
    dplyr::select(gene_name, logFC_groupII) %>% 
    dplyr::mutate(logFC = logFC_groupII) %>% 
    arrange(logFC)
  gprofiler_list[['MYCN Close']] <- mycn_down$gene_name
  #
  #
  message('Running ATRX vs Control...')
  message('Upregulated ')
  atrx_ctrl_up <- .group_filtered_df(big_data,'Cluster_Monocytes',"ATRX", "Up")
  message('Get promoter peaks..')
  atrx_ctrl_up_promoters <- 
    .get_promoter_peaks(
      .create_granges_from_peak_df(atrx_ctrl_up,'gene'))
  message('Get promoter filtered Upregulated')
  atrx_up <- 
    .filter_df_for_promoters(atrx_ctrl_up, 
                             atrx_ctrl_up_promoters) %>% 
    dplyr::select(gene_name, logFC_groupIII) %>% 
    dplyr::mutate(logFC = logFC_groupIII) %>% 
    arrange(logFC)
  gprofiler_list[['ATRX Open']] <- atrx_up$gene_name
  message('Downregulated ')

  #
  atrx_ctrl_down <- .group_filtered_df(big_data,'Cluster_Monocytes',"ATRX", "Down")
  message('Get promoter peaks..')
  atrx_ctrl_down_promoters <- 
    .get_promoter_peaks(
      .create_granges_from_peak_df(atrx_ctrl_down,'gene'))
  message('Get promoter filtered Downregulated')
  atrx_down <- 
    .filter_df_for_promoters(atrx_ctrl_down, 
                             atrx_ctrl_down_promoters) %>% 
    dplyr::select(gene_name, logFC_groupIII) %>% 
    dplyr::mutate(logFC = logFC_groupIII) %>% 
    arrange(logFC)
  gprofiler_list[['ATRX Close']] <- atrx_down$gene_name
  #
  #
  message('Running Sporadic vs Control...')
  message('Upregulated ')
  sporadic_ctrl_up <- .group_filtered_df(big_data,'Cluster_Monocytes',"Sporadic", "Up")
  message('Get promoter peaks..')
  sporadic_ctrl_up_promoters <- 
    .get_promoter_peaks(
      .create_granges_from_peak_df(sporadic_ctrl_up,'gene'))
  message('Get promoter filtered Upregulated')
  sporadic_up <- 
    .filter_df_for_promoters(sporadic_ctrl_up, 
                             sporadic_ctrl_up_promoters) %>% 
    dplyr::select(gene_name, logFC_groupIV) %>% 
    dplyr::mutate(logFC = logFC_groupIV) %>% 
    arrange(logFC)
  gprofiler_list[['Sporadic Open']] <- sporadic_up$gene_name
  message('Downregulated ')
  
  #
  sporadic_ctrl_down <- .group_filtered_df(big_data,'Cluster_Monocytes',"Sporadic", "Down")
  message('Get promoter peaks..')
  sporadic_ctrl_down_promoters <- 
    .get_promoter_peaks(
      .create_granges_from_peak_df(sporadic_ctrl_down,'gene'))
  message('Get promoter filtered Downregulated')
  sporadic_down <- 
    .filter_df_for_promoters(sporadic_ctrl_down, 
                             sporadic_ctrl_down_promoters) %>% 
    dplyr::select(gene_name, logFC_groupIV) %>% 
    dplyr::mutate(logFC = logFC_groupIV) %>% 
    arrange(logFC)
  gprofiler_list[['Sporadic Close']] <- sporadic_down$gene_name
  
  return(gprofiler_list)
}

monocytes_df <- cell_type_go_data(
  big_data, 
  'Cluster_Monocytes') 

monocytes_profile <- 
  run_gprofiler(monocytes_df)
library(gprofiler2)
gostplot(monocytes_profile, 
         capped = TRUE, 
         interactive = FALSE)

monocytes_profile$result %>% 
  filter(source =='GO:BP') %>% 
  arrange('term_size') %>% 
  #head(30) %>% 
  pull(term_id) ->highlight_terms_gobp

c('TF', highlight_terms_tf)
c('GO:CC', highlight_terms_gocc)
c('GO:MF', highlight_terms_gomf)
c('GO:BP', highlight_terms_gobp)

apply(monocytes_profile$result$significant, 1, as.numeric)
as.integer(monocytes_profile$result$significant)

monocytes_profile$result$significant[monocytes_profile$result$significant=='TRUE']


publish_gosttable(
  monocytes_profile,
  highlight_terms = highlight_terms_gobp,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size"),
  filename='/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/images/January/gobp.pdf',
  #ggplot = TRUE
)

run_gsea(testdf)
dotplot(run_msigdb(testdf))


cell_type_bed_files <- function(cell_file) {
  pacman::p_load(plyr,readr, tools)
  mydir = "/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj"
  .create_save_bed <- function(dataframe, outfile) {
    df_gr <- 
      Signac::StringToGRanges(unique(dataframe$gene))
    
    peak_bed <- data.frame(
      seqnames = seqnames(df_gr),
      start = start(df_gr),
      end = end(df_gr),
      strand = strand(df_gr))
    
    peak_bed$strand <- 
      ifelse(peak_bed$start < peak_bed$end, '+','-')
    
    peak_bed$id <- 
      paste0(as.character(peak_bed$seqnames), '_',
             as.character(peak_bed$start),'_',
             as.character(peak_bed$end) )
    peak_bed <- 
      peak_bed %>% 
      select(id, seqnames, start, end, strand )
    write.table(peak_bed, file.path(mydir, outfile), 
                quote = FALSE, sep="\t", row.names = FALSE)
  }
  
  dat_csv1 <-
    read.csv(file.path(mydir, cell_file ), sep = ' ', header = TRUE)
  outstr <- 
    paste0('/group_microenv/',
           tolower(file_path_sans_ext(cell_file)))
  dat_csv1 %>% 
    filter(logFC_groupII < 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr, '_mycn_down.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIII < 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr, '_atrx_down.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIV < 0) %>% 
    .create_save_bed(., 
                     outfile= paste0(outstr, '_sporadic_down.bed'))
  dat_csv1 %>% 
    filter(logFC_groupII > 0) %>% 
    .create_save_bed(., 
                     outfile= paste0(outstr, '_mycn_up.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIII > 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr, '_atrx_up.bed'))
  dat_csv1 %>% 
    filter(logFC_groupIV > 0) %>% 
    .create_save_bed(., 
                     outfile = paste0(outstr,'_sporadic_up.bed'))
}
for (i in c('Monocytes.csv', 'Basophils.csv', 'T-cells.csv',
            'B_cells_memory_34.csv', 'B-cells.csv', '37.csv', 
            '36.csv', '33.csv', '32.csv', 
            '26.csv', '22.csv', '11.csv')) {
  cell_type_bed_files(i) 
}



beds_from_peak_df <- function(big_data_frame, cell_type, path ) {
  
  .get_promoter_peaks <- function(granges) {
    require(ChIPseeker)
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peakAnno <- 
      annotatePeak(granges, 
                   tssRegion=c(-3000, 3000),
                   TxDb=txdb, 
                   annoDb="org.Hs.eg.db")
    message('Promoter peaks for MYCN vs Control')
    promoter_peaks <- 
      peakAnno@anno %>% 
      Repitools::annoGR2DF() %>% 
      dplyr::filter(str_detect(annotation, "Promoter|intron 1 of"))  %>% 
      dplyr::mutate(peak = paste0(chr,'_' ,start,'_', end)) %>% 
      dplyr::pull(peak)
    return(promoter_peaks)
  }
  .create_bed_df <- function(granges, promoters_lst) {
    df <-
      data.frame(seqnames = seqnames(granges),
                 start = start(granges), 
                 end = end(granges), 
                 strand = strand(granges)
      ) %>% 
      dplyr::mutate(id = paste0(seqnames,'_',start, '_', end) ) %>% 
      dplyr::filter(id %in% promoters_lst) %>% 
      dplyr::select(id, seqnames, start, end, strand )%>% 
      dplyr::mutate(strand = ifelse(as.numeric(start) < as.numeric(end), 
                                    '+', '-'))
  }
  .save_bed_files <- function(dataframe, path, string) {
    write.table(dataframe, 
                file = file.path(path, paste0(string, ".bed")), 
                row.names = FALSE, sep = '\t', quote = FALSE)
  }
  
  message('Running MYCN vs Control...')
  mycn_ctrl_up <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupII <= 0.05, 
                  logFC_groupII > 0) %>% 
    pull(gene) %>% Signac::StringToGRanges() 
  mycn_ctrl_down <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupII <= 0.05,
                  logFC_groupII < 0) %>% 
    pull(gene) %>% Signac::StringToGRanges() 
  mycn_ctrl_up_promoters <- .get_promoter_peaks(mycn_ctrl_up)
  mycn_ctrl_up_df <- .create_bed_df(mycn_ctrl_up, mycn_ctrl_up_promoters)
  .save_bed_files(mycn_ctrl_up_df, path, 'monocytes_up_mycn_30_12')
  mycn_ctrl_down_promoters <- .get_promoter_peaks(mycn_ctrl_down)
  mycn_ctrl_down_df <- .create_bed_df(mycn_ctrl_down, mycn_ctrl_down_promoters)
  .save_bed_files(mycn_ctrl_down_df, path, 'monocytes_down_mycn_30_12')
  
  message('Running ATRX vs Control...')
  atrx_ctrl_up <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupIII <= 0.05, 
                  logFC_groupIII > 0) %>% 
    pull(gene) %>% Signac::StringToGRanges() 
  atrx_ctrl_down <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupIII <= 0.05,
                  logFC_groupIII < 0) %>% 
    pull(gene) %>% Signac::StringToGRanges() 
  atrx_ctrl_up_promoters <- .get_promoter_peaks(atrx_ctrl_up)
  atrx_ctrl_up_df <- .create_bed_df(atrx_ctrl_up, atrx_ctrl_up_promoters)
  .save_bed_files(atrx_ctrl_up_df, path, 'monocytes_up_atrx_30_12')
  atrx_ctrl_down_promoters <- .get_promoter_peaks(atrx_ctrl_down)
  atrx_ctrl_down_df <- .create_bed_df(atrx_ctrl_down, atrx_ctrl_down_promoters)
  .save_bed_files(atrx_ctrl_down_df, path, 'monocytes_down_atrx_30_12')
  
  message('Running Sporadic vs Control...')
  sporadic_ctrl_up <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupIV <= 0.05, 
                  logFC_groupIV > 0) %>% 
    pull(gene) %>% Signac::StringToGRanges() 
  sporadic_ctrl_down <- 
    big_data_frame %>% 
    dplyr::filter(cell_type == cell_type,
                  q_groupIV <= 0.05,
                  logFC_groupIV < 0) %>% 
    pull(gene) %>% Signac::StringToGRanges() 
  sporadic_ctrl_up_promoters <- .get_promoter_peaks(sporadic_ctrl_up)
  sporadic_ctrl_up_df <- .create_bed_df(sporadic_ctrl_up, sporadic_ctrl_up_promoters)
  .save_bed_files(sporadic_ctrl_up_df, path, 'monocytes_up_sporadic_30_12')
  sporadic_ctrl_down_promoters <- .get_promoter_peaks(sporadic_ctrl_down)
  sporadic_ctrl_down_df <- .create_bed_df(sporadic_ctrl_down, sporadic_ctrl_down_promoters)
  .save_bed_files(sporadic_ctrl_down_df, path, 'monocytes_down_sporadic_30_12')
}


beds_from_peak_df(
  big_data, 
  'Cluster_Monocytes',
  '/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/group_microenv')

# findMotifsGenome.pl \
#   /Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/group_microenv/monocytes_mycn_df.bed \
#   hg38 \
#   /Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/group_microenv/monocytes_mycn \
#   -size 100 \
#   -mask


# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster_Monocytes',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->monocytes_low_acc2

filter_peaks(
  atrx_vs_mycn,
  "gene_name" ,
  "logFC_groupIII")


big_data %>% 
  filter(cell_type == 'Cluster:Monocytes',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->monocytes_high_acc2

monocytes_high_acc_go2 <- 
  enrichr(monocytes_high_acc2, 
          c('MSigDB_Hallmark_2020'))
monocytes_high_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'High'
monocytes_low_acc_go2 <- 
  enrichr(monocytes_low_acc2, 
          c('MSigDB_Hallmark_2020'))
monocytes_low_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'Low'

monocytes_db2 <- rbind(
  monocytes_low_acc_go2$MSigDB_Hallmark_2020,
  monocytes_high_acc_go2$MSigDB_Hallmark_2020) %>%
  mutate(cellType = 'Cluster:Monocytes') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# B-cells
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:B-cells',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->bcells_low_acc2

big_data %>% 
  filter(cell_type == 'Cluster:B-cells',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->bcells_high_acc2

bcells_high_acc_go2 <- 
  enrichr(bcells_high_acc2, 
          c('MSigDB_Hallmark_2020'))
bcells_high_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'High'
bcells_low_acc_go2 <- 
  enrichr(bcells_low_acc2, 
          c('MSigDB_Hallmark_2020'))
bcells_low_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'Low'

bcells_low_acc_go2$MSigDB_Hallmark_2020$cellType <- 'Cluster:B-cells'
bcell_db2 <- bcells_low_acc_go2$MSigDB_Hallmark_2020 %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# T-cells
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:T-cells',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->tcells_low_acc2

big_data %>% 
  filter(cell_type == 'Cluster:T-cells',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->tcells_high_acc2

tcells_high_acc_go2 <- 
  enrichr(tcells_high_acc2, 
          c('MSigDB_Hallmark_2020'))
tcells_high_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'High'
tcells_low_acc_go2 <- 
  enrichr(tcells_low_acc2, 
          c('MSigDB_Hallmark_2020'))
tcells_low_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'Low'

#tcells_low_acc_go2$MSigDB_Hallmark_2020$cellType <- 'Cluster:T-cells'
tcell_db2 <- rbind(tcells_low_acc_go2$MSigDB_Hallmark_2020,
                   tcells_high_acc_go2$MSigDB_Hallmark_2020) %>% 
  mutate(cellType = 'Cluster:T-cells') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')


# Cluster::22
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:22',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster22_low_acc2

big_data %>% 
  filter(cell_type == 'Cluster:22',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster22_high_acc2

cluster22_low_acc_go2 <- 
  enrichr(cluster22_low_acc2, 
          c('MSigDB_Hallmark_2020'))
cluster22_low_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'Low'
cluster22_low_acc_go2$MSigDB_Hallmark_2020$cellType <- 'Cluster:22' 
cluster22_db2 <- cluster22_low_acc_go2$MSigDB_Hallmark_2020 %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# Cluster:32
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:32',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster32_low_acc2

big_data %>% 
  filter(cell_type == 'Cluster:32',
         q_groupII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupII < 0 ~'Low',
                                   logFC_groupII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster32_high_acc2

cluster32_high_acc_go2 <- 
  enrichr(cluster32_high_acc2, 
          c('MSigDB_Hallmark_2020'))
cluster32_high_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'High'
cluster32_low_acc_go2 <- 
  enrichr(cluster32_low_acc2, 
          c('MSigDB_Hallmark_2020'))
cluster32_low_acc_go2$MSigDB_Hallmark_2020$Accessibility <- 'Low'
#cluster32_low_acc_go2$MSigDB_Hallmark_2020$cellType <- 'Cluster:32' 

cluster32_db2 <- rbind(
  cluster32_low_acc_go$MSigDB_Hallmark_2020,
  cluster32_high_acc_go$MSigDB_Hallmark_2020) %>% 
  mutate(cellType = 'Cluster:32') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

mycn_vs_ctrl <- do.call(rbind, list(monocytes_db2,
                                    bcell_db2, 
                                    tcell_db2, 
                                    cluster22_db2, 
                                    cluster32_db2 )) %>% 
  ggplot(aes(y = Term, x = Accessibility)) +
  geom_point(
    aes(size = log(Odds.Ratio), color= Accessibility)) + 
  labs(color = 'Chromatin accessibility',
       title = 'MYCN vs Control', y = '', x = '', 
       caption = 'MSigDB Hallmark 2020') +
  hrbrthemes::theme_ipsum() +
  #ggpubr::theme_pubclean() +
  facet_wrap(~ cellType, nrow = 1) +
  scale_color_manual(values = c('High' = "steelblue1",'Low' = "firebrick")) +
  scale_fill_manual(values = c('High' = "steelblue2",'Low' = "firebrick")) 

# Section_2 ---------------------------------------------------------------

# Monocytes
# atrx vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:Monocytes',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->monocytes_low_acc3

big_data %>% 
  filter(cell_type == 'Cluster:Monocytes',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->monocytes_high_acc3

monocytes_high_acc_go3 <- 
  enrichr(monocytes_high_acc3, 
          c('MSigDB_Hallmark_2020'))
monocytes_high_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'High'
monocytes_low_acc_go3 <- 
  enrichr(monocytes_low_acc3, 
          c('MSigDB_Hallmark_2020'))
monocytes_low_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'Low'
monocytes_db3 <- rbind(
  monocytes_low_acc_go3$MSigDB_Hallmark_2020,
  monocytes_high_acc_go3$MSigDB_Hallmark_2020) %>%
  mutate(cellType = 'Cluster:Monocytes') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# B-cells
# atrx vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:B-cells',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->bcells_low_acc3

big_data %>% 
  filter(cell_type == 'Cluster:B-cells',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->bcells_high_acc3
bcells_high_acc_go3 <- 
  enrichr(bcells_high_acc3, 
          c('MSigDB_Hallmark_2020'))
bcells_high_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'High'
bcells_low_acc_go3 <- 
  enrichr(bcells_low_acc3, 
          c('MSigDB_Hallmark_2020'))
bcells_low_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'Low'

bcells_low_acc_go3$MSigDB_Hallmark_2020$cellType <- 'Cluster:B-cells'
bcell_db3 <- bcells_low_acc_go3$MSigDB_Hallmark_2020 %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# T-cells
# atrx vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:T-cells',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->tcells_low_acc3

big_data %>% 
  filter(cell_type == 'Cluster:T-cells',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->tcells_high_acc3

tcells_high_acc_go3 <- 
  enrichr(tcells_high_acc3, 
          c('MSigDB_Hallmark_2020'))
tcells_high_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'High'
tcells_low_acc_go3 <- 
  enrichr(tcells_low_acc3, 
          c('MSigDB_Hallmark_2020'))
tcells_low_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'Low'

tcell_db3 <- rbind(tcells_low_acc_go3$MSigDB_Hallmark_2020,
                   tcells_high_acc_go3$MSigDB_Hallmark_2020) %>% 
  mutate(cellType = 'Cluster:T-cells') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')


# Cluster::22
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:22',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIII) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster22_low_acc3

big_data %>% 
  filter(cell_type == 'Cluster:22',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIII)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster22_high_acc3
cluster22_low_acc_go3 <- 
  enrichr(cluster22_low_acc3, 
          c('MSigDB_Hallmark_2020'))
cluster22_low_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'Low'
cluster22_low_acc_go3$MSigDB_Hallmark_2020$cellType <- 'Cluster:22' 
cluster22_db3 <- cluster22_low_acc_go3$MSigDB_Hallmark_2020 %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# Cluster:32
# atrx vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:32',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIII) %>% 
  head(100) %>% 
  pull(gene_name) %>% 
  unique() ->cluster32_low_acc3

big_data %>% 
  filter(cell_type == 'Cluster:32',
         q_groupIII <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIII < 0 ~'Low',
                                   logFC_groupIII > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIII)) %>% 
  head(100) %>% 
  pull(gene_name) %>% 
  unique() ->cluster32_high_acc3

cluster32_high_acc_go3 <- 
  enrichr(cluster32_high_acc3, 
          c('MSigDB_Hallmark_2020'))
cluster32_high_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'High'
cluster32_low_acc_go3 <- 
  enrichr(cluster32_low_acc3, 
          c('MSigDB_Hallmark_2020'))
cluster32_low_acc_go3$MSigDB_Hallmark_2020$Accessibility <- 'Low'
cluster32_low_acc_go3$MSigDB_Hallmark_2020$cellType <- 'Cluster:32' 

cluster32_db3 <- rbind(
  cluster32_low_acc_go3$MSigDB_Hallmark_2020,
  cluster32_high_acc_go3$MSigDB_Hallmark_2020) %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')
#, cluster22_db, cluster32_db 
atrx_vs_ctrl <- do.call(rbind, list(monocytes_db,bcell_db, tcell_db)) %>% 
  ggplot(aes(y = Term, x = Accessibility)) +
  geom_point(
    aes(size = log(Odds.Ratio), color= Accessibility)) + 
  labs(color = 'Chromatin accessibility',
       title = 'ATRX vs Control', y = '', x = '', 
       caption = 'MSigDB Hallmark 2020') +
  hrbrthemes::theme_ipsum() +
  facet_wrap(~ cellType, nrow = 1) +
  #ggpubr::theme_pubclean() +
  scale_color_manual(values = c('High' = "steelblue1",'Low' = "firebrick")) +
  scale_fill_manual(values = c('High' = "steelblue1",'Low' = "firebrick")) 

# Section_3 ---------------------------------------------------------------

# Monocytes
# sporadic vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:Monocytes',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIV) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->monocytes_low_acc4

big_data %>% 
  filter(cell_type == 'Cluster:Monocytes',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIV)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->monocytes_high_acc4

monocytes_high_acc_go4 <- 
  enrichr(monocytes_high_acc4, 
          c('MSigDB_Hallmark_2020'))
monocytes_high_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'High'
monocytes_low_acc_go4 <- 
  enrichr(monocytes_low_acc4, 
          c('MSigDB_Hallmark_2020'))
monocytes_low_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'Low'

monocytes_db4 <- rbind(
  monocytes_low_acc_go4$MSigDB_Hallmark_2020,
  monocytes_high_acc_go4$MSigDB_Hallmark_2020) %>%
  mutate(cellType = 'Cluster:Monocytes') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')


# T-cells
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:T-cells',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIV) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->tcells_low_acc4

big_data %>% 
  filter(cell_type == 'Cluster:T-cells',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIV)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->tcells_high_acc4

tcells_high_acc_go4 <- 
  enrichr(tcells_high_acc4, 
          c('MSigDB_Hallmark_2020'))
tcells_high_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'High'
tcells_low_acc_go4 <- 
  enrichr(tcells_low_acc4, 
          c('MSigDB_Hallmark_2020'))
tcells_low_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'Low'

tcell_db4 <- rbind(tcells_low_acc_go4$MSigDB_Hallmark_2020, 
                   tcells_high_acc_go4$MSigDB_Hallmark_2020) %>% 
  mutate(cellType = 'Cluster:T-cells') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')


# Cluster::11
# mycn vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:11',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIV) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster11_low_acc4

big_data %>% 
  filter(cell_type == 'Cluster:11',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIV)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster11_high_acc4

cluster11_low_acc_go4 <- 
  enrichr(cluster11_low_acc4, 
          c('MSigDB_Hallmark_2020'))
cluster11_low_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'Low'

cluster11_high_acc_go4 <- 
  enrichr(cluster11_high_acc4, 
          c('MSigDB_Hallmark_2020'))
cluster11_high_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'High'

#cluster11_low_acc_go$MSigDB_Hallmark_2020$cellType <- 'Cluster:22' 
cluster11_db4 <- cluster11_low_acc_go4$MSigDB_Hallmark_2020 %>% 
  mutate(cellType = 'Cluster:11') %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')

# Cluster::22
# sporadic vs ctrl enrichment
big_data %>% 
  filter(cell_type == 'Cluster:22',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'Low') %>% 
  arrange( logFC_groupIV) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster22_low_acc4

big_data %>% 
  filter(cell_type == 'Cluster:22',
         q_groupIV <= 0.05) %>% 
  mutate(accessibility = case_when(logFC_groupIV < 0 ~'Low',
                                   logFC_groupIV > 0 ~'High')) %>% 
  filter(accessibility == 'High') %>% 
  arrange(desc(logFC_groupIV)) %>% 
  head(50) %>% 
  pull(gene_name) %>% 
  unique() ->cluster22_high_acc4
cluster22_low_acc_go4 <- 
  enrichr(cluster22_low_acc4, 
          c('MSigDB_Hallmark_2020'))
cluster22_low_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'Low'

cluster22_high_acc_go4 <- 
  enrichr(cluster22_high_acc4, 
          c('MSigDB_Hallmark_2020'))
cluster22_high_acc_go4$MSigDB_Hallmark_2020$Accessibility <- 'High'

cluster22_db4 <- rbind(cluster22_low_acc_go4$MSigDB_Hallmark_2020,
                       cluster22_high_acc_go4$MSigDB_Hallmark_2020) %>% 
  mutate(cellType = 'Cluster:22' ) %>% 
  select(Term, Odds.Ratio, Accessibility, cellType) %>% 
  pivot_longer(!c(Term, Accessibility, cellType), 
               values_to = 'Odds.Ratio')


sporadic_vs_ctrl <- do.call(rbind, 
                            list(monocytes_db4, 
                                 tcell_db4, 
                                 cluster11_db4, 
                                 cluster22_db4 )) %>% 
  ggplot(aes(y = Term, x =  Accessibility)) +
  geom_point(
    aes(size = log(Odds.Ratio), color= Accessibility)) + 
  labs(color = 'Chromatin accessibility',
       title = 'Sporadic vs Control', y = '', x = '', 
       caption = 'MSigDB Hallmark 2020') +
  hrbrthemes::theme_ipsum() +
  #ggpubr::theme_pubclean() +
  facet_wrap(~ cellType, nrow = 1) +
  #facet_grid(Accessibility ~ ., scales = "free") +
  scale_color_manual(values = c('High' = "steelblue1",'Low' = "firebrick")) +
  scale_fill_manual(values = c('High' = "steelblue1",'Low' = "firebrick")) 


# One to One mapped data --------------------------------------------------

filter_peaks_microenv <- 
  function(datafmr, gene_col, logfc_col) {

    mybiglist <- list()
    special_count = c()
    for (row in 1:nrow(datafmr)) {
      gene_name <- 
        datafmr[row, gene_col]
      logFC  <- 
        datafmr[row, logfc_col]
      
      if (exists(gene_name, where=mybiglist)) {
        
        if  (( mybiglist[[gene_name]] >0  & logFC >0) | (
          mybiglist[[gene_name]] < 0  & logFC  < 0)) {
          if (abs(mybiglist[[gene_name]]) < abs(logFC)) {
            message(
              "Current value larger than previous: ",
              paste(gene_name, "\t",mybiglist[[gene_name]], logFC))
            mybiglist[[gene_name]] <- logFC
          } 
          else if (mybiglist[[gene_name]] > logFC) {
            message(
              "Current value smaller than previous",
              paste(gene_name, "\t",mybiglist[[gene_name]], logFC))
            next
          }
        } else if (( mybiglist[[gene_name]] > 0  & logFC < 0) | (
          mybiglist[[gene_name]] < 0  & logFC  > 0)){
          
          message("Special case: ", 
                  gene_name, "\t", 
                  mybiglist[[gene_name]], "\t",
                  logFC )
          special_count = c(special_count,gene_name)
        }
        
      } else {
        mybiglist[[gene_name]] <- logFC
      }
    }
    
    # # mybiglist
    mappings <-
      mapIds(org.Hs.eg.db,
             names(mybiglist),
             'ENTREZID',
             'SYMBOL')
    mapping_df <-
      enframe(mappings) %>%
      unnest()
    # 
    special_mappings <-
      mapIds(org.Hs.eg.db,
             special_count,
             'ENTREZID',
             'SYMBOL')
      
    #
    mydf <-
      enframe(mybiglist) %>%
      unnest() %>%
      left_join(mapping_df, by='name') %>%
      mutate(gene_name = name,
             logFC = value.x,
             entrez = value.y) %>%
      dplyr::select(gene_name, logFC, entrez) %>%
      dplyr::filter(!is.na(entrez))
    
    return(mydf)
  }

# Group II comparisons
comparison2 <- list()

for (i in names(table(big_data$cell_type))) {
   print(i)
  tmp_grpII_df <- big_data %>% 
    dplyr::filter(cell_type == i) %>% 
    dplyr::filter(p_groupII < 0.05 & abs(logFC_groupII) >1 )
  filtered_grpII_df <- filter_peaks_microenv(
    tmp_grpII_df, 'gene_name', 'logFC_groupII' )
  filtered_grpII_df$cell_type <- i
  #
  comparison2[[i]] <- filtered_grpII_df
  #
  
}

kegg2 <- list()
msigdb2 <- list()
reactome2 <- list()
namesl2 <- names(comparison2)

for(i in seq_along(comparison2)) {
  name <- namesl2[[i]]
  print(name)
  #
  msigdb_results2 <-     run_msigdb(comparison2[[name]])
  if (dim(msigdb_results2)[1] != 0) {
    msigdb_results2@result$cell_type <- name
  } else { next }
  msigdb2[[name]] <- msigdb_results2@result %>% 
    arrange(pvalue)  %>% 
    head(20)
  #
  reactomedb_results2 <- run_reactome(comparison2[[name]])
  if (dim(reactomedb_results2)[1] != 0) {
    reactomedb_results2@result$cell_type <- name
  } else { next }
  reactome2[[name]] <- reactomedb_results2@result %>% 
    arrange(pvalue)  %>% 
    head(20)
  #
  kegg_results2 <-       run_kegg(comparison2[[name]])
  if (dim(kegg_results2)[1] != 0) {
    kegg_results2@result$cell_type <- name
  } else { next }
  kegg2[[name]] <- kegg_results2@result %>% 
    arrange(pvalue)  %>% 
    head(20)
}

kegg2_data = do.call(rbind, kegg2)
msigdb2_data <- do.call(rbind, msigdb2)
reactome2_data <- do.call(rbind, reactome2)

msigdb2_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'MYCN vs CTRL', caption= 'MsigDB Hallmark DB')

kegg2_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'MYCN vs CTRL', caption= 'KEGG DB')

reactome2_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'MYCN vs CTRL', caption= 'ReactomeDB')

# Group III comparisons
#  names(table(big_data$cell_type))
msigdb3[['Cluster_Monocytes']] %>% names()
comparison3 <- list()

for (i in names(table(big_data$cell_type))) {
  tmp_grpIII_df <- big_data %>% 
    dplyr::filter(cell_type == i) %>% 
    dplyr::filter(p_groupIII < 0.05 & abs(logFC_groupIII) >1 )
  filtered_grpIII_df <- filter_peaks_microenv(
    tmp_grpIII_df, 'gene_name', 'logFC_groupIII' )
  filtered_grpIII_df$cell_type <- i
  #
  comparison3[[i]] <- filtered_grpIII_df

}

kegg3 <- list()
msigdb3 <- list()
reactome3 <- list()
namesl3 <- names(comparison3)

for(i in seq_along(comparison3)) {
  name <- namesl3[[i]]
  print(name)
  #
  msigdb_results3 <-     run_msigdb(comparison3[[i]])
  if (dim(msigdb_results3)[1] != 0) {
    msigdb_results3@result$cell_type <- name
  } else { next }
  msigdb3[[name]] <- msigdb_results3@result %>% 
    arrange(pvalue)  %>% 
    head(20)
  #
  reactomedb_results3 <- run_reactome(comparison3[[i]])
  if (dim(reactomedb_results3)[1] != 0) {
    reactomedb_results3@result$cell_type <- name
  } else { next }
  reactome3[[name]] <- reactomedb_results3@result %>% 
    arrange(pvalue)  %>% 
    head(20)
  #
  kegg_results3 <-       run_kegg(comparison3[[i]])
  if (dim(kegg_results3)[1] != 0) {
    kegg_results3@result$cell_type <- name
  } else { next }
  kegg3[[name]] <- kegg_results@result %>% 
    arrange(pvalue)  %>% 
    head(20)

}
kegg3_data = do.call(rbind, kegg3)
msigdb3_data <- do.call(rbind, msigdb3)
reactome3_data <- do.call(rbind, reactome3)

msigdb3_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Clusters", y = "Terms",
       title = 'ATRX vs CTRL', caption= 'MsigDB Hallmark DB')

kegg3_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'ATRX vs CTRL', caption= 'KEGG DB')

reactome3_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'ATRX vs CTRL', caption= 'ReactomeDB')

# Group IV comparisons
comparison4 <- list()

for (i in names(table(big_data$cell_type))) {
  
  tmp_grpIV_df <- big_data %>% 
    dplyr::filter(cell_type == i) %>% 
    dplyr::filter(p_groupIV < 0.05 & abs(logFC_groupIV) >1 )
  filtered_grpIV_df <- filter_peaks_microenv(
    tmp_grpIV_df, 'gene_name', 'logFC_groupIV' )
  filtered_grpIV_df$cell_type <- i
  #
  comparison4[[i]] <- filtered_grpIV_df
  #
  print(paste('Sporadic vs CTRL', i))

}
kegg4 <- list()
msigdb4 <- list()
reactome4 <- list()
namesl4 <- names(comparison4)

for(i in seq_along(comparison4)) {
  name <- namesl4[[i]]
  print(name)
  #
  msigdb_results4 <-     run_msigdb(comparison4[[i]])
   if (dim(msigdb_results4)[1] != 0) {
     msigdb_results4@result$cell_type <- name
   } else { next }
  msigdb4[[name]] <- msigdb_results4@result %>% 
    arrange(pvalue)  %>% 
    head(20)
  #
  reactomedb_results4 <- run_reactome(comparison4[[i]])
  if (dim(reactomedb_results4)[1] != 0) {
    reactomedb_results4@result$cell_type <- name
  } else { next }
  reactome4[[name]] <- reactomedb_results4@result %>% 
    arrange(pvalue)  %>% 
    head(20)
  #
  kegg_results4 <-       run_kegg(comparison4[[i]])
  if (dim(kegg_results4)[1] != 0) {
    kegg_results4@result$cell_type <- name
  } else { next }
  kegg4[[name]] <- kegg_results4@result %>% 
    arrange(pvalue)  %>% 
    head(20)
}


kegg4_data = do.call(rbind, kegg4)
msigdb4_data <- do.call(rbind, msigdb4)
reactome4_data <- do.call(rbind, reactome4)

msigdb4_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'Sporadic vs CTRL', caption= 'MsigDB Hallmark DB') +
  theme(axis.text.x = element_text(angle = 90))

kegg4_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'Sporadic vs CTRL', caption= 'KEGG DB') +
  theme(axis.text.x = element_text(angle = 90))

reactome4_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  hrbrthemes::theme_ipsum() +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'Sporadic vs CTRL', caption= 'ReactomeDB') +
  theme(axis.text.x = element_text(angle = 90))



atrx_vs_mycn <- 
  results_with_closest_gene %>% 
  mutate(logFC_groupIII = logFC_groupIII/log(2)) %>% 
  dplyr::filter(p_groupIII < 0.05 & abs(logFC_groupIII) >1 ) %>%
  dplyr::select(gene_name, logFC_groupIII) 


mycn_df_monocytes_test <- filter_peaks_microenv(
  big_data,
  "gene_name" ,
  "logFC_groupIII", 
  'Cluster_Monocytes')

mycn_df_monocytes_msigdb <- run_msigdb(mycn_df_monocytes_test)
mycn_df_monocytes_dotplot <- dotplot(mycn_df_monocytes_msigdb, 
                               showCategory = 20, 
                               title = "Monocytes ATRX vs CTRL MsigDB Hallmark Pathways", 
                               split=".sign") + facet_grid(.~.sign)

mycn_df_monocytes_msigdb@result %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = NES, y = reorder(Description, -NES), fill = FoldChange)) +
  geom_col(position = "identity") +
  scale_fill_manual(values = c('negative' = "springgreen4", 'positive' = "navy"))+
  hrbrthemes::theme_ipsum() +
  #cowplot::theme_cowplot() +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'ATRX vs CTRL Monocytes', subtitle= 'MsigDB Hallmark DB')
#geom_point(aes(size = p.adjust,  color = NES))

pal <- wes_palette("Zissou1", 100, type = "continuous")

kegg3_data %>% 
  mutate(FoldChange = ifelse(NES > 0, 'positive', 'negative')) %>% 
  ggplot(aes(x = cell_type, y = Description)) +
  geom_point(aes(size = -log10(p.adjust),  color = NES)) +
  scale_colour_gradient2(high = 'magenta', low = 'deepskyblue', mid = 'grey40', midpoint= 0 ) +
  ggpubr::theme_pubclean() +
  labs(x = "Normalized Enrichment Score", y = "Terms",
       title = 'ATRX vs CTRL', caption= 'MsigDB Hallmark DB')


# Save Tables -------------------------------------------------------------
# peak, logFC, gene, annotation, distanceToTSS, p_value, q_value, comparison

# Comparison MycN vs Control
tII <- 
  big_data %>% 
  filter(cell_type == "Cluster_T-cells" & significanceII !="No change")
tIII <- 
  big_data %>% 
  filter(cell_type == "Cluster_T-cells" & significanceIII !="No change")
tIV <- 
  big_data %>% 
  filter(cell_type == "Cluster_T-cells" & significanceIV !="No change")

tIIdb <- tII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(tII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
         gene = SYMBOL,
         logFC = logFC_groupII,
         p_value = p_groupII,
         q_value = q_groupII,
         comparison  = "T_MycN_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
         distanceToTSS, p_value, q_value, 
         comparison)
tIIIdb <- tIII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(tIII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIII,
                p_value = p_groupIII,
                q_value = q_groupIII,
                comparison  = "T_Atrx_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
tIVdb <- tIV %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(tIV, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIV,
                p_value = p_groupIV,
                q_value = q_groupIV,
                comparison  = "T_Sporadic_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
t_db <- do.call(rbind, list(tIIdb,tIIIdb,tIVdb))

bII <- 
  big_data %>% 
  filter(cell_type == "Cluster_B-cells" & significanceII !="No change")
bIII <- 
  big_data %>% 
  filter(cell_type == "Cluster_B-cells" & significanceIII !="No change")
bIV <- 
  big_data %>% 
  filter(cell_type == "Cluster_B-cells" & significanceIV !="No change")

bIIdb <- bII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(bII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupII,
                p_value = p_groupII,
                q_value = q_groupII,
                comparison  = "B_MycN_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
bIIIdb <- bIII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(bIII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIII,
                p_value = p_groupIII,
                q_value = q_groupIII,
                comparison  = "B_Atrx_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
bIVdb <- bIV %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(bIV, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIV,
                p_value = p_groupIV,
                q_value = q_groupIV,
                comparison  = "B_Sporadic_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
b_db <- do.call(rbind, list(bIIdb,bIIIdb,bIVdb))

mbII <- 
  big_data %>% 
  filter(cell_type == "Cluster_B_cells_memory_34" & significanceII !="No change")
mbIII <- 
  big_data %>% 
  filter(cell_type == "Cluster_B_cells_memory_34" & significanceIII !="No change")
mbIV <- 
  big_data %>% 
  filter(cell_type == "Cluster_B_cells_memory_34" & significanceIV !="No change")


mbIIdb <- mbII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(mbII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupII,
                p_value = p_groupII,
                q_value = q_groupII,
                comparison  = "MB_MycN_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
mbIIIdb <- mbIII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(mbIII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIII,
                p_value = p_groupIII,
                q_value = q_groupIII,
                comparison  = "MB_Atrx_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)

mbIVdb <- mbIV %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(mbIV, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIV,
                p_value = p_groupIV,
                q_value = q_groupIV,
                comparison  = "B_Sporadic_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
mb_db <- mbIIIdb

eryII <-
  big_data %>% 
  filter(cell_type == "Cluster_22" & significanceII !="No change") 
# eryIII <-
#   big_data %>% 
#   filter(cell_type == "Cluster_22" & significanceIII !="No change")
eryIV <-
  big_data %>% 
  filter(cell_type == "Cluster_22" & significanceIV !="No change")

eryIIdb <- eryII %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(eryII, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupII,
                p_value = p_groupII,
                q_value = q_groupII,
                comparison  = "E_MycN_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
eryIVdb <- eryIV %>% 
  pull(gene) %>% 
  .get_chipseeker_annot(.) %>% 
  left_join(eryIV, by = 'gene')  %>%
  dplyr::mutate(peak = gene,
                gene = SYMBOL,
                logFC = logFC_groupIV,
                p_value = p_groupIV,
                q_value = q_groupIV,
                comparison  = "E_Sporadic_Ctrl") %>% 
  dplyr::select(peak, logFC, gene, annotation, 
                distanceToTSS, p_value, q_value, 
                comparison)
e_db <- do.call(rbind, list(eryIIdb, eryIVdb))

enrich_microenv_db <-
  do.call(rbind, list(t_db, b_db, mb_db, e_db))
write.csv(enrich_microenv_db,
          file = "/media/AGFORTELNY/people/rohit/projects/enrich_microenv.csv", 
          row.names = FALSE, 
          col.names = TRUE,
          quote = TRUE)


