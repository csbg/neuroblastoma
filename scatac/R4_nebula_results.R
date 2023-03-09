pacman::p_load(Seurat, Signac, SingleCellExperiment, 
               tidyverse, nebula, scater)
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
 

# Monocytes ---------------------------------------------------------------

mono_iv <- big_data %>%
  mutate(q_groupII = qvalue::qvalue(big_data$p_groupII)$qvalues) %>% 
  mutate(q_groupIII = qvalue::qvalue(big_data$p_groupIII)$qvalues) %>%
  mutate(q_groupIV = qvalue::qvalue(big_data$p_groupIV)$qvalues) %>% 
  mutate(logFC_groupII = logFC_groupII/log(2),
         logFC_groupIII = logFC_groupIII/log(2),
         logFC_groupIV = logFC_groupIV/log(2)) %>% 
  mutate(significance = case_when(
    (logFC_groupIV < -1 & p_groupIV < 0.05 )~ 'Down',
    (logFC_groupIV > 1 & p_groupIV < 0.05 )~ 'Up',
    TRUE ~ 'No change')) %>% 
  mutate(delabel = case_when(
    significance != 'No change' ~ gene_name)) %>%
  filter(cell_type == "Cluster_Monocytes") %>% 
  ggplot(aes(x = logFC_groupIV, 
             y = -log10(p_groupIV), color = significance)) +
  geom_point(size = 0.5) + #geom_text(aes(label=delabel), size = 2) +
  ylim(0, 10) +  xlim(-10, 10) + 
  labs(x = 'LogFC_Sporadic_Ctrl', y = '-Log10(p_value)') +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red",linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") + 
  ggpubr::theme_pubclean() +
  theme(legend.position="none")

mono_iii <- big_data %>%
  mutate(q_groupII = qvalue::qvalue(big_data$p_groupII)$qvalues) %>% 
  mutate(q_groupIII = qvalue::qvalue(big_data$p_groupIII)$qvalues) %>%
  mutate(q_groupIV = qvalue::qvalue(big_data$p_groupIV)$qvalues) %>% 
  mutate(logFC_groupII = logFC_groupII/log(2),
         logFC_groupIII = logFC_groupIII/log(2),
         logFC_groupIV = logFC_groupIV/log(2)) %>% 
  mutate(significance = case_when(
    (logFC_groupIII < -1 & p_groupIII < 0.05 )~ 'Down',
    (logFC_groupIII > 1 & p_groupIII < 0.05 )~ 'Up',
    TRUE ~ 'No change')) %>% 
  mutate(delabel = case_when(
    significance != 'No change' ~ gene_name)) %>%
  filter(cell_type == "Cluster_Monocytes") %>% 
  ggplot(aes(x = logFC_groupIII, 
             y = -log10(p_groupIII), color = significance)) +
  geom_point(size = 0.5) + #geom_text(aes(label=delabel), size = 2) + 
  ylim(0, 10) + xlim(-10, 10) + 
  labs( x = 'LogFC_ATRX_Ctrl', y = '-Log10(p_value)') +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") + 
  ggpubr::theme_pubclean() +
  theme(legend.position="none")


mono_ii <- big_data %>%
  mutate(q_groupII = qvalue::qvalue(big_data$p_groupII)$qvalues) %>% 
  mutate(q_groupIII = qvalue::qvalue(big_data$p_groupIII)$qvalues) %>%
  mutate(q_groupIV = qvalue::qvalue(big_data$p_groupIV)$qvalues) %>% 
  mutate(logFC_groupII = logFC_groupII/log(2),
         logFC_groupIII = logFC_groupIII/log(2),
         logFC_groupIV = logFC_groupIV/log(2)) %>% 
  mutate(significance = case_when(
    (logFC_groupII < -1 & q_groupII < 0.05 )~ 'Down',
    (logFC_groupII > 1 & q_groupII < 0.05 )~ 'Up',
    TRUE ~ 'No change')) %>% 
  mutate(delabel = case_when(
    significance != 'No change' ~ gene_name)) %>%
  filter(cell_type == "Cluster_Monocytes") %>% 
  ggplot(aes(x = logFC_groupII, 
             y = -log10(q_groupII), color = significance)) +
  geom_point(size = 0.5) + #geom_text(aes(label=delabel), size = 2) + 
  ylim(0, 10) + xlim(-10, 10) + 
  labs( x = 'LogFC_MYCN_Ctrl', y = '-Log10(q_value)') +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red", linetype = "longdash") +
  geom_hline(yintercept=-log10(0.05), col="red", linetype = "longdash") + 
  ggpubr::theme_pubclean() +
  theme(legend.position="none")

mono_ii + mono_iii + mono_iv +
  patchwork::plot_annotation(title = 'Monocytes differential accessible regions')



# scratch -----------------------------------------------------------------


logFC_(Intercept)

summary(result$summary$`logFC_(Intercept)`)
pvalue <- result$summary$`p_(Intercept)`

summary(result$summary$`p_(Intercept)`)
pvalues <- hedenfalk$p

qobj <- qvalue(p = pvalue)

# B-cells -----------------------------------------------------------------

results_with_closest_gene <- 
  read.table('/Volumes/GFS_MBIO_AGFORTELNY/people/rohit/projects/clusterObj/B-cells.csv',
             header= TRUE, sep = ' ')
results_with_closest_gene %>%
  mutate(significance = case_when(
    (logFC_groupIV < -2 & p_groupIV < 0.05 )~ 'Down',
    (logFC_groupIV > 2 & p_groupIV < 0.05 )~ 'Up',
    TRUE ~ 'No change')) %>% 
  mutate(delabel = case_when(
    significance != 'No change' ~ gene_name)) %>%
  ggplot(aes(x = logFC_groupIV, 
             y = -log10(p_groupIV), color = significance)) +
  geom_point(size = 0.5) + #geom_text(aes(label=delabel), size = 2) + 
  labs(title = "B-cell",x = 'LogFC_Sporadic_MYCN', y = '-Log10(p_value)') +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggpubr::theme_pubclean() +
  theme(legend.title = element_blank()) 


results_with_closest_gene %>%
  mutate(significance = case_when(
    (logFC_groupIII < -2 & p_groupIII < 0.05 )~ 'Down',
    (logFC_groupIII > 2 & p_groupIII < 0.05 )~ 'Up',
    TRUE ~ 'No change')) %>% 
  mutate(delabel = case_when(
    significance != 'No change' ~ gene_name)) %>%
  ggplot(aes(x = logFC_groupIII, 
             y = -log10(p_groupIII), color = significance)) +
  geom_point(size = 0.5) + #geom_text(aes(label=delabel), size = 2) + 
  labs(title = "B-cell", x = 'LogFC_ATRX_MYCN', y = '-Log10(p_value)') +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggpubr::theme_pubclean() +
  theme(legend.title = element_blank()) 


results_with_closest_gene %>% 
  mutate(q_groupII = qvalue(results_with_closest_gene$p_groupII)$qvalues) %>% 
  mutate(q_groupIII = qvalue(results_with_closest_gene$p_groupIII)$qvalues) %>%
  mutate(q_groupIV = qvalue(results_with_closest_gene$p_groupIV)$qvalues) %>%
  
  results_with_closest_gene$q_groupII <- qvalue(results_with_closest_gene$p_groupII)$qvalues
results_with_closest_gene$q_groupIII <- qvalue(results_with_closest_gene$p_groupIII)$qvalues
results_with_closest_gene$q_groupIV <- qvalue(results_with_closest_gene$p_groupIV)$qvalues

results_with_closest_gene %>% 
  filter(q_groupIV <0.05) %>% View()

pvalues <- results_with_closest_gene$p_groupII
qobj <- qvalue(p = pvalues)
str(qobj)
summary(qobj$qvalues)  