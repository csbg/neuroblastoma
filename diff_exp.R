library(Seurat)
library(edgeR)
library(Matrix.utils)
library(tidyverse)


nb <- readRDS("data_generated/all_datasets_current/nb_only_RNA.rds")

meta_data <-
  nb@meta.data %>%
  as_tibble(rownames = "cell") %>% 
  left_join(
    read_csv("data_raw/metadata/sample_groups.csv", comment = "#"),
    by = "sample"
  )

counts_agg <- 
  nb$RNA@counts %>% 
  t() %>% 
  aggregate.Matrix(
    meta_data %>% select(group, cluster = integrated_snn_res.0.5)
  ) %>% 
  t()

targets <-
  meta_data %>%
  distinct(group, cluster = integrated_snn_res.0.5) %>% 
  arrange(group, cluster) %>%
  mutate(name = str_c(group, cluster, sep = "_")) %>% 
  column_to_rownames("name")
design <- model.matrix(~group + cluster, targets)

dge <- DGEList(counts_agg, group = colnames(counts_agg))
dge


keep <- filterByExpr(dge)
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge$samples

plotMDS(y, col = rep(1:4, 19))


dge <- estimateDisp(dge, design)
plotBCV(dge)

fit <- glmQLFit(dge, design = design)
qlf <- glmQLFTest(fit, coef = 8)
# qlf <- glmQLFTest(fit, contrast = makeContrasts(levels = design))
topTags(qlf)

FindMarkers()
setactive ()

# Design matrix for patients and groups -----------------------------------

Patient <- factor(rep(1:9, each = 2))
Disease <- factor(rep(c("Healthy", "Disease1", "Disease2"), each = 6), levels=c("Healthy","Disease1","Disease2"))
Treatment <- factor(rep(c("None","Hormone"), 9), levels=c("None","Hormone"))

length(Treatment)

design <- model.matrix(~Patient)

Healthy.Hormone <- Disease == "Healthy" & Treatment == "Hormone"
Disease1.Hormone <- Disease =="Disease1" & Treatment == "Hormone"
Disease2.Hormone <- Disease =="Disease2" & Treatment == "Hormone"
design <- cbind(design, Healthy.Hormone, Disease1.Hormone, Disease2.Hormone)



targets <- 
  meta_data %>%
  distinct(disease = group, patient = sample, cluster = integrated_snn_res.0.5)

Patient <- factor(targets$patient)
Disease <- factor(targets$disease, levels=c("I", "II", "III", "IV"))
Cluster <- factor(targets$cluster)
design <- model.matrix(~Patient)
Healthy.Hormone <- Disease == "I" & Cluster == "1"