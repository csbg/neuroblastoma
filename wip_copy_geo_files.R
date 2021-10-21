library(fs)
library(tidyverse)
source("styling.R")



# Create table with files to be copied ------------------------------------

file_table <-
  read_csv("metadata/sample_groups.csv", comment = "#") %>%
  arrange(group, sample) %>%
  mutate(
    patient =
      rename_patients(sample) %>%
      make.unique() %>%
      recode(M2 = "M2a", M2.1 = "M2b"),
    files = list(c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
  ) %>%
  unnest_longer(files) %>%
  mutate(
    source = str_glue("/mnt/agfortelny/DATA_Neuroblastoma/data_raw/COUNT/{bsf_id}_transcriptome/filtered_feature_bc_matrix/{files}"),
    target = str_glue("/mnt/agfortelny/PROJECTS/Neuroblastoma/data/geo_submission/{patient}_{files}"),
  ) %>%
  select(source, target)

file_table %>% write_delim("geo_files.txt", col_names = FALSE)



# Run bash script ---------------------------------------------------------

# while read -r source target; do
#   echo "copying $source"
#   cp "$source" "$target"
# done < geo_files.txt



# Check md5sums -----------------------------------------------------------

target_sums <- read_delim(
  "/mnt/agfortelny/PROJECTS/Neuroblastoma/data/geo_submission/md5sums.txt",
  delim = " ",
  trim_ws = TRUE,
  col_names = c("md5", "file")
)

source_sums <- read_csv(
  "/mnt/agfortelny/PROJECTS/Neuroblastoma/data/md5sums_bsf.csv",
  comment = "#"
)

file_table %>% 
  extract(source, into = "source", regex = "data_raw/(.+)") %>% 
  extract(target, into = "target", regex = "geo_submission/(.+)") %>% 
  left_join(source_sums, by = c(source = "file")) %>% 
  left_join(target_sums, by = c(target = "file")) %>% 
  filter(md5.x != md5.y)
  