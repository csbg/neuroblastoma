library(fs)
library(tidyverse)

fastq_files <-
  read_csv("metadata/sample_groups.csv", comment = "#") %>%
  mutate(raw_folder = str_glue("/mnt/agfortelny/PROJECTS/Neuroblastoma/data_raw/bsf/RAW/{bsf_id}_Transcriptome")) %>% 
  select(bsf_id, raw_folder) %>% 
  deframe() %>% 
  map_dfr(dir_info, .id = "bsf_id")

fastq_files %>% 
  summarise(size = sum(size))