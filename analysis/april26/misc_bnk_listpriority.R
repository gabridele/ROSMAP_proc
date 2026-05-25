library(dplyr)

bnk_qc <- read_csv("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/BNK_sub_ses_QC1805.csv")

priority_list <- read_csv("~/Desktop/priority_list_xcpd_dl.txt", col_names = FALSE) %>%
  select(sub_ses = X2) %>%
  mutate(priority = "priority")

bnk_qc <- bnk_qc %>%
  left_join(priority_list, by = "sub_ses")

write.csv(bnk_qc, "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/BNK_sub_ses_QC2205.csv", row.names = FALSE)
