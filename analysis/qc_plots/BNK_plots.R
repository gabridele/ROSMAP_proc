library(dplyr)
library(tidyr)
library(readxl)
library(lme4)
library(emmeans)
library(ggplot2)
library(flexplot)
library(patchwork)
library(stringr)
library(ggseg)
library(purrr)
library(tibble)
library(readr)
library(gghalves)

bnk_visual <- read_csv("~/Desktop/BNK_BBRnoBBR_plot/bnk_visual_dice_cp.csv")

# if QC_aug26 is NA, then get the value from visual column (previous QC)
bnk_visual <- bnk_visual %>%
  mutate(QC = if_else(is.na(QC_aug26), visual_rating, QC_aug26)) %>%
  select(-visual_rating)

bnk_dice_bbr <- read_csv("~/Desktop/BNK_BBRnoBBR_plot/qc_metrics_rosmap_BBR.csv") %>%
  mutate(sub_ses = str_c(sub_id, ses_id, sep = "_"))
# discard duplicate sub_ses rows
bnk_dice_bbr <- bnk_dice_bbr %>%
  distinct(sub_ses, .keep_all = TRUE) %>%
  filter(protocol == "BNK20090211") %>%
  select(sub_ses, dice) %>%
  rename(dice_bbr = dice)

bnk_dice_nobbr <- read_csv("~/Desktop/BNK_BBRnoBBR_plot/onlybnk_dice_2605_noBBR.csv") %>%
  rename(dice_nobbr = dice)

BNK_dices <- bnk_visual %>%
  right_join(bnk_dice_bbr, by = "sub_ses") %>%
  left_join(bnk_dice_nobbr, by = "sub_ses")

BNK_dice_long <- BNK_dices %>%
  pivot_longer(cols = c(dice_bbr, dice_nobbr), names_to = "version", values_to = "dice") %>%
  mutate(version = recode(version, dice_bbr = "BBR", dice_nobbr = "No BBR")) %>%
  select(-c(contains("v25"))) %>%
  filter(!is.na(dice))

# bc all no-bbr passed, visual qc performed
BNK_dice_long <- BNK_dice_long %>%
  mutate(
    plot_qc = if_else(version == "No BBR", "pass", QC)
  )

ggplot(BNK_dice_long, aes(x = version, y = dice)) +

  geom_half_violin(
    side = c("l", "r"),
    color = "grey60",
    fill = "grey60",
    alpha = 0.1,
    linewidth = 0.3,
    trim = FALSE
  ) +

  geom_line(
  aes(group = sub_ses),
  color = "grey70",
  alpha = 0.4,
  linewidth = 0.4
) +
# NA points: transparent fill, grey outline
geom_point(
  data = subset(BNK_dice_long, is.na(plot_qc)),
  shape = 21,
  fill = NA,
  color = "grey70",
  stroke = 0.7,
  size = 2.5
) +
# pass/fail points
geom_point(
  data = subset(BNK_dice_long, !is.na(plot_qc)),
  aes(color = plot_qc),
  size = 2.5
) +

scale_color_manual(
  values = c(
    "pass" = "#28A64B",
    "fail" = "#A62843"
  ),
  na.value = "grey70",
  name = "Visual QC",
) +
  labs(
    title = "Dice coefficient comparison of anat and func masks in MNI space",
    x = "Version",
    y = "Dice Coefficient"
  ) +
  theme_minimal() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title = element_text(size=18, face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18)
      )



### === now coverage =================================

# count of NA parcels per ID, bbr vs nobbr

bbr_count <- read_tsv("~/Desktop/BNK_BBRnoBBR_plot/bbr_count_non-coverage_456parcels_priority.tsv") %>%
  rename(sub_ses = sub_ses, bbr_count = zero_count)

nobbr_count <- read_tsv("~/Desktop/BNK_BBRnoBBR_plot/nobbr_count_non-coverage_456parcels_priority.tsv") %>%
  rename(sub_ses = sub_ses, nobbr_count = zero_count)

counts <- bbr_count %>%
  left_join(nobbr_count, by = "sub_ses") %>%
  pivot_longer(cols = c(bbr_count, nobbr_count), names_to = "version", values_to = "count") %>%
  mutate(version = recode(version, bbr_count = "BBR", nobbr_count = "No BBR"))

ggplot(counts, aes(x = version, y = count)) +
  geom_half_violin(
    side = c("l", "r"),
    color = "grey60",
    fill = "grey60",
    alpha = 0.1,
    linewidth = 0.3,
    trim = FALSE
  ) +
  geom_line(
    aes(group = sub_ses),
    color = "grey70",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = subset(counts, is.na(count)),
    shape = 21,
    fill = NA,
    color = "grey70",
    stroke = 0.7,
    size = 2.5
  ) +
  geom_point(
    data = subset(counts, !is.na(count)),
    aes(color = version),
    size = 2.5
  ) +
  scale_color_manual(
    values = c(
      "BBR" = "#298c8c",
      "No BBR" = "#f1a226"
    ),
    na.value = "grey70",
    name = "Version"
  ) +
  labs(
    title = "Count of NA parcels per ID",
    x = "Version",
    y = "Count of NA Parcels"
  ) +
  theme_minimal()

# nobbr: sub-06129174_ses-0, 59

# avg coverage per ID, bbr vs nobbr

avg_bbr <- read_tsv("~/Desktop/BNK_BBRnoBBR_plot/bbr_mean_coverage_456parcels_priority.tsv") %>%
  rename(sub_ses = sub_ses, bbr_avg = mean_coverage_percent)

avg_nobbr <- read_tsv("~/Desktop/BNK_BBRnoBBR_plot/nobbr_mean_coverage_456parcels_priority.tsv") %>%
  rename(sub_ses = sub_ses, nobbr_avg = mean_coverage_percent)

average <- avg_bbr %>%
  left_join(avg_nobbr, by = "sub_ses") %>%
  pivot_longer(cols = c(bbr_avg, nobbr_avg), names_to = "version", values_to = "avg_coverage") %>%
  mutate(version = recode(version, bbr_avg = "BBR", nobbr_avg = "No BBR"))

ggplot(average, aes(x = version, y = avg_coverage)) +
  geom_half_violin(
    side = c("l", "r"),
    color = "grey60",
    fill = "grey60",
    alpha = 0.1,
    linewidth = 0.3,
    trim = FALSE
  ) +
  geom_line(
    aes(group = sub_ses),
    color = "grey70",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = subset(average, is.na(avg_coverage)),
    shape = 21,
    fill = NA,
    color = "grey70",
    stroke = 0.7,
    size = 2.5
  ) +
  geom_point(
    data = subset(average, !is.na(avg_coverage)),
    aes(color = version),
    size = 2.5
  ) +
  scale_color_manual(
    values = c(
      "BBR" = "blue",
      "No BBR" = "orange"
    ),
    na.value = "grey70",
    name = "Version"
  ) +
  labs(
    title = "Average coverage per ID",
    x = "Version",
    y = "Average Coverage (%)"
  ) +
  theme_minimal()
