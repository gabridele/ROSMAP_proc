library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggeffects)
library(lme4)
library(lmerTest)
library(lubridate)
library(emmeans)
library(ggpubr)
library(readxl)
library(grDevices)
library(visreg)

emm_options(pbkrtest.limit = 10000)
emm_options(lmerTest.limit = 50000)

fd_threshold <- 0.25

target_cols <- c(
  "Vis", "SomMot", "DorsAttn",
  "SalVentAttn", "Limbic", "Cont", "Default"
)

network_colors <- c(
  "Vis" = "#9B59B6",
  "SomMot" = "#6C8EBF",
  "Default" = "#D36B78",
  "Limbic" = "#C9D39A",
  "DorsAttn" = "#3C8D2F",
  "SalVentAttn" = "#C84CCF",
  "Cont" = "#E5B53A"
)

# Between-network columns
target_combos <- c(
  "Cont_to_Default", "Cont_to_DorsAttn",
  "Cont_to_Limbic", "Cont_to_SalVentAttn", "Cont_to_SomMot",
  "Cont_to_Vis", "Default_to_DorsAttn", "Default_to_Limbic",
  "Default_to_SalVentAttn", "Default_to_SomMot", "Default_to_Vis",
  "DorsAttn_to_Limbic", "DorsAttn_to_SalVentAttn", "DorsAttn_to_SomMot",
  "DorsAttn_to_Vis", "Limbic_to_SalVentAttn", "Limbic_to_SomMot",
  "Limbic_to_Vis", "SalVentAttn_to_SomMot", "SalVentAttn_to_Vis",
  "SomMot_to_Vis"
)

between_network_colors <- c(
  "Cont_to_Default" = "#DC9059",
  "Cont_to_DorsAttn" = "#91A135",
  "Cont_to_Limbic" = "#D7C46A",
  "Cont_to_SalVentAttn" = "#D78185",
  "Cont_to_SomMot" = "#A9A27D",
  "Cont_to_Vis" = "#C08778",
  "Default_to_DorsAttn" = "#887C54",
  "Default_to_Limbic" = "#CE9F89",
  "Default_to_SalVentAttn" = "#CE5CA3",
  "Default_to_SomMot" = "#A07D9C",
  "Default_to_Vis" = "#B76297",
  "DorsAttn_to_Limbic" = "#83B065",
  "DorsAttn_to_SalVentAttn" = "#826D7F",
  "DorsAttn_to_SomMot" = "#548E77",
  "DorsAttn_to_Vis" = "#6B7373",
  "Limbic_to_SalVentAttn" = "#C990B4",
  "Limbic_to_SomMot" = "#9BB1AD",
  "Limbic_to_Vis" = "#B296A8",
  "SalVentAttn_to_SomMot" = "#9A6DC7",
  "SalVentAttn_to_Vis" = "#B252C3",
  "SomMot_to_Vis" = "#8474BB"
)

covariates_to_run <- c(
  "msex",
  "site",
  "eyes",
  "syn_bin",
  "dcfdx"
)

make_long <- function(data) {
  data %>%
    pivot_longer(
      cols = all_of(target_cols),
      names_to = "network",
      values_to = "within_conn"
    ) %>%
    mutate(
      network = factor(network, levels = target_cols)
    ) %>%
    filter(
      !is.na(ses_num),
      !is.na(mean_FD),
      !is.na(within_conn),
      !is.na(msex),
      !is.na(site),
      !is.na(age_scandate),
      !is.na(distortion_correction),
      !is.na(eyes)
    )
}

site_cols <- c(
  "BNK" = "#0072B2",
  "MG" = "#E69F00",
  "RIRC" = "#009E73",
  "UC" = "#d50700"
)