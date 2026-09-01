# Covariate and motion-effect analyses

This directory contains the statistical analyses used to assess demographic,
acquisition, diagnosis, and motion effects on within- and between-network
functional connectivity.

## Shared configuration and utilities

Entry-point scripts source, in order:

1. `analysis/paths.R` — repository paths, output paths, and the canonical
   framewise-displacement threshold (`FD_THRESHOLD <- 0.25`).
2. `analysis/effect_covs/covs_utils.R` — general packages, network names/colours, network
   pair names, and the list of categorical covariates. And other helpers specific to this folder's scripts.

The local `covs_utils.R` contains repeated programming infrastructure only: baseline
selection, strict FD filtering, long-format conversion, model loops, emmeans
extraction, partial-residual extraction, FD coefficient/FDR tables, and common
plotting/saving helpers. Statistical formulas remain in the individual scripts.

## Input

All model scripts use the prepared connectivity/demographics table returned by:

```r
demos("demos_conn.csv")
```
The table is the result of running `analysis/0_compute_fc_meanconn.py` to get 
the connectivity scores, and `analysis/1_prepare_analysis_table.R` to merge it 
with demographics information.

The scripts expect the seven within-network columns defined by `target_cols`,
the 21 between-network columns defined by `target_combos`, and the covariates
used by the formulas below.

## Motion-exclusion rule

Mean framewise-displacement (FD) threshold for exclusion:

```r
mean_FD < 0.25
```

All post-filtering samples in this directory are created through
`ec_apply_fd_filter()`. Pre-filtering analyses retain
the full available FD range by design.

## Entry-point scripts

| Script | Analysis |
|---|---|
| `effect_covariates_within_network_baseline.R` | Baseline categorical-covariate effects on the 7 within-network measures |
| `effect_covariates_between_network_baseline.R` | Baseline categorical-covariate effects on the 21 between-network measures |
| `effect_covariates_within_network_longitudinal.R` | Longitudinal categorical-covariate effects on within-network connectivity |
| `effect_covariates_between_network_longitudinal.R` | Longitudinal categorical-covariate effects on between-network connectivity |
| `effect_fd_within_network_baseline.R` | Baseline adjusted association between mean FD and within-network connectivity |
| `effect_fd_between_network_baseline.R` | Baseline adjusted association between mean FD and between-network connectivity |
| `effect_fd_within_network_longitudinal.R` | Longitudinal adjusted FD sensitivity analysis for within-network connectivity |
| `effect_fd_between_network_longitudinal.R` | Longitudinal adjusted FD sensitivity analysis for between-network connectivity |
| `plot_covariate_effect_heatmaps.R` | Heatmaps from the longitudinal categorical-covariate result tables |
| `utils.R` | Shared implementation helpers for this directory; not analysis |

Baseline scripts choose the earliest numeric session parsed from `ses_id` for
each participant. Longitudinal scripts retain repeated observations and use a
participant random intercept where shown below. 
Only longitudinal analyses are reported in the manuscript.

## Model specifications

### Baseline

```r
connectivity ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin
```

### Longitudinal

```r
connectivity ~ mean_FD + msex + site + age_scandate +
  eyes + dcfdx + syn_bin + (1 | sub_id)
```

These models are fitted independently to each connectivity
outcome, whether between-network combinations or within-network regions.
Mixed-effects models have participant-specific random intercepts.

Each connectivity-outcome model is fitted once per pre/post-FD dataset and reused for estimated marginal means, population-level fitted-value plots, and partial-residual analyses.

## Multiple-comparison handling

Categorical covariate analyses retain several inferential quantities because they address different multiple-testing questions and should not be
interchanged.

For each categorical covariate, estimated marginal means and pairwise contrasts are obtained using `emmeans`.

The exported result tables contain:

`p_raw`: unadjusted p-value for the pairwise contrast.
`p_tukey`: Tukey-adjusted p-value for pairwise comparisons within a single
fitted connectivity outcome.
`q_across`: Benjamini-Hochberg false-discovery-rate (FDR) adjusted p-value
across all connectivity outcomes and pairwise contrasts for the corresponding covariate, analysis type, and FD-filtering condition.

`sig_tukey`: significance label derived from `p_tukey`.
`sig_across`: significance label derived from `q_across`.

Thus, Tukey adjustment controls multiplicity among factor-level comparisons
within an individual connectivity model, whereas the BH-FDR correction addresses the broader multiplicity introduced by testing the covariate across multiple connectivity outcomes.

Primary significance annotations in the covariate-effect plots and heatmaps use `q_across` / `sig_across`.

Significance labels are:

`*`: q < 0.05
`**`: q < 0.01
`***`: q < 0.001

The Tukey-adjusted values are retained in the output tables so within-outcome pairwise inference remains available separately.

Raw p-values are corrected using the Benjamini-Hochberg procedure across the seven within-network outcomes or the 21 between-network outcomes within the corresponding analysis.


## Predicted values and partial residuals

Longitudinal models include a participant-specific random intercept:

(1 | sub_id)

to account for repeated observations within participants.

For visualization, however, longitudinal fitted-value plots use
population-level fixed-effect predictions. Participant-specific random
intercepts are excluded from the plotted predictions (`re.form` = NA).

Where partial-residual plots are generated, they likewise use fixed-effect
predictions for mixed-effects models. These plots are intended to visualize the adjusted relationship between the covariate of interest and connectivity after accounting for the remaining model covariates.

## Heatmaps

`plot_covariate_effect_heatmaps.R` searches recursively below the configured analysis output directory for:

`all_pairwise_results_withinconn_covs_long.csv`
`all_pairwise_results_betweenconn_covs_long.csv`

Within- and between-network statistics are combined into network-by-network
heatmaps.

Binary covariates (msex, eyes, and syn_bin) are displayed as individual
contrast heatmaps.

Site and diagnosis, which contain more than two levels, are displayed as
lower-triangle contrast grids using the prespecified factor ordering defined in `covs_utils.R`.

Heatmap colour represents the model contrast t-statistic. Pre- and post-FD
versions of each covariate share the same symmetric t-value limits so that
colour magnitude is directly comparable across filtering conditions.

Each populated cell displays the t-statistic. Significance stars appended to the statistic are derived from sig_across, and therefore indicate results surviving the across-connectivity Benjamini-Hochberg FDR correction rather than the within-outcome Tukey correction.

## Suggested run order

The scripts are analysis entry points rather than one monolithic pipeline. For heatmaps, run both longitudinal categorical scripts first so their result CSVs exist, then run `plot_covariate_effect_heatmaps.R`.
The FD scripts are independent of the categorical heatmap step.