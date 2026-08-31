# ROSMAP analysis

This directory contains the study-level analyses used to summarize functional
connectivity, assemble the private analysis table, assess nuisance covariates,
produce QC figures, and explore cortical connectivity gradients.

The code is designed for a mixed environment: preprocessing outputs may live on
a SLURM/HPC filesystem, while figures and exploratory work may also be run on a
local workstation. **No ROSMAP participant-level data are distributed here.**

## Analysis flow

1. `0_compute_fc_meanconn.py` computes parcelwise correlation matrices and
   seven-network within/between-network summaries from XCP-D parcel time series.
2. `1_prepare_analysis_table.R` merges those summaries with restricted metadata,
   scan information, diagnosis, and QC variables and writes
   `outputs/prepared/demos_conn.csv` by default.
3. `effect_covs/` contains baseline and longitudinal nuisance/covariate analyses.
4. `fc_related/` contains FC-matrix and longitudinal DMN visualizations.
5. `qc_plots/` contains sample-count and BBR/no-BBR QC analyses.
6. `gradients/` contains an exploratory gradient workflow and its supporting
   utilities/reference assets.
7. `plot_swimmer.R` visualizes diagnosis trajectories for scans retained by the publication motion rule.

These scripts are intentionally separate rather than presented as a single
push-button workflow because they answer distinct analysis/QC questions and
require different private inputs.

## Motion exclusion rule

The publication rule is **strictly**:

```text
mean_FD < 0.25
```

`FD_THRESHOLD <- 0.25` is defined once in `paths.R`. Every scan-exclusion filter
uses `< FD_THRESHOLD`, not `<=`. Some motion-sensitivity scripts deliberately
fit both a pre-filter model (all usable scans) and a post-filter model
(`mean_FD < 0.25`); the pre-filter model is not an accidental omission.

Values of `0.25` inside gradient parameter grids, plot alpha values, axis limits,
or layout code are unrelated to framewise displacement and should not be changed.

## Configuration

Start from `analysis.env.example`. The main variables are:

- `ROSMAP_DATA`: directory containing private input tables.
- `ROSMAP_OUTPUT`: directory for generated tables/figures.
- `ROSMAP_DEMOS_CSV`: optional override for the downstream prepared `demos_conn.csv`.
- `ROSMAP_PREPARED_CSV`: optional output path for `1_prepare_analysis_table.R`.
- `ROSMAP_ANALYSIS_DIR`: optional explicit path to this directory.

If `ROSMAP_DEMOS_CSV` is not set, analyses first look for
`$ROSMAP_OUTPUT/prepared/demos_conn.csv`, making the output from
`1_prepare_analysis_table.R` the default input to downstream analyses.

## Functional-connectivity calculation

Example:

```bash
python analysis/0_compute_fc_meanconn.py \
  --atlas-tsv /path/to/atlas-4S456Parcels_dseg.tsv \
  --timeseries-dir /path/to/xcpd/timeseries \
  --output-csv /path/to/analysis_outputs/atlas_mean_connectivity456.csv
```

The script saves one 456 x 456 correlation matrix per scan and an average matrix,
and summarizes the seven cortical network labels. By default, correlations are
averaged in Fisher-z space before transformation back to correlation space.

## Private inputs

`1_prepare_analysis_table.R` currently expects these files under `ROSMAP_DATA`:

- `derivatives_list_with_age.csv`
- `OLD_mean_within_conn_demos.csv`
- `ROSMAP_demos2026.xlsx`
- `atlas_mean_connectivity456.csv`
- `variables_ses_specific_may26.xlsx`
- `age_atscan.csv`

The input tables are not included because data access is handled separately by
the study/publication. The preparation script writes both unfiltered descriptive
summaries and explicit `*_postFD.csv` summaries for the final `mean_FD < 0.25`
imaging sample.

## Baseline versus longitudinal models

Baseline scripts explicitly select the earliest numeric session per participant.
Longitudinal scripts retain repeated observations and, where indicated in the
script, use participant-level random intercepts. Model formulas and output names
remain in each analysis file so the statistical specification is auditable.

## Subdirectories

- [`effect_covs/README.md`](effect_covs/README.md): covariate/motion analyses.
- [`fc_related/README.md`](fc_related/README.md): FC visualizations and manifests.
- [`qc_plots/README.md`](qc_plots/README.md): QC input requirements.
- [`gradients/README.md`](gradients/README.md): gradient workflow, provenance,
  external dependencies, and generated atlas resources.

## Dependencies and validation

Python dependencies are listed in `requirements.txt`; R packages are listed in
`R_PACKAGES.md`. The scripts have been made path-portable and Python files can be
syntax-checked with:

```bash
python -m compileall -q analysis
```

R should additionally be parsed/run in the project R environment before release,
because static inspection cannot validate installed-package versions or the
private study inputs.
