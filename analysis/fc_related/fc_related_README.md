# Functional-connectivity analyses

This directory contains the following scripts:

- `plot_average_fc.R` — plots site-level average parcelwise functional-connectivity matrices.
- `plot_longitudinal_dmn.R` — models and visualizes longitudinal Default Mode Network (DMN) within-network connectivity.

## Shared configuration

Both scripts source `analysis/paths.R`.

Relevant shared settings include:

- `FD_THRESHOLD <- 0.25`; post-FD analyses retain scans only when `mean_FD < 0.25`.
- `ROSMAP_OUTPUT`; controls where generated tables and figures are written.
- `ROSMAP_DEMOS_CSV`; optionally overrides the prepared `demos_conn.csv` location.
- `ROSMAP_FC_MATRIX_MANIFEST`; optionally overrides the FC-matrix manifest location.
- `SCHAEFER4S456_ATLAS_TSV`; optionally overrides the atlas-label TSV used by `plot_average_fc.R`.

Run scripts from the repository root unless `ROSMAP_ANALYSIS_DIR` is set explicitly.

---

## Average FC matrices

### Script

```text
plot_average_fc.R
```

### Purpose

Plots average parcelwise FC matrices using a common parcel order and a fixed connectivity scale of `[-1, 1]`, allowing direct visual comparison between matrices.

### Matrix manifest

By default, `analysis/paths.R` points to:

```text
analysis/fc_related/fc_matrix_manifest.csv
```

This file is intentionally local/gitignored. It must contain two columns:

```csv
name,path
BNK,/path/to/BNK/avg_fc_matrix.npy
UC,/path/to/UC/avg_fc_matrix.npy
MG,/path/to/MG/avg_fc_matrix.npy
RIRC,/path/to/RIRC/avg_fc_matrix.npy
```

- `name` is the matrix label used in figures and output filenames.
- `path` points to the corresponding NumPy `.npy` matrix.
- Relative paths are resolved relative to the manifest file itself.

The manifest can be located elsewhere by setting:

```bash
export ROSMAP_FC_MATRIX_MANIFEST=/path/to/fc_matrix_manifest.csv
```

### Atlas labels

The script requires a TSV containing one row per parcel and at least:

- `network_label`
- `atlas_name`

The atlas TSV path is configured by `SCHAEFER4S456_ATLAS_TSV` in `analysis/paths.R` and can be overridden with the environment variable of the same name.

The number of TSV rows must exactly match the matrix dimension.

### ROI ordering

By default:

```r
REORDER_BY_NETWORK <- FALSE
```

so the atlas/dseg parcel order is preserved. Set this to `TRUE` only if the supplied TSV is not already grouped by network and grouped network blocks are desired in the figure.

### Outputs

The script writes labelled SVG and clean PNG versions of each matrix under the configured output directory. It also writes:

```text
fc_matrices/average_fc_matrix_summary.csv
```

with the number of finite values, mean, SD, median, minimum and maximum for each matrix. Diagonal/self-connections are excluded from these summaries by default.

### Dependencies

R:

- `dplyr`
- `ggplot2`
- `readr`
- `reticulate`

Python, accessed through `reticulate`:

- `numpy`

---

## Longitudinal DMN trajectories

### Script

```text
plot_longitudinal_dmn.R
```

### Purpose

Models longitudinal within-network connectivity for the Default Mode Network (`Default`) and visualizes participant-specific fitted trajectories before and after FD filtering.

### Input

The script uses the prepared analysis table resolved by:

```r
demos("demos_conn.csv")
```

Required variables are:

- `sub_id`
- `years_from_baseline`
- `mean_FD`
- `msex`
- `site`
- `age_bl`
- `eyes`
- `dcfdx`
- `syn_bin`
- `Default`

Rows with missing model variables are excluded. Participants must have at least two eligible observations in the corresponding pre- or post-FD dataset.

The post-FD sample additionally requires:

```r
mean_FD < FD_THRESHOLD
```

with `FD_THRESHOLD = 0.25` in `analysis/paths.R`.

### Factor coding

Explicit factor ordering is used:

- `msex`: `female`, `male`
- `site`: `BNK`, `UC`, `MG`, `RIRC`
- `eyes`: `closed`, `open`
- `dcfdx`: `NCI`, `MCI`, `AD`, `other`
- `syn_bin`: `not SyN`, `SyN`

Unexpected levels cause an error rather than being silently accepted.

### Model

The same specification is fitted separately to the pre- and post-FD datasets:

```r
within_conn ~
  years_from_baseline +
  mean_FD +
  msex +
  site +
  age_bl +
  eyes +
  dcfdx +
  syn_bin +
  (1 + years_from_baseline | sub_id)
```

The model includes participant-specific random intercepts and random slopes for follow-up time. `lme4::isSingular()` is checked after fitting and singular fits trigger a warning.

### Subject-specific fitted trajectories

The plotted curves are **conditional participant-specific predictions**. Predictions use:

```r
re.form = NULL
```

so fitted participant random intercepts and slopes are included.

For each participant, predictions span that participant's observed follow-up interval. While time varies:

- `mean_FD` is held at the participant's mean across modelled visits;
- `age_bl` is held at the participant's baseline-age value;
- `msex`, `site`, `eyes`, `dcfdx`, and `syn_bin` are held at their first observed values in the model data.

The resulting curves are therefore standardized subject-specific fitted trajectories.

Plots grouped by last-visit diagnosis or site-change status use those variables descriptively after prediction; the prediction grid itself does not change diagnosis or site over follow-up.

### Trajectory direction

Participants are descriptively classified as `up`, `flat`, or `down` according to fitted change between the beginning and end of their predicted trajectory.

The `flat` interval is defined as plus/minus one standard deviation of participant fitted changes within the corresponding pre- or post-FD dataset.

### Outputs

Tables:

```text
fc_related/dmn_longitudinal_sample_summary.csv
fc_related/dmn_model_coefficients_pre_fd.csv
fc_related/dmn_model_coefficients_post_fd.csv
fc_related/dmn_subject_trajectory_summary.csv
```

Main figures:

```text
fc_related/predicted_dmn_pre_fd_direction.pdf
fc_related/predicted_dmn_post_fd_direction.pdf
```

Additional pre-FD descriptive figures group the fitted trajectories by first-visit site, site-change status, and last-visit diagnosis.

### Dependencies

- `dplyr`
- `ggplot2`
- `lme4`
- `lmerTest`
- `readr`
- `tidyr`

---

## Suggested run order

The scripts are independent.

Before running `plot_average_fc.R`, ensure the FC manifest and atlas TSV resolve to valid local files. Before running `plot_longitudinal_dmn.R`, ensure the prepared `demos_conn.csv` table has been generated and is resolvable through `demos()`.
