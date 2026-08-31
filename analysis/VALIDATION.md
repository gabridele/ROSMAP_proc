# Validation notes

Validation performed during the publication-readiness cleanup of the supplied
`analysis/` folder.

## Motion-threshold audit

The canonical exclusion is `mean_FD < 0.25`, defined as `FD_THRESHOLD <- 0.25`
in `paths.R`.

A repository-wide search of executable analysis files found no remaining
`mean_FD <= 0.25`, `mean_FD > 0.25`, or `mean_FD >= 0.25` exclusion expressions.
The explicit post-FD filters are in:

- `1_prepare_analysis_table.R` (publication sample summaries)
- `plot_swimmer.R` (publication imaging cohort visualization)
- `qc_plots/counts.R`
- baseline/longitudinal within- and between-network covariate sensitivity scripts
- baseline/longitudinal within- and between-network FD sensitivity scripts
- `fc_related/plot_longitudinal_dmn.R`

Pre-FD sensitivity models intentionally remain unfiltered. Gradient-method values
of `0.25` are unrelated to framewise displacement.

## Python

- `python -m compileall -q analysis` passes.
- Every code cell in `qc_plots/qc_bbr_vs_nobbr_alignment.ipynb` compiles as
  Python and the notebook contains zero stored outputs.
- `0_compute_fc_meanconn.py` passed a synthetic two-scan, 456-parcel smoke test:
  it produced 29 summary columns (ID + 7 within-network + 21 between-network),
  individual 456 x 456 matrices, and a 456 x 456 average matrix with unit diagonal.

The execution environment used for this cleanup does not contain `nibabel` /
`nilearn`, so the gradient-image CLIs could be syntax-checked but not executed.
Install `requirements.txt` in the publication environment before running them.

## R

R is not installed in the cleanup environment, so the R scripts could not be
parsed by `Rscript` or executed against the private study inputs. Static checks
found balanced delimiters in every R file and no remaining personal absolute
paths or stale `analysis/april26` helper references.

Before release, run every publication-relevant R entry point in the actual
project R environment and record exact package versions (ideally with `renv`).

## Data-dependent checks still required

- Confirm final participant/session counts after `mean_FD < 0.25` against the
  manuscript tables when available.
- Confirm the four FC matrices supplied to the gradient manifest were themselves
  generated from the intended post-QC/post-FD scan set.
- Run the BBR/no-BBR QC scripts with the private QC tables/images.
- Verify provenance/licensing/citation requirements for third-party atlas and
  Margulies gradient assets under `gradients/`.
