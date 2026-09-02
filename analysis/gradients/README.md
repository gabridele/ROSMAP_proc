# Functional-connectivity gradients

This directory contains the scripts used to derive functional-connectivity gradients from site-level average ROSMAP connectivity matrices.

The gradient implementation in `compute_gradients.R` and `util_gradients.R` was adapted from the BioFINDER gradient-analysis workflow developed by Jorrit Rittmo:

https://github.com/DeMONLab-BioFINDER/fc_changes_follow_gradients

The numerical gradient implementation has intentionally been kept close to the original fork. Repository-specific modifications primarily concern input paths, local configuration, and ROSMAP atlas/reference preparation rather than the underlying gradient methodology.

## Directory structure

```text
gradients/
├── compute_gradients.R
├── util_gradients.R
├── prepa_grad_atlas.R
├── crop_cortical_atlas.py
│
├── atlas_data/
│   ├── Schaefer2018_400Parcels_7Networks_order.txt
│   ├── Schaefer2018_400Parcels_order.txt
│   ├── org_mask_yeo_400.txt
│   ├── schaef400_ggseg2.rds
│   └── atlas-4S456Parcels/
│
└── marg/
    ├── gradient1.nii.gz
    ├── gradient2.nii.gz
    ├── gradient3.nii.gz
    ├── gradient4.nii.gz
    ├── gradient5.nii.gz
    └── volumetric_margcsv.py
```

---

## Main analysis

The main script is:

```text
compute_gradients.R
```

It performs the following steps:

1. loads Schaefer-400 parcel and Yeo-7-network metadata;
2. loads the Margulies reference gradients;
3. loads site-level average FC matrices for BNK, UC, MG, and RIRC;
4. restricts the FC matrices to the first 400 cortical parcels;
5. computes gradients independently for each site;
6. aligns gradient directions to the Margulies reference;
7. generates gradient visualizations;
8. exports gradient values and variance explained.

The numerical functions used for gradient estimation are defined in:

```text
util_gradients.R
```

---

## Site-level FC matrices

The site-level average FC matrices are configured using the local file:

```text
analysis/gradients/gradient_fc_sheet.csv
```

It must contain:

```csv
name,path
BNK,/path/to/BNK/avg_fc_matrix.npy
UC,/path/to/UC/avg_fc_matrix.npy
MG,/path/to/MG/avg_fc_matrix.npy
RIRC,/path/to/RIRC/avg_fc_matrix.npy
```

The default location is defined in `analysis/paths.R`:

```r
GRADIENT_FC_SHEET <- Sys.getenv(
  "ROSMAP_GRADIENT_FC_SHEET",
  unset = file.path(
    ANALYSIS_DIR,
    "gradients",
    "gradient_fc_sheet.csv"
  )
)
```

A different sheet can therefore be supplied with:

```bash
export ROSMAP_GRADIENT_FC_SHEET=/path/to/gradient_fc_sheet.csv
```

Absolute paths are recommended in the local sheet.

The matrices are loaded with NumPy through `reticulate`.

Matrices are restricted as follows:

```r
matrix[1:400, 1:400]
```

corresponding to the 400 cortical-only Schaefer parcels.

---

## Margulies reference gradients

Gradient orientation is referenced to the canonical Margulies connectivity gradients.

The volumetric Margulies gradient images are stored under:

```text
analysis/gradients/marg/
```

as:

```text
gradient1.nii.gz
gradient2.nii.gz
gradient3.nii.gz
gradient4.nii.gz
gradient5.nii.gz
```

`volumetric_margcsv.py` parcellates these images using the same 400-parcel cortical atlas used by the ROSMAP gradient workflow.

Running:

```bash
python analysis/gradients/marg/volumetric_margcsv.py
```

generates:

```text
analysis/gradients/marg/volumetric_marg_xcpd_400.csv
```

The resulting table contains one row per cortical parcel and columns including:

```text
gradient1
gradient2
gradient3
gradient4
gradient5
```

The primary gradient analysis uses the first three reference gradients:

```r
marg_gradients <- read_csv(
  MARGULIES_GRADIENTS_400
) %>%
  select(
    gradient1,
    gradient2,
    gradient3
  )
```

Its default location can be overridden with:

```bash
export ROSMAP_MARGULIES_GRADIENTS_400=/path/to/reference.csv
```

---

## Atlas preparation

### Cortical 400-parcel atlas

The FC matrices originate from the 4S456 atlas representation. The gradient analysis uses the cortical 400-parcel component.

`crop_cortical_atlas.py` can be used to create a cortical-only atlas by removing parcels without network labels.

The script accepts the atlas TSV and NIfTI as explicit command-line inputs and can optionally verify the expected input and output parcel counts.

For example:

```bash
python analysis/gradients/crop_cortical_atlas.py \
  --atlas-tsv /path/to/atlas-4S456Parcels_dseg.tsv \
  --atlas-nifti /path/to/atlas-4S456Parcels_dseg.nii.gz \
  --output-nifti /path/to/atlas-400Parcels_dseg.nii.gz \
  --output-tsv /path/to/atlas-400Parcels_dseg.tsv \
  --expected-input-parcels 456 \
  --expected-output-parcels 400
```

Unless explicitly requested, retained parcel IDs are not renumbered.

### Schaefer/Yeo plotting metadata

`prep_grad_atlas.R` prepares metadata used by the inherited gradient visualization functions.

It generates:

```text
Schaefer2018_400Parcels_order.txt
org_mask_yeo_400.txt
schaef400_ggseg2.rds
```

The parcel-order file provides the region names used by `util_gradients.R`.

`org_mask_yeo_400.txt` maps each parcel to one of the seven Yeo networks:

1. Visual
2. Somatomotor
3. Dorsal Attention
4. Salience/Ventral Attention
5. Limbic
6. Control
7. Default

`prep_grad_atlas.R` uses:

```text
Schaefer2018_400Parcels_7Networks_order.txt
```

as the source parcel-order file.
---

## Primary gradient specification

The primary gradient analysis is:

```r
grad_list_individ <- get_gradients(
  connectome_ests = list(
    bnk400xcpd = bnk400xcpd,
    uc400xcpd = uc400xcpd,
    mg400xcpd = mg400xcpd,
    rirc400xcpd = rirc400xcpd
  ),
  reference_gradients = marg_gradients,
  n_gradients = c(1, 2, 3),
  threshold = 0.0,
  similarity_method = "cosine",
  on_affinity = FALSE,
  method = "pca",
  visualize = TRUE
)
```

Thus, the primary analysis uses:

- PCA decomposition;
- gradients 1–3;
- threshold parameter `0.0`;
- no affinity-matrix transformation;
- sign alignment to the Margulies reference gradients;
- site-level average FC matrices;
- 400 cortical parcels.

### Treatment of negative connectivity

The inherited implementation applies:

```r
L[L < 0] <- 0
```

after the thresholding step.

Negative functional-connectivity values are therefore set to zero before gradient decomposition.

This behavior is retained from the forked implementation rather than being redefined during repository cleanup.

### Thresholding

Sparsification is implemented in the inherited `zero_out_mat()` function:

```r
thresholds <- colQuantiles(
  abs(mat),
  probs = thresh
)

mat[
  abs(mat) <= thresholds[col(mat)]
] <- 0
```

The primary analysis passes:

```r
threshold = 0.0
```

This README describes the implementation as executed rather than redefining the historical thresholding behavior.

### Gradient alignment

Gradient components can have arbitrary sign.

After decomposition, `align_gradients()` compares each derived component with the corresponding Margulies reference gradient:

```r
if (
  cor(
    original[, i],
    derived[, i]
  ) < 0
) {
  derived[, i] <- -derived[, i]
}
```

Components with negative correlation to the corresponding reference gradient are sign-flipped.

The primary analysis does not reorder components; it aligns their sign to the corresponding Margulies gradients.

---

## Outputs

Primary numerical outputs are written under:

```text
$ROSMAP_OUTPUT/gradients/
```

including:

```text
primary_gradients.csv
primary_variance_explained.csv
```

`primary_gradients.csv` contains the reference and site-derived gradient values produced by `get_gradients()`.

`primary_variance_explained.csv` contains the variance explained by the selected gradient components for each site.

When:

```r
visualize = TRUE
```

the inherited visualization code also generates gradient brain plots under:

```text
analysis/gradients/plots/
```

for the Margulies reference and each site-derived gradient.

---

## Dependencies

### R

The primary R workflow uses packages including:

- `tidyverse`
- `SCORPIUS`
- `conflicted`
- `pbapply`
- `reticulate`
- `matrixStats`
- `scales`
- `patchwork`
- `ggside`
- `ggpmisc`
- `sf`

The atlas-preparation helper additionally uses:

- `ggsegSchaefer`
- `ggseg.formats`

Some packages are loaded inside the inherited `util_gradients.R` functions rather than at the beginning of `compute_gradients.R`.

### Python

`compute_gradients.R` accesses NumPy through `reticulate`.

Required Python packages for the active workflow include:

- `numpy`
- `pandas`
- `nibabel`
- `nilearn`
- `tqdm`

`crop_cortical_atlas.py` additionally uses `nibabel`.

---

## Recommended run order

From the repository root:

### 1. Prepare atlas metadata if required

```bash
Rscript analysis/gradients/prep_grad_atlas.R
```

This step is unnecessary if the required parcel-order, Yeo-mask and ggseg geometry files already exist.

### 2. Create the volumetric Margulies reference

```bash
python analysis/gradients/marg/volumetric_margcsv.py
```

Verify that this produces:

```text
analysis/gradients/marg/volumetric_marg_xcpd_400.csv
```

with 400 rows and the expected gradient columns.

### 3. Create the local FC sheet

Create:

```text
analysis/gradients/gradient_fc_sheet.csv
```

with paths to the BNK, UC, MG and RIRC site-average FC matrices.

### 4. Compute gradients

```bash
Rscript analysis/gradients/compute_gradients.R
```

---