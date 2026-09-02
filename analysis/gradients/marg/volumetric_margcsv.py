from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm
from nilearn.maskers import NiftiLabelsMasker


# ------------------------------------------------------------
# Atlas files
# ------------------------------------------------------------

xcpd_400_atlas_file = (
    "analysis/gradients/"
    "atlas_data/atlas-4S456Parcels/atlas-400Parcels_space-MNI152NLin6Asym_res-2_dseg.nii.gz"
)

# ------------------------------------------------------------
# Create maskers
# ------------------------------------------------------------

maskers = {
    "xcpd_400": NiftiLabelsMasker(
        labels_img=xcpd_400_atlas_file,
        smoothing_fwhm=None,
        standardize=False,
        strategy="mean",
        verbose=0,
    )
}


# ------------------------------------------------------------
# Gradient images
# ------------------------------------------------------------

gradient_dir = Path(
    "/Users/ga0034de/github_dir/ROSMAP_proc/"
    "analysis/gradients/marg"
)

nifti_images = sorted(gradient_dir.glob("*.nii.gz"))

if not nifti_images:
    raise FileNotFoundError(
        f"No .nii.gz files were found in {gradient_dir}"
    )

print("Gradient images:")

for image in nifti_images:
    print(" ", image.name)


# ------------------------------------------------------------
# Extract each gradient with each atlas
# ------------------------------------------------------------

atlas_results = {}

for atlas_name, masker in maskers.items():

    print(f"\nProcessing atlas: {atlas_name}")

    gradient_values = {}
    expected_n_rois = None

    for img_path in tqdm(
        nifti_images,
        desc=f"Extracting {atlas_name}",
    ):

        # Correct call: use the masker object directly
        values = masker.fit_transform(str(img_path))
        values = np.asarray(values).ravel()

        # gradient1.nii.gz -> gradient1
        gradient_name = img_path.name.removesuffix(".nii.gz")

        if expected_n_rois is None:
            expected_n_rois = len(values)

        elif len(values) != expected_n_rois:
            raise ValueError(
                f"Inconsistent number of parcels for {atlas_name}: "
                f"{img_path.name} returned {len(values)}, but the "
                f"first gradient returned {expected_n_rois}."
            )

        gradient_values[gradient_name] = values

    # One row per atlas parcel
    atlas_df = pd.DataFrame(gradient_values)

    atlas_df.insert(
        0,
        "parcel_id",
        np.arange(1, len(atlas_df) + 1),
    )

    atlas_df.insert(
        0,
        "atlas",
        atlas_name,
    )

    atlas_results[atlas_name] = atlas_df

    print(
        f"{atlas_name}: extracted "
        f"{len(atlas_df)} parcels and "
        f"{len(gradient_values)} gradients"
    )


# ------------------------------------------------------------
# Inspect results
# ------------------------------------------------------------

for atlas_name, dataframe in atlas_results.items():
    print(f"\n{atlas_name}")
    print(dataframe.shape)
    print(dataframe.head())
    dataframe.to_csv(
        f"/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/marg/volumetric_marg_{atlas_name}.csv",
        index=False,
    )