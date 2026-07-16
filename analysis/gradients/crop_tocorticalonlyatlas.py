from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

ATLAS_TSV_PATH = Path(
    "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/"
    "atlas_data/atlas-4S456Parcels/atlas-4S456Parcels_dseg.tsv"
)

ATLAS_NIFTI_PATH = Path(
    "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/"
    "atlas-4S456Parcels/atlas-4S456Parcels_space-MNI152NLin6Asym_res-2_dseg.nii.gz"
)

OUTPUT_NIFTI_PATH = Path(
    "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/atlas_data/"
    "atlas-4S456Parcels/atlas-400Parcels_space-MNI152NLin6Asym_res-2_dseg.nii.gz"
)

OUTPUT_TSV_PATH = Path(
    "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/gradients/"
    "atlas_data/atlas-4S456Parcels/atlas-400Parcels_dseg.tsv"
)

EXPECTED_N_PARCELS = 456

# False preserves the original parcel numbers.
# True renumbers retained parcels as 1, 2, 3, ...
RELABEL_CONSECUTIVELY = False


# =============================================================================
# Helpers
# =============================================================================

def find_parcel_id_column(atlas_df: pd.DataFrame) -> str | None:
    """
    Find the TSV column containing integer NIfTI parcel values.
    """

    candidates = [
        "index",
        "parcel_id",
        "label_id",
        "value",
        "id",
    ]

    for candidate in candidates:
        if candidate in atlas_df.columns:
            return candidate

    return None


def load_atlas_table(tsv_path: Path) -> tuple[pd.DataFrame, str]:
    """
    Load the atlas TSV and identify the parcel-ID column.
    """

    atlas_df = pd.read_csv(tsv_path, sep="\t")

    if "network_label" not in atlas_df.columns:
        raise ValueError(
            "Atlas TSV must contain a 'network_label' column."
        )

    parcel_id_column = find_parcel_id_column(atlas_df)

    if parcel_id_column is None:
        print(
            "Warning: no explicit parcel-ID column was found. "
            "Assuming TSV row 1 corresponds to NIfTI label 1, "
            "row 2 to label 2, and so on."
        )

        parcel_id_column = "index"
        atlas_df.insert(
            0,
            parcel_id_column,
            np.arange(1, len(atlas_df) + 1),
        )

    atlas_df[parcel_id_column] = pd.to_numeric(
        atlas_df[parcel_id_column],
        errors="raise",
    ).astype(int)

    if atlas_df[parcel_id_column].duplicated().any():
        duplicates = atlas_df.loc[
            atlas_df[parcel_id_column].duplicated(),
            parcel_id_column,
        ].tolist()

        raise ValueError(
            f"Duplicate parcel IDs in TSV: {duplicates[:20]}"
        )

    return atlas_df, parcel_id_column


def remove_nan_network_parcels(
    atlas_df: pd.DataFrame,
    parcel_id_column: str,
    atlas_nifti_path: Path,
    output_nifti_path: Path,
    output_tsv_path: Path,
    relabel_consecutively: bool = False,
) -> tuple[nib.Nifti1Image, pd.DataFrame]:
    """
    Replace parcels with missing network labels by background value 0.
    """

    valid_mask = atlas_df["network_label"].notna()

    valid_df = atlas_df.loc[valid_mask].copy()
    removed_df = atlas_df.loc[~valid_mask].copy()

    valid_ids = valid_df[parcel_id_column].to_numpy(dtype=int)
    removed_ids = removed_df[parcel_id_column].to_numpy(dtype=int)

    print("\n--- Atlas label summary ---")
    print(f"Total TSV parcels:       {len(atlas_df)}")
    print(f"Retained parcels:        {len(valid_df)}")
    print(f"Removed NaN parcels:     {len(removed_df)}")

    if len(removed_ids) > 0:
        print(
            "Removed NIfTI label values",
            removed_ids,
        )

    # Load the label atlas.
    atlas_img = nib.load(atlas_nifti_path)

    # Preserve the original integer datatype.
    atlas_data = np.asanyarray(atlas_img.dataobj).copy()

    nifti_labels = np.unique(atlas_data)
    nifti_nonzero_labels = nifti_labels[nifti_labels != 0].astype(int)

    print(f"Nonzero labels in NIfTI: {len(nifti_nonzero_labels)}")

    # Check for IDs listed in the TSV but absent from the NIfTI.
    missing_from_nifti = np.setdiff1d(
        atlas_df[parcel_id_column].to_numpy(dtype=int),
        nifti_nonzero_labels,
    )

    if len(missing_from_nifti) > 0:
        print(
            "Warning: TSV parcel IDs absent from the NIfTI",
            missing_from_nifti[:20],
        )

    # Check for labels present in the NIfTI but absent from the TSV.
    missing_from_tsv = np.setdiff1d(
        nifti_nonzero_labels,
        atlas_df[parcel_id_column].to_numpy(dtype=int),
    )

    if len(missing_from_tsv) > 0:
        print(
            "Warning: NIfTI labels absent from the TSV",
            missing_from_tsv[:20],
        )

    # Set NaN-network parcels to background.
    cleaned_data = atlas_data.copy()
    cleaned_data[np.isin(cleaned_data, removed_ids)] = 0

    if relabel_consecutively:
        print("\nRelabeling retained parcels consecutively...")

        # Use a new array so that old and new IDs cannot collide.
        relabeled_data = np.zeros(
            cleaned_data.shape,
            dtype=np.int32,
        )

        old_to_new = {
            int(old_id): int(new_id)
            for new_id, old_id in enumerate(valid_ids, start=1)
        }

        for old_id, new_id in old_to_new.items():
            relabeled_data[cleaned_data == old_id] = new_id

        cleaned_data = relabeled_data

        valid_df["original_index"] = valid_df[parcel_id_column]
        valid_df[parcel_id_column] = np.arange(
            1,
            len(valid_df) + 1,
        )

    # Create output image while retaining spatial metadata.
    output_header = atlas_img.header.copy()
    output_header.set_data_dtype(cleaned_data.dtype)

    cleaned_img = nib.Nifti1Image(
        cleaned_data,
        affine=atlas_img.affine,
        header=output_header,
    )

    # Explicitly preserve qform and sform information.
    cleaned_img.set_qform(
        atlas_img.get_qform(),
        code=int(atlas_img.header["qform_code"]),
    )

    cleaned_img.set_sform(
        atlas_img.get_sform(),
        code=int(atlas_img.header["sform_code"]),
    )

    output_nifti_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_tsv_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    nib.save(cleaned_img, output_nifti_path)

    valid_df.to_csv(
        output_tsv_path,
        sep="\t",
        index=False,
    )

    final_labels = np.unique(cleaned_data)
    final_labels = final_labels[final_labels != 0]

    print("\n--- Output ---")
    print(f"Output parcels: {len(final_labels)}")
    print(f"NIfTI: {output_nifti_path}")
    print(f"TSV:   {output_tsv_path}")

    return cleaned_img, valid_df


# =============================================================================
# Main
# =============================================================================

def main() -> None:

    atlas_df, parcel_id_column = load_atlas_table(
        ATLAS_TSV_PATH
    )

    if len(atlas_df) != EXPECTED_N_PARCELS:
        raise ValueError(
            f"Atlas TSV has {len(atlas_df)} parcels; "
            f"expected {EXPECTED_N_PARCELS}."
        )

    cleaned_img, cleaned_table = remove_nan_network_parcels(
        atlas_df=atlas_df,
        parcel_id_column=parcel_id_column,
        atlas_nifti_path=ATLAS_NIFTI_PATH,
        output_nifti_path=OUTPUT_NIFTI_PATH,
        output_tsv_path=OUTPUT_TSV_PATH,
        relabel_consecutively=RELABEL_CONSECUTIVELY,
    )

    print("\n--- Networks retained ---")

    network_counts = (
        cleaned_table["network_label"]
        .value_counts()
        .sort_index()
    )

    for network, count in network_counts.items():
        print(f"{network}: {count} parcels")


if __name__ == "__main__":
    main()