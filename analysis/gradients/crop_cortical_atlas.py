"""Create a cortical-only label atlas by removing parcels without network labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


def find_parcel_id_column(atlas_df: pd.DataFrame) -> str | None:
    """Return the first recognized integer parcel-ID column, if present."""
    for candidate in ("index", "parcel_id", "label_id", "value", "id"):
        if candidate in atlas_df.columns:
            return candidate
    return None


def load_atlas_table(tsv_path: Path) -> tuple[pd.DataFrame, str]:
    """Load the atlas TSV and identify/create the NIfTI parcel-ID column."""
    atlas_df = pd.read_csv(tsv_path, sep="\t")

    if "network_label" not in atlas_df.columns:
        raise ValueError("Atlas TSV must contain a 'network_label' column.")

    parcel_id_column = find_parcel_id_column(atlas_df)
    if parcel_id_column is None:
        parcel_id_column = "index"
        atlas_df.insert(0, parcel_id_column, np.arange(1, len(atlas_df) + 1))

    atlas_df[parcel_id_column] = pd.to_numeric(
        atlas_df[parcel_id_column], errors="raise"
    ).astype(int)

    if atlas_df[parcel_id_column].duplicated().any():
        duplicates = atlas_df.loc[
            atlas_df[parcel_id_column].duplicated(), parcel_id_column
        ].tolist()
        raise ValueError(f"Duplicate parcel IDs in TSV: {duplicates[:20]}")

    return atlas_df, parcel_id_column


def create_cortical_atlas(
    atlas_df: pd.DataFrame,
    parcel_id_column: str,
    atlas_nifti_path: Path,
    output_nifti_path: Path,
    output_tsv_path: Path,
    relabel_consecutively: bool = False,
) -> tuple[nib.Nifti1Image, pd.DataFrame]:
    """Set parcels with missing network labels to background and export TSV/NIfTI."""
    network = atlas_df["network_label"]
    valid_mask = network.notna() & network.astype(str).str.strip().ne("")

    valid_df = atlas_df.loc[valid_mask].copy()
    removed_df = atlas_df.loc[~valid_mask].copy()

    valid_ids = valid_df[parcel_id_column].to_numpy(dtype=int)
    removed_ids = removed_df[parcel_id_column].to_numpy(dtype=int)

    atlas_img = nib.load(atlas_nifti_path)
    atlas_data = np.asanyarray(atlas_img.dataobj).copy()

    nifti_nonzero = np.unique(atlas_data)
    nifti_nonzero = nifti_nonzero[nifti_nonzero != 0].astype(int)

    missing_from_nifti = np.setdiff1d(
        atlas_df[parcel_id_column].to_numpy(dtype=int), nifti_nonzero
    )
    missing_from_tsv = np.setdiff1d(
        nifti_nonzero, atlas_df[parcel_id_column].to_numpy(dtype=int)
    )

    if len(missing_from_nifti):
        print("Warning: TSV IDs absent from NIfTI:", missing_from_nifti[:20])
    if len(missing_from_tsv):
        print("Warning: NIfTI IDs absent from TSV:", missing_from_tsv[:20])

    cleaned_data = atlas_data.copy()
    cleaned_data[np.isin(cleaned_data, removed_ids)] = 0

    if relabel_consecutively:
        relabeled = np.zeros(cleaned_data.shape, dtype=np.int32)
        old_to_new = {
            int(old_id): int(new_id)
            for new_id, old_id in enumerate(valid_ids, start=1)
        }
        for old_id, new_id in old_to_new.items():
            relabeled[cleaned_data == old_id] = new_id
        cleaned_data = relabeled
        valid_df["original_index"] = valid_df[parcel_id_column]
        valid_df[parcel_id_column] = np.arange(1, len(valid_df) + 1)

    header = atlas_img.header.copy()
    header.set_data_dtype(cleaned_data.dtype)
    cleaned_img = nib.Nifti1Image(cleaned_data, atlas_img.affine, header)
    cleaned_img.set_qform(
        atlas_img.get_qform(), code=int(atlas_img.header["qform_code"])
    )
    cleaned_img.set_sform(
        atlas_img.get_sform(), code=int(atlas_img.header["sform_code"])
    )

    output_nifti_path.parent.mkdir(parents=True, exist_ok=True)
    output_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(cleaned_img, output_nifti_path)
    valid_df.to_csv(output_tsv_path, sep="\t", index=False)

    final_labels = np.unique(cleaned_data)
    final_labels = final_labels[final_labels != 0]

    print(f"Input TSV parcels: {len(atlas_df)}")
    print(f"Retained cortical parcels: {len(valid_df)}")
    print(f"Removed parcels: {len(removed_df)}")
    print(f"Nonzero labels in output NIfTI: {len(final_labels)}")
    print(f"NIfTI: {output_nifti_path}")
    print(f"TSV:   {output_tsv_path}")

    return cleaned_img, valid_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas-tsv", required=True, type=Path)
    parser.add_argument("--atlas-nifti", required=True, type=Path)
    parser.add_argument("--output-nifti", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument(
        "--expected-input-parcels",
        type=int,
        default=None,
        help="Optional validation of the input TSV row count.",
    )
    parser.add_argument(
        "--expected-output-parcels",
        type=int,
        default=None,
        help="Optional validation of the retained cortical parcel count.",
    )
    parser.add_argument(
        "--relabel-consecutively",
        action="store_true",
        help="Renumber retained labels from 1..N instead of preserving IDs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    atlas_df, parcel_id_column = load_atlas_table(args.atlas_tsv)

    if (
        args.expected_input_parcels is not None
        and len(atlas_df) != args.expected_input_parcels
    ):
        raise ValueError(
            f"Atlas TSV has {len(atlas_df)} rows; expected "
            f"{args.expected_input_parcels}."
        )

    _, cleaned_table = create_cortical_atlas(
        atlas_df=atlas_df,
        parcel_id_column=parcel_id_column,
        atlas_nifti_path=args.atlas_nifti,
        output_nifti_path=args.output_nifti,
        output_tsv_path=args.output_tsv,
        relabel_consecutively=args.relabel_consecutively,
    )

    if (
        args.expected_output_parcels is not None
        and len(cleaned_table) != args.expected_output_parcels
    ):
        raise ValueError(
            f"Cortical atlas contains {len(cleaned_table)} rows; expected "
            f"{args.expected_output_parcels}."
        )

    print("\nNetworks retained:")
    print(cleaned_table["network_label"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
