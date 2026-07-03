from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.maskers import NiftiLabelsMasker


def load_and_validate_labels(labels_path: Path) -> pd.DataFrame:
    """Load parcel labels and validate the expected columns."""
    df = pd.read_table(labels_path)

    required_columns = {"index", "label"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(
            f"Label file is missing required columns: {sorted(missing_columns)}"
        )

    df = (
        df.sort_values("index")
        .reset_index(drop=True)
        .rename(columns={"label": "name"})
    )

    if df["index"].duplicated().any():
        duplicated = df.loc[df["index"].duplicated(), "index"].tolist()
        raise ValueError(f"Duplicate parcel indices found: {duplicated}")

    return df


def extract_timeseries(
    atlas_path: Path,
    labels_path: Path,
    bold_path: Path,
) -> pd.DataFrame:
    """Extract parcel time series and align columns to the full label list."""
    atlas_img = nib.load(str(atlas_path))
    df_parcel_labels = load_and_validate_labels(labels_path)

    node_labels = df_parcel_labels["name"].tolist()
    n_nodes = len(node_labels)

    # Map atlas parcel value -> final output column index
    full_parcel_mapper = {
        int(parcel_value): col_idx
        for col_idx, parcel_value in enumerate(df_parcel_labels["index"].tolist())
    }

    atlas_values = np.unique(atlas_img.get_fdata())
    atlas_values = atlas_values[atlas_values != 0].astype(int)

    masker_lut = (
        df_parcel_labels.loc[
            df_parcel_labels["index"].isin(atlas_values),
            ["index", "name"],
        ]
        .reset_index(drop=True)
    )

    if masker_lut.empty:
        raise ValueError("No atlas labels in the image matched the label table.")

    masker = NiftiLabelsMasker(
        labels_img=atlas_img,
        lut=masker_lut,
        background_label=0,
        mask_img=None,
        smoothing_fwhm=None,
        standardize=False,
        resampling_target=None,
        keep_masked_labels=True,
        strategy="mean",
    )

    timeseries_arr = masker.fit_transform(str(bold_path))

    if timeseries_arr.ndim == 1:
        timeseries_arr = timeseries_arr[None, :]

    masker_parcel_mapper = {
        col_idx: atlas_value
        for col_idx, atlas_value in masker.region_ids_.items()
        if col_idx != "background"
    }

    final_arr = np.full(
        (timeseries_arr.shape[0], n_nodes),
        np.nan,
        dtype=timeseries_arr.dtype,
    )

    for col_idx, atlas_value in masker_parcel_mapper.items():
        atlas_value = int(atlas_value)
        full_col_idx = full_parcel_mapper.get(atlas_value)
        if full_col_idx is not None:
            final_arr[:, full_col_idx] = timeseries_arr[:, col_idx]

    return pd.DataFrame(final_arr, columns=node_labels)


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "Usage: python script.py <atlas> <labels> <denoised_bold> <output>",
            file=sys.stderr,
        )
        return 1

    atlas_path = Path(sys.argv[1])
    labels_path = Path(sys.argv[2])
    bold_path = Path(sys.argv[3])
    output_path = Path(sys.argv[4])

    df = extract_timeseries(
        atlas_path=atlas_path,
        labels_path=labels_path,
        bold_path=bold_path,
    )
    df.to_csv(output_path, sep="\t", na_rep="n/a", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())