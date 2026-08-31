#!/usr/bin/env python3
"""Compute parcelwise and network-level functional connectivity for 4S456.

The script consumes XCP-D-derived parcel timeseries and the Schaefer-Tien 4S456
label table. For each scan it:

1. computes a 456 x 456 Pearson-correlation matrix;
2. optionally stores that subject/session matrix;
3. summarizes mean within-network and between-network connectivity for parcels
   assigned to one of the seven cortical networks; and
4. writes one tabular row per subject/session.

Correlations are averaged in Fisher-z space by default and transformed back to
correlation space. Parcels without a ``network_label`` (for example non-cortical
regions) remain in the saved full FC matrices but are excluded from the network
summaries.

The atlas TSV is not distributed with ROSMAP_proc. Obtain it and pass it explicitly
with ``--atlas-tsv``.
"""

from __future__ import annotations

import argparse
import re
from glob import glob
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

EXPECTED_N_PARCELS = 456
DEFAULT_TIMESERIES_PATTERN = "*.tsv"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--atlas-tsv",
        required=True,
        type=Path,
        help="AtlasPack atlas-4S456Parcels_dseg.tsv label table.",
    )
    parser.add_argument(
        "--timeseries-dir",
        required=True,
        type=Path,
        help="Directory containing parcel timeseries TSV files.",
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_TIMESERIES_PATTERN,
        help=f"Glob pattern inside --timeseries-dir (default: {DEFAULT_TIMESERIES_PATTERN!r}).",
    )
    parser.add_argument(
        "--output-csv",
        required=True,
        type=Path,
        help="Destination CSV for within- and between-network summaries.",
    )
    parser.add_argument(
        "--matrix-dir",
        type=Path,
        default=None,
        help=(
            "Directory for subject/session FC matrices and avg_fc_matrix.npy. "
            "Defaults to <output-csv parent>/fc_matrices_456."
        ),
    )
    parser.add_argument(
        "--expected-parcels",
        type=int,
        default=EXPECTED_N_PARCELS,
        help=f"Expected number of timeseries columns (default: {EXPECTED_N_PARCELS}).",
    )
    parser.add_argument(
        "--raw-correlation-average",
        action="store_true",
        help="Average correlations directly instead of using Fisher-z averaging.",
    )
    return parser.parse_args()


def load_network_labels(atlas_path: Path, expected_n_parcels: int) -> np.ndarray:
    """Load and validate the parcel-level network labels (dseg.tsv)"""
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Atlas TSV not found: {atlas_path}")

    atlas_df = pd.read_csv(atlas_path, sep="\t")
    if "network_label" not in atlas_df.columns:
        raise ValueError("Atlas TSV must contain a 'network_label' column.")
    if len(atlas_df) != expected_n_parcels:
        raise ValueError(
            f"Atlas has {len(atlas_df)} rows; expected {expected_n_parcels}."
        )
    return atlas_df["network_label"].to_numpy()


def create_valid_parcel_mask(network_labels: np.ndarray) -> np.ndarray:
    """Return a mask for parcels assigned to a cortical network."""
    valid_mask = ~pd.isna(network_labels)
    print("\nAtlas label summary")
    print(f"  total parcels: {len(network_labels)}")
    print(f"  network-labelled parcels: {int(valid_mask.sum())}")
    print(f"  excluded from network summaries: {int((~valid_mask).sum())}")
    return valid_mask


def extract_sub_ses_id(file_path: str | Path) -> str:
    """Extract ``sub-*`` and ``ses-*`` entities from a BIDS-like filename."""
    name = Path(file_path).name
    match = re.search(r"(sub-[^_]+)_(ses-[^_]+)", name)
    if match is None:
        raise ValueError(
            "Could not extract a subject/session ID from filename: " f"{name}"
        )
    return f"{match.group(1)}_{match.group(2)}"


def collect_timeseries_files(timeseries_dir: Path, pattern: str) -> list[Path]:
    """Collect input timeseries files"""
    if not timeseries_dir.is_dir():
        raise NotADirectoryError(f"Timeseries directory not found: {timeseries_dir}")

    files = [Path(p) for p in sorted(glob(str(timeseries_dir / pattern), recursive=True))]
    if not files:
        raise FileNotFoundError(
            f"No timeseries matched {pattern!r} in {timeseries_dir}."
        )
    return files


def load_timeseries(file_path: Path, expected_n_parcels: int) -> np.ndarray:
    """Load one TSV as a timepoints-by-parcels numeric array."""
    ts_df = pd.read_csv(file_path, sep="\t")
    if ts_df.shape[1] != expected_n_parcels:
        raise ValueError(
            f"{file_path} has {ts_df.shape[1]} columns; "
            f"expected {expected_n_parcels}."
        )

    ts = ts_df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if ts.shape[0] < 2:
        raise ValueError(f"Need at least two timepoints: {file_path}")
    return ts


def compute_fc_matrix(timeseries: np.ndarray) -> np.ndarray:
    """Compute a parcel-by-parcel Pearson correlation matrix."""
    return np.corrcoef(timeseries, rowvar=False)


def mean_correlations(values: np.ndarray, use_fisher_z: bool) -> float:
    """Average finite correlation values, optionally in Fisher-z space."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan

    if use_fisher_z:
        values = np.clip(values, -0.999999, 0.999999)
        return float(np.tanh(np.mean(np.arctanh(values))))
    return float(np.mean(values))


def compute_mean_within_network(
    fc_valid: np.ndarray,
    network_labels_valid: np.ndarray,
    network: str,
    use_fisher_z: bool,
) -> float:
    """Compute mean off-diagonal connectivity within one network."""
    idx = np.flatnonzero(network_labels_valid == network)
    if idx.size < 2:
        return np.nan

    submatrix = fc_valid[np.ix_(idx, idx)]
    off_diagonal = submatrix[~np.eye(idx.size, dtype=bool)]
    return mean_correlations(off_diagonal, use_fisher_z=use_fisher_z)


def compute_mean_between_network(
    fc_valid: np.ndarray,
    network_labels_valid: np.ndarray,
    network_a: str,
    network_b: str,
    use_fisher_z: bool,
) -> float:
    """Compute mean connectivity between two distinct networks."""
    idx_a = np.flatnonzero(network_labels_valid == network_a)
    idx_b = np.flatnonzero(network_labels_valid == network_b)
    if idx_a.size == 0 or idx_b.size == 0:
        return np.nan

    submatrix = fc_valid[np.ix_(idx_a, idx_b)]
    return mean_correlations(submatrix.ravel(), use_fisher_z=use_fisher_z)


def average_corr_matrices(
    matrices: Sequence[np.ndarray],
    use_fisher_z: bool,
) -> np.ndarray:
    """Average same-sized correlation matrices and restore a unit diagonal."""
    if not matrices:
        raise ValueError("No correlation matrices were supplied for averaging.")

    stack = np.stack(matrices, axis=0).astype(float, copy=False)
    if use_fisher_z:
        stack = np.arctanh(np.clip(stack, -0.999999, 0.999999))
        average = np.tanh(np.nanmean(stack, axis=0))
    else:
        average = np.nanmean(stack, axis=0)

    np.fill_diagonal(average, 1.0)
    return average


def main() -> None:
    """Run the network-connectivity workflow."""
    args = parse_args()
    use_fisher_z = not args.raw_correlation_average

    matrix_dir = args.matrix_dir or args.output_csv.parent / f"fc_matrices_{args.expected_parcels}"
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    matrix_dir.mkdir(parents=True, exist_ok=True)

    network_labels = load_network_labels(args.atlas_tsv, args.expected_parcels)
    valid_mask = create_valid_parcel_mask(network_labels)
    network_labels_valid = network_labels[valid_mask]
    networks = sorted(pd.unique(network_labels_valid).tolist())

    print("\nNetworks retained")
    for network in networks:
        n_parcels = int(np.sum(network_labels_valid == network))
        print(f"  {network}: {n_parcels} parcels")

    files = collect_timeseries_files(args.timeseries_dir, args.pattern)
    print(f"\nProcessing {len(files)} timeseries files")

    results: list[dict[str, object]] = []
    matrices: list[np.ndarray] = []

    for file_path in files:
        scan_id = extract_sub_ses_id(file_path)
        print(f"  {scan_id}: {file_path.name}")

        timeseries = load_timeseries(file_path, args.expected_parcels)
        fc = compute_fc_matrix(timeseries)
        matrices.append(fc)
        np.save(matrix_dir / f"{scan_id}_fc_matrix.npy", fc)

        # Apply the atlas mask once; all subsequent indices refer to this
        # network-labelled parcel space.
        fc_valid = fc[np.ix_(valid_mask, valid_mask)]
        row: dict[str, object] = {"sub_ses": scan_id}

        for network in networks:
            row[network] = compute_mean_within_network(
                fc_valid,
                network_labels_valid,
                network,
                use_fisher_z,
            )

        for i, network_a in enumerate(networks):
            for network_b in networks[i + 1 :]:
                row[f"{network_a}_to_{network_b}"] = compute_mean_between_network(
                    fc_valid,
                    network_labels_valid,
                    network_a,
                    network_b,
                    use_fisher_z,
                )

        results.append(row)

    pd.DataFrame(results).to_csv(args.output_csv, index=False)
    np.save(
        matrix_dir / "avg_fc_matrix.npy",
        average_corr_matrices(matrices, use_fisher_z=use_fisher_z),
    )

    print("\nDone")
    print(f"  Fisher-z averaging: {use_fisher_z}")
    print(f"  summary CSV: {args.output_csv}")
    print(f"  matrix directory: {matrix_dir}")


if __name__ == "__main__":
    main()
