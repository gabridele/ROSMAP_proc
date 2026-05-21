"""
Compute mean within-network functional connectivity for 4S456 atlas timeseries.

This script:
1. Loads parcel network labels from the atlas TSV.
2. Removes parcels with missing network labels.
3. Collects timeseries files from folders A, B, and C.
   - A and B are both included fully.
   - C is included only when files are not already present in A or B.
4. Computes subject-level FC matrices.
5. Computes mean within-network connectivity for each network.
6. Saves results to a CSV.

Output:
    One row per sub_ses
    First column = sub_ses ID.
    Remaining columns = mean within-network connectivity per network.
    NB: Some networks may get Nan mean within connectivity because there was no coverage for all the parcels in that network. This happens specifically for the limbic network, mostly due to dropout.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from glob import glob

# =============================================================================
# Configuration
# =============================================================================

ATLAS_PATH = Path(
    "/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/"
    "atlas-4S456Parcels/atlas-4S456Parcels_dseg.tsv"
)

# supposedly all same input folder

OUTPUT_CSV = Path("/Users/ga0034de/github_dir/ROSMAP_proc/analysis/april26/mean_connectivity_190526.csv")

EXPECTED_N_PARCELS = 456

# Set to False if you want raw, untransformed correlations.
USE_FISHER_Z = True


# =============================================================================
# Loading atlas labels
# =============================================================================

def load_network_labels(atlas_path: Path) -> np.ndarray:
    """
    Load parcel-level network labels from the atlas TSV.

    Parameters
    ----------
    atlas_path : Path
        Path to atlas TSV file. Must contain a column named 'network_label'.

    Returns
    -------
    np.ndarray
        Array of network labels, one per parcel.

    Raises
    ------
    ValueError
        If the file does not contain the required 'network_label' column.
    """
    atlas_df = pd.read_csv(atlas_path, sep="\t")

    if "network_label" not in atlas_df.columns:
        raise ValueError("Atlas file must contain a 'network_label' column.")

    return atlas_df["network_label"].to_numpy()


def create_valid_parcel_mask(network_labels: np.ndarray) -> np.ndarray:
    """
    Create mask identifying parcels with valid network labels.

    Parcels with NaN network labels are excluded from both the timeseries
    and the FC matrix.

    Parameters
    ----------
    network_labels : np.ndarray
        Full array of parcel network labels.

    Returns
    -------
    np.ndarray
        Boolean mask where True means the parcel has a valid label.
    """
    valid_mask = ~pd.isna(network_labels)

    print("\n--- Atlas label summary ---")
    print(f"Total parcels: {len(network_labels)}")
    print(f"Valid labeled parcels: {valid_mask.sum()}")
    print(f"Removed NaN parcels: {(~valid_mask).sum()}")

    removed_idx = np.where(~valid_mask)[0]
    if len(removed_idx) > 0:
        print("Removed parcel indices, first 20:", removed_idx[:20])

    return valid_mask


# =============================================================================
# File handling
# =============================================================================

def extract_sub_ses_id(file_path: Path) -> str:
    """
    Extract subject/session ID from a timeseries filename.

    Example
    -------
    Input:
        sub-001_ses-01_task-rest_timeseries.csv

    Output:
        sub-001_ses-01

    Parameters
    ----------
    file_path : Path
        Path to timeseries file.

    Returns
    -------
    str
        Subject/session identifier.
    """
    # get filename from abs path
    file_path = file_path.split("/")[-1]
    parts = file_path.split("_")
    return "_".join(parts[:2])


# def collect_timeseries_files(
#     folder_a: Path,
#     folder_b: Path,
#     folder_c: Path,
# ) -> list[Path]:
#     """
#     Collect timeseries files from folders A, B, and C.

#     Rules
#     -----
#     - Include all CSV files from folder A.
#     - Include all CSV files from folder B.
#     - Include CSV files from folder C only if their subject/session ID
#       is not already present in A or B.

#     Deduplication is based on subject/session ID, not exact filename. This is
#     safer for neuroimaging workflows where filenames can differ after the
#     second underscore.

#     Parameters
#     ----------
#     folder_a : Path
#         First input folder.
#     folder_b : Path
#         Second input folder.
#     folder_c : Path
#         Third input folder that may contain duplicates.

#     Returns
#     -------
#     list[Path]
#         Combined list of timeseries files to process.
#     """
#     files_a = sorted(folder_a.glob("sub-*_ses-*/*nifti/sub-*/ses-*/func/*space-MNI152NLin6Asym_seg-4S456Parcels_stat-mean_timeseries.tsv"))
#     files_b = sorted(folder_b.glob("sub-*_ses-*/*nifti/sub-*/ses-*/func/*space-MNI152NLin6Asym_seg-4S456Parcels_stat-mean_timeseries.tsv"))
#     files_c_all = sorted(folder_c.glob("sub-*_ses-*/*nifti/sub-*/ses-*/func/*space-MNI152NLin6Asym_seg-4S456Parcels_stat-mean_timeseries.tsv"))
    
#     ids_ab = {extract_sub_ses_id(f) for f in files_a + files_b}

#     files_c = [
#         f for f in files_c_all
#         if extract_sub_ses_id(f) not in ids_ab
#     ]

#     skipped_c = [
#         f for f in files_c_all
#         if extract_sub_ses_id(f) in ids_ab
#     ]

#     print("\n--- File summary ---")
#     print(f"Folder A files: {len(files_a)}")
#     print(f"Folder B files: {len(files_b)}")
#     print(f"Folder C files total: {len(files_c_all)}")
#     print(f"Folder C files included: {len(files_c)}")
#     print(f"Folder C files skipped as duplicates: {len(skipped_c)}")
#     print(f"Total files to process: {len(files_a) + len(files_b) + len(files_c)}")

#     if skipped_c:
#         print("\nSkipped duplicate files from C, first 10:")
#         for f in skipped_c[:10]:
#             print(f"  {f.name}")

#     return files_a + files_b + files_c


# =============================================================================
# FC computation
# =============================================================================

def load_and_filter_timeseries(
    file_path: Path,
    valid_mask: np.ndarray,
    expected_n_parcels: int,
) -> np.ndarray:
    """
    Load a timeseries CSV and remove unlabeled parcels.

    Parameters
    ----------
    file_path : Path
        Path to timeseries CSV.
    valid_mask : np.ndarray
        Boolean mask identifying valid labeled parcels.
    expected_n_parcels : int
        Expected number of columns before filtering.

    Returns
    -------
    np.ndarray
        Filtered timeseries array with shape:
        timepoints x valid_parcels.

    Raises
    ------
    ValueError
        If the input file does not have the expected number of parcels.
    """
    ts = pd.read_csv(file_path, sep="\t").to_numpy()
    # shape of pandas dataframe is (timepoints, parcels)
    print(f"Shape of {file_path}: {ts.shape}")
    if ts.shape[1] != expected_n_parcels:
        raise ValueError(
            f"{file_path} has {ts.shape[1]} columns, "
            f"expected {expected_n_parcels}."
        )

    ts = ts[:, valid_mask]

    return ts


def compute_fc_matrix(timeseries: np.ndarray) -> np.ndarray:
    """
    Compute parcel-by-parcel functional connectivity matrix.

    Parameters
    ----------
    timeseries : np.ndarray
        Timeseries array with shape timepoints x parcels.

    Returns
    -------
    np.ndarray
        Pearson correlation matrix with shape parcels x parcels.
    """
    return np.corrcoef(timeseries, rowvar=False)


def compute_mean_within_network(
    fc: np.ndarray,
    network_labels: np.ndarray,
    network: str,
    use_fisher_z: bool = True,
) -> float:
    """
    Compute mean within-network connectivity for one network.

    Parameters
    ----------
    fc : np.ndarray
        Functional connectivity matrix.
    network_labels : np.ndarray
        Network labels after removing NaN parcels.
    network : str
        Network label to summarize.
    use_fisher_z : bool, optional
        If True, apply Fisher z transform before averaging.
        If False, average raw correlations.

    Returns
    -------
    float
        Mean within-network connectivity.

        If use_fisher_z=True, this returns the mean after applying the Fisher z transform and then back-transforming to correlation space.
        If use_fisher_z=False, this returns the mean raw correlation.
    """
    idx = np.where(network_labels == network)[0]

    if len(idx) < 2:
        return np.nan

    submat = fc[np.ix_(idx, idx)]

    # Keep only off-diagonal values.
    # Diagonal values are self-correlations and always equal 1.
    off_diag_mask = ~np.eye(len(idx), dtype=bool)
    edge_values = submat[off_diag_mask]
    edge_values = pd.to_numeric(edge_values, errors="coerce")
    if use_fisher_z:
        # Avoid infinite values if any correlations are exactly -1 or 1.
        edge_values = np.clip(edge_values, -0.999999, 0.999999)

        # Fisher z transform.
        # z-transform is used because you are averaging correlation coefficients. 
        # It is statistically cleaner than averaging raw correlations directly.

        # To turn this off, set USE_FISHER_Z = False at the top of the script.
        edge_values = np.arctanh(edge_values)
        return np.tanh(np.nanmean(edge_values))
    return np.nanmean(edge_values)

#### compute between-network connectivity
def compute_mean_between_network(
    fc: np.ndarray,
    network_labels: np.ndarray,
    network_a: str,
    network_b: str,
    use_fisher_z: bool = True,
) -> float:

    """
    Compute mean between-network connectivity for two networks.

    Parameters
    ----------
    fc : np.ndarray
        Functional connectivity matrix.
    network_labels : np.ndarray
        Network labels after removing NaN parcels.
    network_a : str
        First network label.
    network_b : str
        Second network label.
    use_fisher_z : bool, optional
        If True, apply Fisher z transform before averaging.
        If False, average raw correlations.

    Returns
    -------
    float
        Mean between-network connectivity.

        If use_fisher_z=True, this returns the mean after applying the Fisher z transform and then back-transforming to correlation space.
        If use_fisher_z=False, this returns the mean raw correlation.
    """
    idx_a = np.where(network_labels == network_a)[0]
    idx_b = np.where(network_labels == network_b)[0]

    if len(idx_a) == 0 or len(idx_b) == 0:
        return np.nan

    submat = fc[np.ix_(idx_a, idx_b)]
    edge_values = submat.flatten()
    edge_values = pd.to_numeric(edge_values, errors="coerce")
    
    if use_fisher_z:
        edge_values = np.clip(edge_values, -0.999999, 0.999999)
        edge_values = np.arctanh(edge_values)
        return np.tanh(np.nanmean(edge_values))
    
    return np.nanmean(edge_values)

# =============================================================================
# Main workflow
# =============================================================================

def main() -> None:
    """
    Run the full within-network connectivity pipeline.
    """
    network_labels = load_network_labels(ATLAS_PATH)

    if len(network_labels) != EXPECTED_N_PARCELS:
        raise ValueError(
            f"Atlas has {len(network_labels)} parcels, "
            f"expected {EXPECTED_N_PARCELS}."
        )

    valid_mask = create_valid_parcel_mask(network_labels)

    network_labels_valid = network_labels[valid_mask]
    networks = sorted(pd.unique(network_labels_valid))

    print("\n--- Networks retained ---")
    for network in networks:
        n_parcels = np.sum(network_labels_valid == network)
        print(f"{network}: {n_parcels} parcels")

    files = glob("/Users/ga0034de/Desktop/timeseries_1905/*.tsv")

    results = []

    print("\n--- Processing timeseries files ---")

    for file_path in files:
        print(f"Processing: {file_path}")

        ts = load_and_filter_timeseries(
            file_path=file_path,
            valid_mask=valid_mask,
            expected_n_parcels=EXPECTED_N_PARCELS,
        )

        if ts.shape[1] != len(network_labels_valid):
            raise ValueError(
                f"After filtering, {file_path} has {ts.shape[1]} parcels, "
                f"but atlas labels have {len(network_labels_valid)} parcels."
            )

        fc = compute_fc_matrix(ts)

        print(f"  FC shape after filtering: {fc.shape}")

        row = {
            "timeseries": extract_sub_ses_id(file_path)
        }

        for network in networks:
            row[network] = compute_mean_within_network(
                fc=fc,
                network_labels=network_labels_valid,
                network=network,
                use_fisher_z=USE_FISHER_Z,
            )
        for i in range(len(networks)):
            for j in range(i+1, len(networks)):
                network_a = networks[i]
                network_b = networks[j]
                row[f"{network_a}_to_{network_b}"] = compute_mean_between_network(
                    fc=fc,
                    network_labels=network_labels_valid,
                    network_a=network_a,
                    network_b=network_b,
                    use_fisher_z=USE_FISHER_Z,
                )
        
        results.append(row)

    output_df = pd.DataFrame(results)
    output_df.to_csv(OUTPUT_CSV, index=False)

    print("\n--- Done ---")
    print(f"Fisher z used: {USE_FISHER_Z}")
    print(f"Saved output to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()