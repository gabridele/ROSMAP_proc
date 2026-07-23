#!/usr/bin/env python3

import os
import glob
import argparse
import numpy as np
import nibabel as nib
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def build_hyperslab(files, i=20, k=50, t=10):
    """
    Build a hyperslab matrix.

    Output shape:
        n_images x n_voxels_along_j

    Each row:
        one EPI image

    Each column:
        voxel coordinate j
    """

    lines = []
    filenames = []

    for f in files:
        img = nib.load(f)
        data = img.dataobj

        if len(img.shape) == 4:
            if t >= img.shape[3]:
                raise IndexError(
                    f"Requested t={t}, but file has only {img.shape[3]} timepoints: {f}"
                )

            if i >= img.shape[0] or k >= img.shape[2]:
                raise IndexError(
                    f"Requested i={i}, k={k}, but image shape is {img.shape}: {f}"
                )

            line = np.asarray(data[i, :, k, t])

        elif len(img.shape) == 3:
            if i >= img.shape[0] or k >= img.shape[2]:
                raise IndexError(
                    f"Requested i={i}, k={k}, but image shape is {img.shape}: {f}"
                )

            line = np.asarray(data[i, :, k])

        else:
            raise ValueError(f"Unsupported image dimensions {img.shape}: {f}")

        lines.append(line)
        filenames.append(os.path.basename(f))

    H = np.vstack(lines)

    return H, filenames


def make_customdata(filenames, n_j):
    """
    Hover metadata for Plotly heatmaps.

    customdata[row, col] contains:
        image index
        filename
        voxel j
    """

    n_images = len(filenames)

    image_indices = np.arange(n_images)
    j_coords = np.arange(n_j)

    customdata = np.empty((n_images, n_j, 3), dtype=object)

    customdata[:, :, 0] = np.repeat(image_indices[:, None], n_j, axis=1)
    customdata[:, :, 1] = np.repeat(
        np.array(filenames, dtype=object)[:, None],
        n_j,
        axis=1
    )
    customdata[:, :, 2] = np.repeat(j_coords[None, :], n_images, axis=0)

    return customdata


def global_percentiles(H, low_pct=2.5, high_pct=97.5):
    """
    Combine every value in H into one vector and compute percentiles.
    """

    values = np.asarray(H).ravel()
    values = values[np.isfinite(values)]

    if values.size == 0:
        raise ValueError("No finite values found for percentile calculation.")

    low = np.percentile(values, low_pct)
    high = np.percentile(values, high_pct)

    return low, high


def align_files_by_name(files_1, files_2):
    """
    Align files by simplified basename.

    Edit file_key() if your filenames differ in another way.
    """

    def file_key(path):
        name = os.path.basename(path)

        name = name.replace("_v1.3", "")
        name = name.replace("_v13", "")
        name = name.replace("_v1", "")
        name = name.replace("_BBR", "")
        name = name.replace("_noBBR", "")
        name = name.replace("_nobbr", "")

        return name

    dict_1 = {file_key(f): f for f in files_1}
    dict_2 = {file_key(f): f for f in files_2}

    common_keys = sorted(set(dict_1) & set(dict_2))

    missing_from_2 = sorted(set(dict_1) - set(dict_2))
    missing_from_1 = sorted(set(dict_2) - set(dict_1))

    if len(common_keys) == 0:
        raise ValueError("No matching files found between the two datasets.")

    if missing_from_2:
        print(f"Warning: {len(missing_from_2)} dataset 1 files missing from dataset 2")

    if missing_from_1:
        print(f"Warning: {len(missing_from_1)} dataset 2 files missing from dataset 1")

    files_1_aligned = [dict_1[key] for key in common_keys]
    files_2_aligned = [dict_2[key] for key in common_keys]

    return files_1_aligned, files_2_aligned


def add_percentile_contours(
    fig,
    H_target,
    low,
    high,
    row,
    col,
    label,
    line_color="black"
):
    """
    Overlay dotted contour lines on H_target at low and high intensity values.

    These values are usually computed from the opposite dataset.
    """

    fig.add_trace(
        go.Contour(
            z=H_target,
            contours=dict(
                start=low,
                end=low,
                size=1,
                coloring="none",
                showlabels=False
            ),
            line=dict(
                color=line_color,
                width=1,
                dash="dot"
            ),
            showscale=False,
            name=f"{label} 2.5%",
            hovertemplate=(
                f"{label} lower percentile<br>"
                f"Intensity: {low:.4f}<br>"
                "Image index: %{y}<br>"
                "Voxel j: %{x}<extra></extra>"
            )
        ),
        row=row,
        col=col
    )

    fig.add_trace(
        go.Contour(
            z=H_target,
            contours=dict(
                start=high,
                end=high,
                size=1,
                coloring="none",
                showlabels=False
            ),
            line=dict(
                color=line_color,
                width=1,
                dash="dot"
            ),
            showscale=False,
            name=f"{label} 97.5%",
            hovertemplate=(
                f"{label} upper percentile<br>"
                f"Intensity: {high:.4f}<br>"
                "Image index: %{y}<br>"
                "Voxel j: %{x}<extra></extra>"
            )
        ),
        row=row,
        col=col
    )

def normalize_to_range(H, vmin, vmax):
    Hn = (H.astype(float) - vmin) / (vmax - vmin)
    Hn = np.clip(Hn, 0, 1)
    Hn[~np.isfinite(Hn)] = 0
    return Hn

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create two interactive hyperslab plots with global 2.5% and 97.5% "
            "percentile boundaries from each dataset overlaid on the opposite dataset."
        )
    )

    parser.add_argument(
        "--dataset1_dir",
        required=True,
        help="Folder containing dataset 1, for example BBR."
    )

    parser.add_argument(
        "--dataset2_dir",
        required=True,
        help="Folder containing dataset 2, for example no BBR."
    )

    parser.add_argument(
        "--dataset1_label",
        default="BBR",
        help="Label for dataset 1."
    )

    parser.add_argument(
        "--dataset2_label",
        default="no BBR",
        help="Label for dataset 2."
    )

    parser.add_argument(
        "--output",
        default="hyperslab_two_datasets_percentiles.html",
        help="Output interactive HTML file."
    )

    parser.add_argument("--i", type=int, default=20)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--t", type=int, default=10)

    parser.add_argument(
        "--pattern",
        default="*.nii*",
        help="Glob pattern for NIfTI files."
    )

    parser.add_argument(
        "--cmap",
        default="YlGnBu",
        help="Plotly colorscale for both heatmaps."
    )

    parser.add_argument(
        "--low_pct",
        type=float,
        default=2.5,
        help="Lower percentile."
    )

    parser.add_argument(
        "--high_pct",
        type=float,
        default=97.5,
        help="Upper percentile."
    )

    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Fixed color scale minimum. If omitted, combined low_pct percentile is used."
    )

    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Fixed color scale maximum. If omitted, combined high_pct percentile is used."
    )

    parser.add_argument(
        "--align_by_name",
        action="store_true",
        help="Align files by simplified filename instead of sorted order."
    )

    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--width", type=int, default=1600)

    args = parser.parse_args()

    files_1 = sorted(glob.glob(os.path.join(args.dataset1_dir, args.pattern)))
    files_2 = sorted(glob.glob(os.path.join(args.dataset2_dir, args.pattern)))

    if len(files_1) == 0:
        raise FileNotFoundError(f"No files found in dataset 1 folder: {args.dataset1_dir}")

    if len(files_2) == 0:
        raise FileNotFoundError(f"No files found in dataset 2 folder: {args.dataset2_dir}")

    if args.align_by_name:
        files_1, files_2 = align_files_by_name(files_1, files_2)
    else:
        if len(files_1) != len(files_2):
            raise ValueError(
                "Different number of files.\n"
                f"{args.dataset1_label}: {len(files_1)} files\n"
                f"{args.dataset2_label}: {len(files_2)} files\n\n"
                "Either make sure folders contain matched files, or rerun with:\n"
                "    --align_by_name"
            )

    print(f"Using {len(files_1)} files from each dataset")

    print(f"Building hyperslab for {args.dataset1_label}...")
    H1, filenames_1 = build_hyperslab(
        files_1,
        i=args.i,
        k=args.k,
        t=args.t
    )

    print(f"Building hyperslab for {args.dataset2_label}...")
    H2, filenames_2 = build_hyperslab(
        files_2,
        i=args.i,
        k=args.k,
        t=args.t
    )

    if H1.shape != H2.shape:
        raise ValueError(
            "Hyperslab shapes do not match.\n"
            f"{args.dataset1_label}: {H1.shape}\n"
            f"{args.dataset2_label}: {H2.shape}"
        )

    n_images, n_j = H1.shape

    print(f"Hyperslab shape: {H1.shape}")

    print("Computing global percentiles from all hyperslab values...")

    H1_low, H1_high = global_percentiles(
        H1,
        low_pct=args.low_pct,
        high_pct=args.high_pct
    )

    H2_low, H2_high = global_percentiles(
        H2,
        low_pct=args.low_pct,
        high_pct=args.high_pct
    )

    print(
        f"{args.dataset1_label}: "
        f"{args.low_pct}% = {H1_low:.4f}, "
        f"{args.high_pct}% = {H1_high:.4f}"
    )

    print(
        f"{args.dataset2_label}: "
        f"{args.low_pct}% = {H2_low:.4f}, "
        f"{args.high_pct}% = {H2_high:.4f}"
    )

    if args.vmin is None or args.vmax is None:
        combined = np.concatenate([H1.ravel(), H2.ravel()])
        combined = combined[np.isfinite(combined)]

        auto_vmin = np.percentile(combined, args.low_pct)
        auto_vmax = np.percentile(combined, args.high_pct)

        vmin = auto_vmin if args.vmin is None else args.vmin
        vmax = auto_vmax if args.vmax is None else args.vmax
    else:
        vmin = args.vmin
        vmax = args.vmax

    print(f"Shared color scale: vmin={vmin:.4f}, vmax={vmax:.4f}")

    customdata_1 = make_customdata(filenames_1, n_j)
    customdata_2 = make_customdata(filenames_2, n_j)

    # # ----- RGB overlap hyperslab -----
    # H1n = normalize_to_range(H1, vmin, vmax)
    # H2n = normalize_to_range(H2, vmin, vmax)

    # # Dataset 1, e.g. BBR, is green
    # # Dataset 2, e.g. no-BBR, is magenta
    # rgb = np.zeros((H1.shape[0], H1.shape[1], 3), dtype=float)
    # rgb[:, :, 0] = H2n  # red channel: dataset 2
    # rgb[:, :, 1] = H1n  # green channel: dataset 1
    # rgb[:, :, 2] = H2n  # blue channel: dataset 2

    # rgb_uint8 = (rgb * 255).astype(np.uint8)
    # -------------------------------

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=(
            f"{args.dataset1_label} hyperslab + {args.dataset2_label} percentiles",
            f"{args.dataset2_label} hyperslab + {args.dataset1_label} percentiles",
        ),
        shared_yaxes=True,
        horizontal_spacing=0.06,
        column_widths=[0.5, 0.5]
    )

    fig.add_trace(
        go.Heatmap(
            z=H1,
            colorscale=args.cmap,
            zmin=vmin,
            zmax=vmax,
            customdata=customdata_1,
            hovertemplate=(
                f"Dataset: {args.dataset1_label}<br>"
                "Image index: %{customdata[0]}<br>"
                "Filename: %{customdata[1]}<br>"
                "Voxel j: %{customdata[2]}<br>"
                "Intensity: %{z}<extra></extra>"
            ),
            colorbar=dict(
                title="Intensity",
                x=1.03,
                y=0.5,
                len=0.7,
                thickness=14
            )
        ),
        row=1,
        col=1
    )

    fig.add_trace(
        go.Heatmap(
            z=H2,
            colorscale=args.cmap,
            zmin=vmin,
            zmax=vmax,
            customdata=customdata_2,
            hovertemplate=(
                f"Dataset: {args.dataset2_label}<br>"
                "Image index: %{customdata[0]}<br>"
                "Filename: %{customdata[1]}<br>"
                "Voxel j: %{customdata[2]}<br>"
                "Intensity: %{z}<extra></extra>"
            ),
            showscale=False
        ),
        row=1,
        col=2
    )

    # Overlay dataset 2 percentile contours onto dataset 1 heatmap
    add_percentile_contours(
        fig=fig,
        H_target=H1,
        low=H2_low,
        high=H2_high,
        row=1,
        col=1,
        label=f"{args.dataset2_label} percentiles",
        line_color="black"
    )

    # Overlay dataset 1 percentile contours onto dataset 2 heatmap
    add_percentile_contours(
        fig=fig,
        H_target=H2,
        low=H1_low,
        high=H1_high,
        row=1,
        col=2,
        label=f"{args.dataset1_label} percentiles",
        line_color="black"
    )
    
    # fig.add_trace(
    #     go.Image(
    #         z=rgb_uint8,
    #         customdata=customdata_1,
    #         hovertemplate=(
    #             "Overlap view<br>"
    #             "Image index: %{customdata[0]}<br>"
    #             "Filename: %{customdata[1]}<br>"
    #             "Voxel j: %{customdata[2]}<extra></extra>"
    #         )
    #     ),
    #     row=1,
    #     col=3
    # )
    fig.update_layout(
        title=(
            f"Hyperslab percentile comparison: "
            f"{args.dataset1_label} vs {args.dataset2_label} "
            f"(i={args.i}, k={args.k}, t={args.t})"
        ),
        width=args.width,
        height=args.height,
        margin=dict(l=80, r=230, t=100, b=80),
        hovermode="closest",
        hoverlabel=dict(
            bgcolor="white",
            font_size=11,
            namelength=-1
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.03,
            xanchor="center",
            x=0.5
        )
    )

    fig.update_xaxes(title_text="Voxel coordinate j", row=1, col=1)
    fig.update_xaxes(title_text="Voxel coordinate j", row=1, col=2)
    #fig.update_xaxes(title_text="Voxel coordinate j", row=1, col=3)
    fig.update_yaxes(title_text="EPI image index", row=1, col=1)

    fig.write_html(args.output)

    print(f"Saved interactive plot to: {args.output}")


if __name__ == "__main__":
    main()