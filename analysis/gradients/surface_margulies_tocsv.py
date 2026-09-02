import nibabel as nib
import numpy as np
import pandas as pd

from neuromaps.datasets import fetch_annotation
from netneurotools.datasets import fetch_schaefer2018
from netneurotools.interface import vertices_to_parcels


# --------------------------------------------------
# Settings
# --------------------------------------------------

N_PARCELS = 400
N_NETWORKS = 7
N_GRADIENTS = 5

ATLAS_KEY = f"{N_PARCELS}Parcels{N_NETWORKS}Networks"
OUTPUT_FILE = (
    f"margulies_gradients1-{N_GRADIENTS}_"
    f"schaefer{N_PARCELS}_{N_NETWORKS}net.csv"
)


# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def load_surface_data(files):
    """Load left- and right-hemisphere GIFTI surface data."""

    return tuple(
        np.asarray(nib.load(str(filename)).agg_data()).squeeze()
        for filename in files
    )


def clean_label(label):
    """Convert byte labels to ordinary strings."""

    if isinstance(label, bytes):
        return label.decode("utf-8")

    return str(label)


# --------------------------------------------------
# Download the Schaefer atlas
# --------------------------------------------------

schaefer = fetch_schaefer2018(version="fslr32k")

if ATLAS_KEY not in schaefer:
    raise KeyError(
        f"{ATLAS_KEY!r} was not found. "
        f"Available atlases include: {list(schaefer.keys())}"
    )

parcellation_files = schaefer[ATLAS_KEY]


# --------------------------------------------------
# Parcellate gradients 1 through 5
# --------------------------------------------------

output = None
reference_labels = None

for gradient_number in range(1, N_GRADIENTS + 1):

    description = f"fcgradient{gradient_number:02d}"

    print(f"Processing {description}...")

    # Download left- and right-hemisphere gradient maps
    gradient_files = fetch_annotation(
        source="margulies2016",
        desc=description,
        space="fsLR",
        den="32k",
        return_single=True
    )

    gradient_vertex_data = load_surface_data(gradient_files)

    # Average vertex values within each Schaefer parcel
    parcel_values, parcel_keys, parcel_labels = vertices_to_parcels(
        vert_data=gradient_vertex_data,
        parc_file=parcellation_files,
        hemi="both"
    )

    # Combine left and right hemispheres
    values = np.concatenate(
        [np.asarray(values).reshape(-1) for values in parcel_values]
    )

    labels = [
        clean_label(label)
        for hemisphere_labels in parcel_labels
        for label in hemisphere_labels
    ]

    # Build metadata during the first iteration
    if output is None:

        left_count = len(parcel_values[0])
        right_count = len(parcel_values[1])

        output = pd.DataFrame({
            "parcel_id": np.arange(1, len(values) + 1),
            "region": labels,
            "hemisphere": (
                ["left"] * left_count +
                ["right"] * right_count
            )
        })

        reference_labels = labels

    # Confirm every gradient uses the same parcel ordering
    elif labels != reference_labels:
        raise RuntimeError(
            f"Parcel ordering changed while processing {description}."
        )

    output[f"G{gradient_number}"] = values


# --------------------------------------------------
# Validate and export
# --------------------------------------------------

if len(output) != N_PARCELS:
    raise ValueError(
        f"Expected {N_PARCELS} parcels, but obtained {len(output)}."
    )

output.to_csv(OUTPUT_FILE, index=False)

print(f"\nSaved: {OUTPUT_FILE}")
print(output.head())
print(output.shape)