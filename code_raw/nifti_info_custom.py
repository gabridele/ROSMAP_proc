"""
Modifying this function taken from CuBIDS repo to tailor it to a datalad subdataset structure, i.e., where each subject is a subdataset within a superdataset.
The original function is located at: https://github.com/PennLINC/CuBIDS

"""
import datalad.api as dl
from pathlib import Path
import nibabel as nb
import numpy as np
import json
import subprocess

def img_to_new_ext(img_path, new_ext):
    """Convert an image file path to a new extension.

    Parameters
    ----------
    img_path : str
        The file path of the image to be converted.
    new_ext : str
        The new extension to be applied to the image file path.

    Returns
    -------
    str
        The file path with the new extension applied.

    Examples
    --------
    >>> img_to_new_ext('/path/to/file_image.nii.gz', '.tsv')
    '/path/to/file_events.tsv'

    >>> img_to_new_ext('/path/to/file_image.nii.gz', '.tsv.gz')
    '/path/to/file_physio.tsv.gz'

    >>> img_to_new_ext('/path/to/file_image.nii.gz', '.json')
    '/path/to/file_image.json'

    Notes
    -----
    The hardcoded suffix associated with each extension may not be comprehensive.
    BIDS has been extended a lot in recent years.
    """
    # handle .tsv edge case
    if new_ext == ".tsv":
        # take out suffix
        return img_path.rpartition("_")[0] + "_events" + new_ext
    elif new_ext == ".tsv.gz":
        return img_path.rpartition("_")[0] + "_physio" + new_ext
    else:
        return img_path.replace(".nii.gz", "").replace(".nii", "") + new_ext

def unlock_dataset(path):
    """
    Unlock a dataset using datalad.

    Parameters
    ----------
    path : str
        Path to the dataset to unlock.

    Raises
    ------
    Exception
        If there is an error unlocking the dataset.
    """
    for sub in Path(path).rglob("sub-*"):
        for json in Path(sub).rglob(f"ses-*/*/*.json"):
            try:
                dl.unlock(dataset=sub, path=str(json))
            except Exception:
                print("Error unlocking file: ", json)

def add_nifti_info(path):
        """
        Add information from NIfTI files to their corresponding JSON sidecars.

        This method processes all NIfTI files in the BIDS directory specified by `self.path`.
        It extracts relevant metadata from each NIfTI file and updates the corresponding JSON
        sidecar files with this information. If `self.force_unlock` is set, it will unlock
        the dataset using `datalad` before processing the files.

        Metadata added to the JSON sidecars includes:
        - Obliquity
        - Voxel sizes (dimensions 1, 2, and 3)
        - Matrix dimensions (sizes of dimensions 1, 2, and 3)
        - Number of volumes (for 4D images)
        - Image orientation

        If `self.use_datalad` is set, the changes will be saved using `datalad`.

        Raises
        ------
        Exception
            If there is an error loading a NIfTI file or parsing a JSON sidecar file.

        Notes
        -----
        - This method assumes that the NIfTI files are organized in a BIDS-compliant
            directory structure.
        - The method will skip any files in hidden directories (directories starting with a dot).
        - If a JSON sidecar file does not exist for a NIfTI file, it will be skipped.
        """

        # loop through all niftis in the bids dir
        for sub in Path(path).rglob("sub-*"):
            for path_ in Path(sub).rglob("ses-*/*/*"):
                print("Processing: ", path_)
                print("Sub: ", sub)
                
                # ignore all dot directories
                if "/." in str(path_):
                    continue

                if str(path_).endswith(".nii") or str(path_).endswith(".nii.gz"):
                    try:
                        img = nb.load(str(path_))
                    except Exception:
                        print("Empty Nifti File: ", str(path_))
                        continue

                    # get important info from niftis
                    obliquity = np.any(nb.affines.obliquity(img.affine) > 1e-4)
                    voxel_sizes = img.header.get_zooms()
                    matrix_dims = img.shape
                    # add nifti info to corresponding sidecars​
                    sidecar = img_to_new_ext(str(path_), ".json")
                    print("Sidecar: ", sidecar)
                    if Path(sidecar).exists():
                        try:
                            with open(sidecar) as f:
                                data = json.load(f)
                        except Exception:
                            print("Error parsing this sidecar: ", sidecar)

                        if "Obliquity" not in data.keys():
                            data["Obliquity"] = str(obliquity)
                        if "VoxelSizeDim1" not in data.keys():
                            data["VoxelSizeDim1"] = float(voxel_sizes[0])
                        if "VoxelSizeDim2" not in data.keys():
                            data["VoxelSizeDim2"] = float(voxel_sizes[1])
                        if "VoxelSizeDim3" not in data.keys():
                            data["VoxelSizeDim3"] = float(voxel_sizes[2])
                        if "Dim1Size" not in data.keys():
                            data["Dim1Size"] = matrix_dims[0]
                        if "Dim2Size" not in data.keys():
                            data["Dim2Size"] = matrix_dims[1]
                        if "Dim3Size" not in data.keys():
                            data["Dim3Size"] = matrix_dims[2]
                        if "NumVolumes" not in data.keys():
                            if img.ndim == 4:
                                data["NumVolumes"] = matrix_dims[3]
                            elif img.ndim == 3:
                                data["NumVolumes"] = 1
                        if "ImageOrientation" not in data.keys():
                            orient = nb.orientations.aff2axcodes(img.affine)
                            orient = [str(orientation) for orientation in orient]
                            joined = "".join(orient) + "+"
                            data["ImageOrientation"] = joined
                        print("Data: ", data)
                        git_log1 = subprocess.run(["git", "log", "-1", "--pretty=format:%H %s"], capture_output=True, text=True, cwd=sub)
                        latest_commit_msg = git_log1.stdout.strip()
                        print("Latest commit: ", latest_commit_msg)

                        if "Added nifti info to sidecars" not in latest_commit_msg:
                            try:
                                dl.unlock(dataset=sub, path=str(sidecar))
                            except Exception:
                                print("Error unlocking file: ", sidecar)
                        else:
                            continue

                        with open(sidecar, "w") as file:
                            json.dump(data, file, indent=4)

            if Path(sub).is_dir():
                try:
                    dl.save(dataset=sub, message="Added nifti info to sidecars")
                except Exception:
                    print("Error saving: ", sub)

add_nifti_info("raw")
