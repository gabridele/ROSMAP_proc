import nibabel as nb
import json
import os
import numpy as np
import datalad.api as dl
import sys
from pathlib import Path

def get_json(path):
    """
    Get the repetition time from the JSON sidecar file.
    
    Args:
        path (str): Path to the NIfTI file.
    returns:
        repetition_time (float): Repetition time from the JSON file.
    """

    print('path:', path)
    # get corresponding json file
    json_path = str(path).replace(".nii.gz", ".json")
    print(json_path)
    
    repetition_time = None
    
    try:
        # extract TR from json
        with open(json_path, 'r') as json_file:
            data = json.load(json_file)
            repetition_time = data.get('RepetitionTime')
            if repetition_time is not None:
                print(f"Repetition Time: {repetition_time}")
            else:
                print("Key 'repetitiontime' not found in JSON file.")
    except FileNotFoundError:
        print(f"JSON file not found: {json_path}")
    except json.JSONDecodeError:
        print(f"Error decoding JSON file: {json_path}")
    return repetition_time
# nake sure its passing right type

def strict_load(fname, subdset):
    """ Load an image that will write out the same as the input image.
    
    Args:
        fname (str): Path to the NIfTI file.
    Returns:
        strict_img (Nifti1Image): Loaded NIfTI image.
        zoom (float): Voxel size in the z-direction.
    """
    
    orig = nb.load(fname)
    strict_img = orig.__class__(np.array(orig.dataobj.get_unscaled()), orig.affine, orig.header)
    strict_img.header.set_slope_inter(orig.dataobj.slope, orig.dataobj.inter)
    header = strict_img.header
    zooms = header.get_zooms()[:3]
    return strict_img, zooms[2]

def set_tr(img, tr):
    """
    This function modifies the TR in the NIfTI header to match the value
    provided in the JSON sidecar file, due to mismatch found during bids validation.
    The TR is set to the value in seconds.

    Args:
        img (Nifti1Image): NIfTI image.
        tr (float): Repetition time to set in the header.
    Returns:
        img (Nifti1Image): NIfTI image with updated TR.
    """

    # get header
    header = img.header
    # get the current zooms up to 3 dimensions (fourth entry is for TR). And append new TR to it
    zooms = header.get_zooms()[:3] + (tr,)
    # set the new zooms
    header.set_zooms(zooms)
    return img


def main(path):
    
    subdset = Path(path).parts[0]
    print(subdset)
    repetition_time = get_json(path, subdset)
    img = strict_load(path, subdset)
    fixed_img = set_tr(img, repetition_time)
    fixed_img.to_filename(path)
    dl.save(path=path, dataset=subdset, message=f'changed TR in {path} header to match json sidecar')

if __name__ == "__main__":

# Check if script is executed with correct number of args
    if len(sys.argv) != 2:
        print("Correct syntax: python [this_script.py] [path]")
        sys.exit(1)
    
    path = str(sys.argv[1])

    main(path)