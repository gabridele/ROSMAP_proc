import os
from glob import glob
import pandas
from pathlib import Path
import sys
import zipfile
import numpy as np
import nibabel as nib
sys.path.insert(0,'./')
import datalad.api as dl
import json
import nibabel.processing as nli
"""
script taken and adapted from https://github.com/PennLINC/cros-map/tree/main/bids
"""

def convert_zipped_fmri(fpath, tmpdir, outfile):
    """
    Convert a zipped fMRI file to NIfTI format.
    Args:
        fpath (str): Path to the zipped fMRI file.
        tmpdir (str): Temporary directory for extraction.
        outfile (str): Output path for the converted NIfTI file.
    Returns:
        str: Path to the converted NIfTI file or 'error' if conversion fails.
    """
    # extract and concatenate
    try:
        with zipfile.ZipFile(fpath, 'r') as zip_ref:
            zip_ref.extractall(tmpdir)
    except:
        return 'error'
    fls = sorted(glob(os.path.join(tmpdir,'*.img')))
    img = nib.concat_images(fls,check_affines=False,axis=3)
    img = nib.Nifti1Image(img.get_fdata(),img.affine,img.header)    

    # save image
    img.to_filename(outfile)

    # clean up
    for fl in fls:
        hdr = fl.split('.')[0]+'.hdr'
        os.remove(fl)
        os.remove(hdr)

    return outfile

# def write_bnk_json(pth,base_json):
#     with open(pth, 'w') as fp:
#         json.dump(base_json, fp,sort_keys=True, indent=4)

def reorient_and_qform2sform(fname):
    '''
    BIDS validator was having trouble with information missing in the qform
    and sform. In reality, the qform info was present, but not the sform. 
    This code will copy the information over and populate the qform and
    sform such that they are "correct".

    this code lifted and lightly adapted from https://github.com/nipreps/niworkflows/blob/1adb83fa05824cd8ef14ef7da07254fab3a248c5/niworkflows/interfaces/images.py#L476'''

    error = False
    message = None

    orig_img = nib.load(fname)
    reoriented = nib.as_closest_canonical(orig_img)

    # Set target shape information (target is self in this case)
    target_zooms = orig_img.header.get_zooms()[:3]
    target_shape = np.array(orig_img.shape)[:3]
    target_span = target_shape * target_zooms

    zooms = np.array(reoriented.header.get_zooms()[:3])
    shape = np.array(reoriented.shape[:3])

    # Reconstruct transform from orig to reoriented image
    ornt_xfm = nib.orientations.inv_ornt_aff(
        nib.io_orientation(orig_img.affine), orig_img.shape
    )
    # Identity unless proven otherwise
    target_affine = reoriented.affine.copy()
    conform_xfm = np.eye(4)

    xyz_unit = reoriented.header.get_xyzt_units()[0]
    if xyz_unit == "unknown":
        # Common assumption; if we're wrong, unlikely to be the only thing that breaks
        xyz_unit = "mm"

    # Set a 0.05mm threshold to performing rescaling
    atol_gross = {"meter": 5e-5, "mm": 0.05, "micron": 50}[xyz_unit]
    # if 0.01 > difference > 0.001mm, freesurfer won't be able to merge the images
    atol_fine = {"meter": 1e-6, "mm": 0.001, "micron": 1}[xyz_unit]

    # Update zooms => Modify affine
    # Rescale => Resample to resized voxels
    # Resize => Resample to new image dimensions
    update_zooms = not np.allclose(zooms, target_zooms, atol=atol_fine, rtol=0)
    rescale = not np.allclose(zooms, target_zooms, atol=atol_gross, rtol=0)
    resize = not np.all(shape == target_shape)
    resample = rescale or resize
    if resample or update_zooms:
        # Use an affine with the corrected zooms, whether or not we resample
        if update_zooms:
            scale_factor = target_zooms / zooms
            target_affine[:3, :3] = reoriented.affine[:3, :3] @ np.diag(
                scale_factor
            )

        if resize:
            # The shift is applied after scaling.
            # Use a proportional shift to maintain relative position in dataset
            size_factor = target_span / (zooms * shape)
            # Use integer shifts to avoid unnecessary interpolation
            offset = (
                reoriented.affine[:3, 3] * size_factor - reoriented.affine[:3, 3]
            )
            target_affine[:3, 3] = reoriented.affine[:3, 3] + offset.astype(int)

        conform_xfm = np.linalg.inv(reoriented.affine) @ target_affine

        # Create new image
        data = reoriented.dataobj
        if resample:
            data = nli.resample_img(reoriented, target_affine, target_shape).dataobj
        reoriented = reoriented.__class__(data, target_affine, reoriented.header)

    #     out_name = fname

    transform = ornt_xfm.dot(conform_xfm)
    if not np.allclose(orig_img.affine.dot(transform), target_affine):
        message = "Original and target affines are not similar"
        print(message)
        error = True

    return reoriented, error, message


def set_xyzt_units(img, xyz='mm', t='sec'):
    '''Was getting an error about TR units not being correct in the image 
    header. Copied this code directly from Luke Chang's solution at
    https://neurostars.org/t/bids-validator-giving-error-for-tr/2538'''
    header = img.header.copy()
    header.set_xyzt_units(xyz=xyz, t=t)
    return img.__class__(img.get_fdata().copy(), img.affine, header) 

base_json = {'Modality': 'MR',
             'MagneticFieldStrength': 1.5,
             'Manufacturer': 'GE',
             'ManufacturersModelName': 'SIGNA_EXCITE',
             'InstitutionName': 'Bannockburn_Radiology_Center',
             'DeviceSerialNumber': '00000847GURNEEMR',
             'StationName': 'BRCMR',
             'PatientPosition': 'HFS',
             'SoftwareVersions': '11',
             'TaskName':'rest',
             'Instructions':'Eyes closed',
             'EchoTime': 0.0033,
             'RepetitionTime': 2,
             'FlipAngle': 85,
             'PulseSequenceType': '2D Spiral GRE Resting State fMRI',
             'SliceThickness':5,
            'SliceEncodingDirection':'k',
            'SequenceName':'sprlio',
            'SliceTiming':[0.0,0.08,0.16,0.24,0.32,0.4,0.48,0.56,0.64,0.72,0.8,0.88,0.96,1.04,1.12,1.2,1.28,1.36,1.44,1.52,1.6,1.68,1.76,1.84,1.92,2.0]
            }

def write_bnk_json(pth,base_json):
    with open(pth, 'w') as fp:
        json.dump(base_json, fp,sort_keys=True, indent=4)

def main():
    fmri_files = os.path.join(os.getcwd(), 'sub*', 'ses-*', 'func', '*.zip')

    zip_files = sorted(glob(fmri_files))

    for ii, fmri in enumerate(zip_files):
        file_path = Path(fmri)
        basename = os.path.basename(file_path)

        directory_path = file_path.parent
        
        sub_id = file_path.parents[2].name

        bids_name = basename.str.replace('.zip', '.nii.gz', regex=False)
        dl.unlock(dataset=sub_id, path=file_path)
        nii = convert_zipped_fmri(fmri, directory_path, bids_name)

        if nii == 'error':
            print('FAILURE!',fmri)
            continue
        # fix nii
        reoriented, error, message = reorient_and_qform2sform(nii)
        if error:
            raise ValueError(message)
        reoriented = set_xyzt_units(reoriented, xyz='mm', t='sec')
        os.remove(nii)
        reoriented.to_filename(nii)
        
        # dl.save(dataset=sub_id, message=f"Converted rs-fmri scan from analyze format to nifti in {directory_path}")
        jpth = file_path.split('.')[0]+'.json'
        write_bnk_json(jpth,base_json)
        # get rid of old files
        os.remove(fmri)

if __name__ == "__main__":
    
    main()