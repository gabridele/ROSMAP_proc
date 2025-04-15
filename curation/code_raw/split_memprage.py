import nibabel as nib
import numpy as np
import pandas as pd
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

## script to split multi-echo mprage files into single-echo 3D volumes

def split_memprage(input_file):
    """
    Split a 4D NIfTI file into multiple 3D volumes.
    Args:
        input_file (str): Path to the NIfTI file.
    """
    input_file = input_file.lstrip('/')
    # Load the 4D NIfTI image
    nii_img = nib.load(input_file)
    data = nii_img.get_fdata()

    # Get the number of volumes in the fourth dimension
    num_volumes = data.shape[3]
    # Split into 3D volumes
    volumes = [data[:, :, :, i] for i in range(num_volumes)]

    # Save each 3D volume separately
    for i, vol in enumerate(volumes, start=1):
        vol_img = nib.Nifti1Image(vol, nii_img.affine)
        input_name = input_file.split('/')[-1]
        base_name = input_name.split('_T1w')[0]
        output_dir = '/'.join(input_file.split('/')[:-1])
        nib.save(vol_img, f"{output_dir}/{base_name}_echo-{i}_T1w.nii.gz")
        print(f"Saved {output_dir}/{base_name}_echo-{i}_T1w.nii.gz")
    print(f"Split {input_file} into {num_volumes} 3D volumes successfully!")

def main():
    df = pd.read_csv('code/CuBIDS/v0c_validation.tsv', sep='\t')

    memprage_df = df[df.iloc[:, 1] == 'T1W_FILE_WITH_TOO_MANY_DIMENSIONS'].copy()
    for file in memprage_df.iloc[:, 0]:
        file = file.lstrip('/')
        try:
            subprocess.run(["datalad", "unlock", file], check=True)
            print(f"Unlocked {file}")
        except subprocess.CalledProcessError as e:
            print(f"Error unlocking {file}: {e}")
    
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(split_memprage, memprage_df.iloc[:, 0])
    

if __name__ == '__main__':
    main()