import os
import shutil
import datalad.api as dl
import glob
from pathlib import Path

def datalad_clone_and_copy():
    """
    Clone and copy datasets using DataLad.
    This function processes subject directories in the 'raw' directory,
    checks for corresponding datasets in the 'sourcedata/batch_1' directory,
    and copies them into the appropriate locations.
    """
    # Define the base directories
    raw = os.path.abspath('raw')  # Absolute path to the 'raw' directory
    source_base = os.path.abspath('sourcedata/batch_1')  # Absolute path to the source directory
    tmpdir = os.path.join(raw, 'tmpdir')  # Temporary directory for processing

    # Iterate through all sub-* directories in the 'raw' directory
    for sub_dir in sorted(os.listdir(raw)):
        sub_path = os.path.join(raw, sub_dir)  # Full path to the subject directory
        
        # Check if the current item is a directory and starts with 'sub-'
        if os.path.isdir(sub_path) and sub_dir.startswith('sub-'):
            source_dataset = os.path.join(source_base, sub_dir)  # Path to the source dataset
            cloned_path = os.path.join(tmpdir, sub_dir)  # Path to the temporary cloned dataset
            
            # Check if the source dataset exists
            if os.path.exists(source_dataset):
                print(f'copying {source_dataset} into {sub_dir}')  # Log the copy operation

                # Collect all non-hidden files/directories to copy
                dirs_to_copy = []
                for f in os.listdir(source_dataset):
                    if f.startswith('.'):
                        continue  # Skip hidden files
                    full_path = os.path.join(source_dataset, f)  # Full path to the file/directory
                    dirs_to_copy.append(full_path)

                # Perform the copy operation if there are files/directories to copy
                if dirs_to_copy:
                    dl.copy_file(dirs_to_copy, f"{sub_path}/{f}", recursive=True, dataset=sub_path)

            else:
                # Log if the source dataset is not found
                print(f'Skipping {sub_dir}, source dataset not found at {source_dataset}')

if __name__ == '__main__':
    datalad_clone_and_copy()
