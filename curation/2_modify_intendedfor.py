import os
import pandas as pd
import json
from concurrent.futures import ThreadPoolExecutor
import datalad.api as dl
import glob
from pathlib import Path

# Function to modify 'INTENDED_FOR' entries
def modify_intended_for(input_file):
    """
    Modify the 'INTENDED_FOR' entry in the JSON file, in order to point to the right bold acq.
    Args:
        input_file (str): Path to the NIfTI file.
    """

    # Load the JSON file
    # Skip processing if 'BNK' is in the input file name
    
    try:
        with open(input_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"File not found: {input_file}")
        return
    except json.JSONDecodeError:
        print(f"Error decoding JSON in file: {input_file}")
        return
    
    # Modify the 'INTENDED_FOR' entry with correct func scan
    if "IntendedFor" in data:
        if 'BNK' in input_file:
            data.pop("IntendedFor", None)
            print(f"Skipping file {input_file} as BNK scans do not need fmap correction")
            return
        try:
            input_dir = os.path.relpath(os.path.dirname(input_file), start=os.getcwd())
            func_dir = os.path.abspath(os.path.join(input_dir, "../func/"))
            data["IntendedFor"] = [
                os.path.relpath(os.path.join(func_dir, file), start=os.path.join(os.getcwd(), input_dir.split("ses-")[0]))

                for file in os.listdir(func_dir) if file.endswith("_bold.nii.gz")
            ]
            print(f"Modified 'IntendedFor' entry in {input_file}")
            print(data["IntendedFor"])
        except FileNotFoundError:
            print("Directory '../func' not found.")
            return

    #Save the modified JSON file
    try:
        with open(input_file, 'w') as f:
            json.dump(data, f, indent=4)
    except Exception as e:
        print(f"Error writing to file {input_file}: {e}")

def main():
    fmap_dirs = os.path.join(os.getcwd(), 'sub*', 'ses-*', 'fmap', '*.json')

    fmaps = sorted(glob(fmap_dirs))

    for ii, fmap in enumerate(fmaps):
        file_path = Path(fmap)
        dl.unlock(file_path)

        modify_intended_for(file_path)

if __name__ == "__main__":
    main()
