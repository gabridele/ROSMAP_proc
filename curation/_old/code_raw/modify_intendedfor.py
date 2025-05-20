import os
import pandas as pd
import json
from concurrent.futures import ThreadPoolExecutor
import datalad.api as dl


# Function to modify 'INTENDED_FOR' entries
def modify_intended_for(input_file):
    """
    Modify the 'INTENDED_FOR' entry in the JSON file, in order to point to the right bold acq.
    Args:
        input_file (str): Path to the NIfTI file.
    """

    input_file = input_file.lstrip('/')
    # Replace the extension of the file from .nii.gz to .json
    if input_file.endswith(".nii.gz"):
        input_file = input_file.replace(".nii.gz", ".json")

    # Load the JSON file
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
    # Read the TSV file
    tsv_file = "code/CuBIDS/v2_validation.tsv"
    data = pd.read_csv(tsv_file, sep='\t')

    intended_for = data[data.iloc[:, 1] == 'INTENDED_FOR'].copy()
    
    # Unlock each file listed in the DataFrame
    for file in intended_for.iloc[:, 0]:
        file = file.lstrip('/')
        if file.endswith(".nii.gz"):
            file = file.replace(".nii.gz", ".json")
            dl.unlock(file)
        
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(modify_intended_for, intended_for.iloc[:, 0])

if __name__ == "__main__":
    main()
