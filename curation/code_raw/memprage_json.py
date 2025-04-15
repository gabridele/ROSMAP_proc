import os
import json
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import datalad.api as dl

def edit_json(nii_path):
    """
    Create new JSON sidecar file for a newly created NIfTI file to update metadata. Done after splitting multi-echo MEMPRAGE into single-echo files.
    
    Args:
        nii_path (str): Path to the NIfTI file.
    """
    # get path portions of filename
    # Remove leading slash if present
    json_file = nii_path.lstrip('/')
    json_path = str(json_file).replace(".nii.gz", ".json")
    
    input_name = nii_path.split('/')[-1]
    base_name = input_name.split('_T1w')[0]
    output_dir = '/'.join(nii_path.split('/')[:-1]).lstrip('/')
    
    # iterate thru the four echoes
    for i in range(1, 5):
        echo_name = f"{output_dir}/{base_name}_echo-{i}_T1w.nii.gz"
        json_echo = echo_name.replace(".nii.gz", ".json")
        
        try:
            # Check if echo_name exists
            if not os.path.exists(echo_name):
                print(f"Error: File {echo_name} does not exist.")
                continue

            # Open the JSON file safely
            with open(json_path, 'r') as f:
                metadata = json.load(f)

            # Update metadata based on the file name
            if f"_echo-{i}_" in json_echo:
                if i == 1:
                    metadata['NumVolumes'] = 1
                elif i == 2:
                    metadata['EchoTime'] = 0.00355
                elif i == 3:
                    metadata['EchoTime'] = 0.00541
                elif i == 4:
                    metadata['EchoTime'] = 0.00727
                metadata['NumVolumes'] = 1
                print(f"Updated metadata for {json_echo}")
            else:
                print(f"Error: Unrecognized echo type in {json_echo}")
                return

            metadata['NumVolumes'] = 1

            with open(json_echo, 'w') as f:
                json.dump(metadata, f, indent=4)
            print(f"Created {json_echo} successfully.")

        except FileNotFoundError:
            print(f"Error: JSON file not found for {json_file}")
        except json.JSONDecodeError:
            print(f"Error: Failed to decode JSON file {json_path}")
        except Exception as e:
            print(f"Unexpected error while processing {json_file}: {e}")

def main():
    df = pd.read_csv('code/CuBIDS/v0c_validation.tsv', sep='\t')

    memprage_df = df[df.iloc[:, 1] == 'T1W_FILE_WITH_TOO_MANY_DIMENSIONS'].copy()
    
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(edit_json, memprage_df.iloc[:, 0])
    
    # datalad remove multi-echo files and sidecar
    for file in memprage_df.iloc[:, 0]:
        file = file.lstrip('/')
        dl.remove(file)
        json_file = file.replace(".nii.gz", ".json")
        dl.remove(json_file)

if __name__ == '__main__':
    main()