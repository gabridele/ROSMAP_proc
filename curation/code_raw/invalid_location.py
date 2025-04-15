import pandas as pd
import os
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import datalad.api as dl
import subprocess

def rename_file(file_path):
    """
    Rename a file by correcting the session identifier in its filename.

    Args:
        file_path (str): The full path to the file to be renamed.
    """
    file_path = file_path.lstrip('/')  # Remove leading slash if present
    filename = os.path.basename(file_path)  # Extract the filename
    directory = os.path.dirname(file_path)  # Extract the directory path
    print(f"Renaming: {file_path}")
    
    # Extract the correct session ID from the directory path
    correct_ses_match = re.search(r'/ses-(\d+)/', directory)
    # Extract the incorrect session ID from the filename
    incorrect_ses_match = re.search(r'ses-(\d+)', filename)
    
    # Proceed if both correct and incorrect session IDs are found
    if correct_ses_match and incorrect_ses_match:
        correct_ses = correct_ses_match.group(1)  # Correct session ID
        incorrect_ses = incorrect_ses_match.group(1)  # Incorrect session ID
        
        # Replace the incorrect session ID in the filename with the correct one
        new_filename = filename.replace(f'ses-{incorrect_ses}', f'ses-{correct_ses}')
        new_file_path = os.path.join(directory, new_filename)  # Construct the new file path
        print(f"Renaming: {file_path} -> {new_file_path}")

        try:
            # Rename the file
            os.rename(file_path, new_file_path)
            print(f"Renamed: {file_path} -> {new_file_path}")
        except Exception as e:
            # Handle errors during renaming
            print(f"Error renaming {file_path}: {e}")

def main():
    # Load the TSV file containing validation results
    df = pd.read_csv("code/CuBIDS/v0c_validation.tsv", sep="\t")

    # Filter for rows where the second column indicates 'INVALID_LOCATION'
    invalid_df = df[df.iloc[:, 1] == "INVALID_LOCATION"].copy()
    
    # Unlock each file listed in the filtered DataFrame
    for file_path in invalid_df.iloc[:, 0]:
        file_path = file_path.lstrip('/')  # Remove leading slash if present
        try:
            # Use DataLad to unlock the file
            subprocess.run(["datalad", "unlock", file_path], check=True)
            print(f"Unlocked: {file_path}")
        except subprocess.CalledProcessError as e:
            # Handle errors during unlocking
            print(f"Error unlocking {file_path}: {e}")

    # Rename files in parallel using a thread pool
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(rename_file, invalid_df.iloc[:, 0])  # Map the rename function to each file

if __name__ == "__main__":
    main()