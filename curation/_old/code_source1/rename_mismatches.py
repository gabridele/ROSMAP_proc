import os
import pandas as pd
import re
import datalad.api as dl

# Load the CSV file containing mismatched entries
df = pd.read_csv("code/filtered_false_entries.csv") 

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    # Extract relevant fields from the row
    sub, ses, acq, protocol, group = row["sub"], row["ses"], row["acq"], row["Protocol"], row['ScannerGroup']
    sub = str(sub).zfill(8)  # Zero-pad the subject ID to 8 digits

    # Extract the date part from the acquisition field using regex
    match = re.search(r"^[A-Z]{2}(\d{8})", acq)
    if not match:
        # Skip if no valid date is found in the acquisition field
        print(f"Skipping {acq} - No valid date found")
        continue
    
    # Debugging: Print the current row and extracted values
    print(row)
    acq_date = match.group(1)  # Extracted date from the acquisition field
    acq_prefix = acq[:2]  # Extract the prefix (first two characters) of the acquisition field
    print(acq_prefix)
    print(acq_date)
    protocol_date = str(protocol)  # Ensure the protocol date is a string
    print(protocol)

    # Update the acquisition prefix if it does not match the scanner group
    if acq_prefix != group:
        acq = group + acq[2:]  # Replace the prefix with the scanner group
        print(f"Updated acq: {acq}")

    # Skip renaming if the acquisition date matches the protocol date
    if acq_date == protocol_date:
        continue
    else:
        # Construct the old and new file prefixes
        old_prefix = f"sub-{sub}_ses-{ses}_task-rest_acq-{acq}_bold"
        new_prefix = old_prefix.replace(acq_date, protocol_date)
        
        # Define the folder path and root path for the files
        folder_path = f"sub-{sub}/ses-{ses}/func/"
        root_path = "/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_1"
        
        # Iterate through file extensions to rename both `.nii.gz` and `.json` files
        for ext in [".nii.gz", ".json"]:
            old_file = os.path.join(folder_path, old_prefix + ext)  # Old file path
            new_file = os.path.join(folder_path, new_prefix + ext)  # New file path
            print(old_file, new_file)  # Debugging: Print old and new file paths
            complete_filepath = os.path.join(root_path, old_file)  # Full path to the old file
            print(complete_filepath)  # Debugging: Print the complete file path

            try:
                # Check if the file exists before attempting to rename
                if not os.path.exists(complete_filepath):
                    print(f"Skipping {complete_filepath} - File does not exist")
                    continue
                
                # Unlock the file using Datalad
                dl.unlock(dataset=f"sub-{sub}", path=complete_filepath)
                
                # Rename the file
                os.rename(old_file, new_file)
                print(f"Renamed: {old_file} -> {new_file}")
                    
            except OSError as e:
                # Handle errors during the renaming process
                print(f"Error renaming {old_file} to {new_file}: {e}")
            except FileNotFoundError:
                # Handle cases where the file is not found
                print(f"File not found: {old_file}")
        
        # Save changes
        dl.save(dataset=f"sub-{sub}", message=f"Renamed files whose acqs had mismatch in protocol and scanner group")