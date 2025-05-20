import os
import json
import pandas as pd

def extract_json_info(bids_root, output_csv):
    data = []

    # Iterate over batch folders (e.g., batch_1, batch_2, batch_3)
    for batch_folder in os.listdir(bids_root):
        
        if not batch_folder.startswith("batch_"):
            continue  # Skip non-batch folders
        
        batch_path = os.path.join(bids_root, batch_folder)
        if not os.path.isdir(batch_path):
            continue  # Ensure it's a directory

        batch_name = os.path.basename(batch_path)  # Extract batch folder name
        print(batch_path)
        # Iterate over subject folders
        for sub_folder in os.listdir(batch_path):
            
            if not sub_folder.startswith("sub-"):
                continue  # Skip non-subject folders
            print(sub_folder)
            sub_id = sub_folder.replace("sub-", "")  # Extract subject ID
            sub_path = os.path.join(batch_path, sub_folder)
            
            for ses_folder in os.listdir(sub_path):
                if ses_folder.startswith("."):
                    continue
                ses_id = "N/A"  # Default session ID if it's not found

                if batch_name == "batch_1" and ses_folder.startswith("ses-"):
                    # For batch 1, ses_id is extracted directly from the 'ses-*' format
                    ses_id = ses_folder.replace("ses-", "")
                else:
                    # For batch 2 and 3, ses_id is the second part between underscores in format ######_##_#########
                    parts = ses_folder.split("_")
                    if len(parts) > 1:  # Ensure there are at least two underscores
                        ses_id = parts[1]  # Extract the session number from between underscores

                ses_path = os.path.join(sub_path, ses_folder)
                
                if not os.path.isdir(ses_path):
                    continue  # Ensure it's a directory
                print("ses_id:", ses_id)
                print("ses_path:", ses_path)
                # Iterate over JSON files inside the session folder
                for root, _, files in os.walk(ses_path):
                    for file in files:
                        if file.endswith(".json"):
                            json_path = os.path.join(root, file)
                            
                            try:
                                with open(json_path, 'r', encoding='utf-8') as f:
                                    json_data = json.load(f)

                                # Extract ScanDate and AcquisitionTime if available
                                scan_date = json_data.get("ScanDate", "N/A")
                                acq_time = json_data.get("AcquisitionTime", "N/A")

                                # Append data
                                data.append([batch_name, sub_id, ses_id, file, scan_date, acq_time])

                            except Exception as e:
                                print(f"Error reading {json_path}: {e}")

    # Save data to CSV
    df = pd.DataFrame(data, columns=["batch", "sub_id", "ses_id", "json_file", "ScanDate", "AcquisitionTime"])
    df.to_csv(output_csv, index=False)
    print(f"CSV saved: {output_csv}")


bids_root = "/home/gabridele/backup/ROSMAP_proc/sourcedata/"  
output_csv = "/home/gabridele/backup/extracted_scan_info1.csv" 
extract_json_info(bids_root, output_csv)
