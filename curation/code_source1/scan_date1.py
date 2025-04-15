import os
import re
import json
import pandas as pd
import datalad.api as dl

#### code to extract ScanDate param from qc.csv file and add it into json file

def update_scan_date(bids_root, csv_path):
    df = pd.read_csv(csv_path, dtype={"sub": str, "ses": str, "ScanDate": str})
    
    # Iterate over subject folders
    for sub_folder in os.listdir(bids_root):
        sub_path = os.path.join(bids_root, sub_folder)

        if not os.path.isdir(sub_path) or not sub_folder.startswith("sub-"):
            continue  # Skip non-subject directories
            
        sub_id = sub_folder.replace("sub-", "")  # Extract subject ID
        print(sub_id)

        # Iterate over session folders
        for ses_folder in os.listdir(sub_path):
            ses_path = os.path.join(sub_path, ses_folder)
            
            if not os.path.isdir(ses_path) or not ses_folder.startswith("ses-"):
                continue  # Skip if not a directory
            
            parts = ses_folder.split("_")
            ses_id = [x.split("-")[1] for x in parts]
            print(ses_id)
            # Convert DataFrame columns to string
            df["sub"] = df["sub"].astype(float).astype(int).astype(str).str.zfill(8)
            df["ses"] = df["ses"].astype(str).str.strip()

            # Convert session ID properly (remove leading zeros)
            sub_id = str(sub_id).strip()
            if isinstance(ses_id, list):
                ses_id = ses_id[0]  # Extract first element if it's a list

            ses_id = str(int(ses_id))  # "00" → "0" to match DataFrame

            # Debugging: Print types
            print(f"Type of sub_id: {type(sub_id)}, Value: {sub_id}")
            print(f"Type of ses_id: {type(ses_id)}, Value: {ses_id}")
            match = df[(df["sub"].eq(sub_id)) & (df["ses"].eq(ses_id))]
            # Perform matching with explicit Series comparison
            if not match.empty:  # Ensure match exists
                scan_date = str(match.iloc[0]["ScanDate"])  # Convert to string
                if "." in scan_date:  # Only split if decimal exists
                    scan_date = scan_date.split(".")[0]
            formatted_date = f"{scan_date[:4]}-{scan_date[4:6]}-{scan_date[6:]}"
            print('date:', formatted_date)
            # if 1 == 21: 
            #     raise Exception('hi')
            # Process all JSON files in the session folder
            for root, _, files in os.walk(ses_path):
                for file in files:
                    if file.endswith(".json"):
                        json_path = os.path.join(root, file)
                        dl.unlock(dataset=sub_path, path=json_path)
                        try:
                            with open(json_path, 'r', encoding='utf-8') as f:
                                data = json.load(f)
                            
                            # Update or add ScanDate
                            data["ScanDate"] = formatted_date

                            # Write back to file
                            with open(json_path, 'w', encoding='utf-8') as f:
                                json.dump(data, f, indent=4)
                            print(f"Updated {json_path} with ScanDate: {formatted_date}")
                        except Exception as e:
                            print(f"Error processing {json_path}: {e}")
        dl.save(dataset=sub_path, message=f"Added ScanDate in json sidecars")

# Example usage
bids_root = "/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_1"  
csv_path = "/home/gabridele/backup/qc.csv"
update_scan_date(bids_root, csv_path)