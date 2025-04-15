import os
import re
import json
import datalad.api as dl

def update_scan_date(bids_root):
    # Iterate over subject folders
    for sub_folder in os.listdir(bids_root):
        sub_path = os.path.join(bids_root, sub_folder)
        print('sub_path:', sub_path)
        if not os.path.isdir(sub_path) or not sub_folder.startswith("sub-"):
            continue  # Skip non-subject directories
        
        # Iterate over session folders
        for ses_folder in os.listdir(sub_path):
            ses_path = os.path.join(sub_path, ses_folder)
            
            if not os.path.isdir(ses_path):
                continue  # Skip if not a directory
            
            # Split folder name and extract session ID and date
            parts = ses_folder.split("_")
            if len(parts) < 3:
                continue  # Skip if naming doesn't match expected format
            
            date = parts[0][:6]  # First 6 digits are the date
            
            print('parts:', parts)
            print('date:', date)
            
            formatted_date = f"20{date[:2]}-{date[2:4]}-{date[4:]}" # Convert to YYYY-MM-DD[Z]
            print(formatted_date)
            
            # Process all JSON files in the session folder
            for root, _, files in os.walk(ses_path):
                for file in files:
                    if file.endswith(".json"):
                        json_path = os.path.join(root, file)
                        print(json_path)
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
bids_root = "/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2"  # Change this to your actual BIDS dataset root
update_scan_date(bids_root)