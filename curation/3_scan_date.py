import os
import re
import json
import pandas as pd
import datalad.api as dl

#### code to extract ScanDate param from qc.csv file and add it into json file

def update_scan_date(bids_root, csv_path):
    df = pd.read_csv(csv_path, dtype={"sub_id": str, "new_ses": str, "rush_scandate": str})
    
    # Filter DataFrame for rows where 'new_path' ends with '.json' and strip '../raw'
    df["new_path"] = df["new_path"].str.strip()
    json_funcs = df[df["new_path"].str.endswith(".json")]["new_path"].str.replace("../raw", "", regex=False)
    
    ## because bnk jsons are not in the df, as they were created after bidsification
    bnk_funcs = df[df["new_path"].str.endswith(".zip")]["new_path"].str.replace("../raw", "", regex=False)
    bnk_funcs = bnk_funcs.str.replace(".zip", ".json", regex=False)

    # Combine json_funcs and bnk_funcs
    json_paths = pd.concat([json_funcs, bnk_funcs], ignore_index=True)
    # Iterate over filtered paths
    for json_rel_path in json_paths:
        json_abs_path = os.path.join(bids_root, json_rel_path.lstrip("/"))  # Construct absolute path
        if os.path.exists(json_abs_path):  # Ensure the path exists
            try:
                with open(json_abs_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                # Extract ScanDate from DataFrame
                match = df[df["new_path"].str.endswith(json_rel_path)]
                if not match.empty:
                    scan_date = str(match.iloc[0]["rush_scandate"])
                    if "." in scan_date:
                        scan_date = scan_date.split(".")[0]
                    
                    formatted_date = f"{scan_date[:4]}-{scan_date[4:6]}-{scan_date[6:]}"
                    
                    # Update or add ScanDate
                    data["ScanDate"] = formatted_date

                    # Write back to file
                    with open(json_abs_path, 'w', encoding='utf-8') as f:
                        json.dump(data, f, indent=4)
                    print(f"Updated {json_abs_path} with ScanDate: {formatted_date}")
            except Exception as e:
                print(f"Error processing {json_abs_path}: {e}")

# Example usage
bids_root = "/home/gabridele/backup/ROSMAP_proc/raw"  
csv_path = "/home/gabridele/backup/ROSMAP_proc/sourcedata/code/dates_scraping/9_logfile.csv"
update_scan_date(bids_root, csv_path)