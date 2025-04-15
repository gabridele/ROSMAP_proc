import os
import pandas as pd

def inspect_and_rename_sessions(csv_path, log_path="/home/gabridele/backup/duplicates_log.txt"):
    # Step 1: Load the CSV containing session details
    df = pd.read_csv(csv_path)

    # Step 2: Group by subject ID and ScanDate, then sort by ScanDate and AcquisitionTime
    df['ScanDate'] = pd.to_datetime(df['ScanDate'])  # Convert ScanDate to datetime format
    df['AcquisitionTime'] = pd.to_datetime(df['AcquisitionTime'], errors='coerce')  # Handle any invalid times
    
    # Sort by ScanDate and AcquisitionTime to get chronological order
    df_sorted = df.sort_values(by=["sub_id", "ScanDate", "AcquisitionTime"])
    
    # Step 3: Assign session numbers based on chronological order of sessions per subject
    session_number = 1  # Start session numbering from 1
    session_dict = {}  # Keep track of sessions by sub_id and ScanDate
    df_sorted['new_ses_id'] = None  # Create a new column for the updated session IDs

    with open(log_path, 'w') as log_file:
        for idx, row in df_sorted.iterrows():
            session_key = (row['sub_id'], row['ScanDate'])  # Group by subject and scan date

            # Check if we've already assigned a session number for this subject and scan date
            if session_key in session_dict:
                # If it's a duplicate session, log it and keep the same session number
                original_session_num = session_dict[session_key]
                log_msg = f"Duplicate session found for subject {row['sub_id']} on {row['ScanDate']} at {row['AcquisitionTime']}. " \
                          f"Original session: ses-{str(original_session_num).zfill(2)} (Session #{original_session_num}), " \
                          f"Duplicate session: ses-{str(session_number).zfill(2)} (Session #{session_number})\n"
                log_file.write(log_msg)
                df_sorted.at[idx, 'new_ses_id'] = f"ses-{str(original_session_num).zfill(2)}"  # Assign the original session number
            else:
                # Assign a new session number for the first occurrence of this session
                df_sorted.at[idx, 'new_ses_id'] = f"ses-{str(session_number).zfill(2)}"
                session_dict[session_key] = session_number  # Track this session as assigned
                session_number += 1  # Increment session number for the next unique session

    print(f"Duplicates (if any) logged to {log_path}")

    inspection_csv_path = csv_path.replace(".csv", "_inspection.csv")
    df_sorted.to_csv(inspection_csv_path, index=False)

    print(f"Inspection CSV saved as: {inspection_csv_path}")
    print(f"\nHere is a preview of the sorted sessions and their new IDs:\n")
    print(df_sorted[['batch', 'sub_id', 'ses_id', 'new_ses_id', 'ScanDate', 'AcquisitionTime']])

output_csv = "/home/gabridele/backup/extracted_scan_info1.csv" 
inspect_and_rename_sessions(output_csv)