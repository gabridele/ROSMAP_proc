import os
import subprocess
import pandas as pd

# Load the dataframe
df_main = pd.read_csv("code/dates_scraping/9_logfile.csv")


for index, row in df_main.iterrows():
    src_path = row['relative_path']
    dest_path = row['new_path']

    if not os.path.isfile(src_path):
        print(f"Source file not found: {src_path}")
        continue

    # Create destination directory if it doesn't exist
    dest_dir = os.path.dirname(dest_path)
    os.makedirs(dest_dir, exist_ok=True)

    rsync_cmd = [
        'rsync', '-av', '--copy-links',  # Use --copy-links to dereference symlinks
        src_path,
        dest_path
    ]

    try:
        subprocess.run(rsync_cmd, check=True)
        print(f"[OK] {src_path} → {dest_path}")
    except subprocess.CalledProcessError as e:
        print(f"[ERR] rsync failed for {src_path} → {dest_path}: {e}")
