import os
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

# txt files containing the paths, aka ground truth, with scandate, protocol and visit info
txts = ["code/dates_scraping/bannockburn_directory_structure.txt", "code/dates_scraping/uc_directory_structure.txt", "code/dates_scraping/mg_directory_structure.txt", "code/dates_scraping/rirc_directory_structure.txt"]

df_main = pd.read_csv("code/dates_scraping/2_file_list_scandatesBIDSfile.csv")
df_main = df_main.copy()

#do some clean up
df_main['protocol'] = df_main['protocol'].astype(str).str.replace('.0', '', regex=False)
df_main['protocol'] = df_main['protocol'].str.replace('-', '', regex=False)
df_main['sub_id'] = df_main['sub_id'].astype(str).str.replace('.0', '', regex=False)
df_main['sub_id'] = df_main['sub_id'].str.lstrip('0')
df_main['visit_number'] = df_main['visit_number'].replace('0.0', '0')
df_main['visit_number'] = df_main['visit_number'].astype(str).str.replace('.0', '', regex=False)
df_main['raw_scandate'] = df_main['raw_scandate'].astype(str).replace('-', '', regex=False)
df_main['scandate_log210326'] = df_main['scandate_log210326'].astype(str).str.replace('.0', '', regex=False)
df_main['qc_scandate'] = df_main['qc_scandate'].astype(str).str.replace('.0', '', regex=False)

def process_paths_to_dataframe(paths):
    """
    Process a list of paths to extract required information and return a DataFrame.

    Args:
        paths (list): List of file paths.

    Returns:
        pd.DataFrame: DataFrame with columns ['sub', 'visit_number', 'site', 'protocol', 'ScanDate'].
    """
    processed_data = []
    for path in paths:
        parts = path.split('/')[-3:]  # [site, protocol, 'ScanDate_VisitNumber_SubID']
        if len(parts) < 3:
            continue
        site = parts[0].upper()
        if site == "BANNOCKBURN":
            site = "BNK"
        if parts[1].isdigit() and len(parts[1]) == 6:
            protocol = '20' + parts[1]  # This is usually '210326' -> '20210326'
        else:
            protocol = parts[1]  # Use as-is if not 6 digits
        
        try:
            scan_date, visit_number, sub = parts[2].split('_')
            sub = sub.lstrip('0')
            visit_number = visit_number.lstrip('0') or '0' 
        except ValueError:
            print(f"Skipping path with unexpected format: {path}")
            continue

        scandate = '20' + scan_date  # Convert to full year
        
        processed_data.append([sub, visit_number, site, protocol, scandate])

    return pd.DataFrame(processed_data, columns=['sub', 'visit_number', 'site', 'protocol', 'ScanDate'])


dfs = {}
for fil in txts:
    with open(fil, "r") as file:
        paths = [line.strip() for line in file.readlines()]
        df_temp = process_paths_to_dataframe(paths)
        df_name = fil.split('/')[-1].replace('_directory_structure.txt', '_df')
        dfs[df_name] = df_temp

def process_chunk_txt(chunk, df_txt):
    row = chunk.iloc[0]  # all rows in chunk refer to same session
    batch = str(row['batch_number'])
    subdir = str(row['sub_id'])
    visit_num = str(row['visit_number'])
    protocol = str(row['protocol'])
    site = str(row['site'])

    match = df_txt[
        (df_txt['sub'] == subdir) &
        (df_txt['visit_number'] == visit_num) &
        (df_txt['protocol'] == protocol) &
        (df_txt['site'] == site)
    ]

    if not match.empty:
        # Update df_txt once per folder/session
        if not pd.isna(df_txt.loc[match.index, 'match']).all():
            df_txt.loc[match.index, 'match'] = df_txt.loc[match.index, 'match'].astype(str) + f', x{batch}'
        else:
            df_txt.loc[match.index, 'match'] = f'x{batch}'
        protocol_rush = match.iloc[0]['protocol']
        scandate = match.iloc[0]['ScanDate']
        for i in chunk.index:
            df_main.at[i, 'rush_scandate'] = scandate
            df_main.at[i, 'rush_folder'] = scandate + '_' + visit_num + '_' + subdir
            df_main.at[i, 'rush_visit'] = visit_num
            df_main.at[i, 'rush_site'] = site
            df_main.at[i, 'rush_protocol'] = protocol_rush

    return df_main, df_txt

# Filter out rows where the batch_number column is 'batch_1' because I dont have visit number for these. gonna do a separate script for these
df_main_no1 = df_main[df_main['batch_number'] != 'batch_1']

df_bannockburn = dfs['bannockburn_df']
df_uc = dfs['uc_df']
df_mg = dfs['mg_df']
df_rirc = dfs['rirc_df']

# Save the processed DataFrames to CSV files
for name, df in dfs.items():
    output_file = f"code/dates_scraping/{name}.csv"
    df.to_csv(output_file, index=False)
    print(f"Saved {name} to {output_file}")
    
df_bnk = pd.read_csv("code/dates_scraping/bannockburn_df.csv", sep=',', dtype=str, header=0)
df_uc = pd.read_csv("code/dates_scraping/uc_df.csv", sep=',', dtype=str, header=0)
df_mg = pd.read_csv("code/dates_scraping/mg_df.csv", sep=',', dtype=str, header=0)
df_rirc = pd.read_csv("code/dates_scraping/rirc_df.csv", sep=',', dtype=str, header=0)

df_bnk['match'] = pd.NA
df_uc['match'] = pd.NA
df_mg['match'] = pd.NA
df_rirc['match'] = pd.NA

grouped = df_main_no1.groupby('folder_name')

for folder_name, group in grouped:
    site = group.iloc[0]['site']
    
    if site == 'BNK':
        df_main, df_bnk = process_chunk_txt(group, df_bnk)
    elif site == 'UC':
        df_main, df_uc = process_chunk_txt(group, df_uc)
    elif site == 'MG':
        df_main, df_mg = process_chunk_txt(group, df_mg)
    elif site == 'RIRC':
        df_main, df_rirc = process_chunk_txt(group, df_rirc)

# for chunk in chunks:
#     df_main, df_bnk = process_chunk_txt(chunk, df_bnk)
df_bnk.to_csv("code/dates_scraping/bannockburn_df_updated.csv", index=False)
df_uc.to_csv("code/dates_scraping/uc_df_updated.csv", index=False)
df_mg.to_csv("code/dates_scraping/mg_df_updated.csv", index=False)
df_rirc.to_csv("code/dates_scraping/rirc_df_updated.csv", index=False)

df_main = df_main.sort_values(by=df_main.columns[0], ascending=True)
df_main.to_csv("code/dates_scraping/3_file_list_rushdates.csv", index=False)

## all folders from rush that match, have an x in the match column in the file

## what about batch_1 having session with multiple protocols inside? or different sessions ....

## also, should merge modalities from diff batches but same sub and session together. because some duplicates might just be different modalities downloaded at differten times




