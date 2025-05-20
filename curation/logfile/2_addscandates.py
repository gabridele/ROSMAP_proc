import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# adds scandates from bids file and qc. only for batch_1 supposedly because these are the files Jake had when he created the csv's

# Load data
df = pd.read_csv(f"code/dates_scraping/1_file_list.tsv", sep='\t', dtype=str)

df['scandate_log210326'] = np.nan
df['quarantine_log210326'] = np.nan
df['qc_scandate'] = np.nan
df['rush_folder'] = np.nan
df['rush_scandate'] = np.nan
df['rush_visit'] = np.nan
df['rush_site'] = np.nan
df['rush_protocol'] = np.nan

df_20210326 = pd.read_csv(f"code/dates_scraping/BIDS_data_log_20210326.csv", sep=',', dtype=str, header=0)
df_qc = pd.read_csv(f"code/dates_scraping/qc.csv", sep=',', dtype=str, header=0)

# Filter df to include only rows with 'batch_1'
df_batch_1 = df[df['batch_number'] == 'batch_1'].copy()

# Function to process a chunk of rows
def process_chunk1(chunk):
    for i, row in chunk.iterrows():
        subdir = row['sub_id']
        sesdir = row['session'].replace('.0', '')

        # Check if subdir and sesdir match any row in df_20210326
        match = df_20210326[
            (df_20210326['subdir'].str.replace('sub-', '', regex=False) == subdir) &
            (df_20210326['Session'].str.replace('.0', '') == sesdir)
        ]

        if not match.empty:
            scandate = match.iloc[0]['ScanDate'].replace('.0', '')
            quarantine = match.iloc[0]['quarantine'] 
            chunk.at[i, 'scandate_log210326'] = scandate
            chunk.at[i, 'quarantine_log210326'] = quarantine
    return chunk

def process_chunk2(chunk):
    for i, row in chunk.iterrows():
        subdir = row['sub_id'].lstrip('0')
        sesdir = row['session'].replace('.0', '')

        # Check if subdir and sesdir match any row in df_20210326
        match = df_qc[
            (df_qc['sub'].str.replace('sub-', '', regex=False) == subdir) &
            (df_qc['ses'].str.replace('.0', '') == sesdir)
        ]

        if not match.empty:
            scandate = match.iloc[0]['ScanDate'].replace('.0', '')
            chunk.at[i, 'qc_scandate'] = scandate
    return chunk

# Split df_batch_1 into chunks
num_chunks = 48
chunks = np.array_split(df_batch_1, num_chunks)

with ProcessPoolExecutor() as executor:
    results = executor.map(process_chunk1, chunks)
df_batch_1 = pd.concat(results)
updated_idx = df_batch_1.index[df_batch_1['scandate_log210326'].notna()]
df.loc[updated_idx, 'scandate_log210326'] = df_batch_1.loc[updated_idx, 'scandate_log210326']
df.loc[updated_idx, 'quarantine_log210326'] = df_batch_1.loc[updated_idx, 'quarantine_log210326']

with ProcessPoolExecutor() as executor:
    results = executor.map(process_chunk2, chunks)
df_batch_1 = pd.concat(results)
updated_idx = df_batch_1.index[df_batch_1['qc_scandate'].notna()]
df.loc[updated_idx, 'qc_scandate'] = df_batch_1.loc[updated_idx, 'qc_scandate']

# Apply specific rules for sub_id 59497970
for i, row in df.iterrows():
    if row['sub_id'] == '59497970':
        if row['session'] == '2' and row['site'] == 'UC' and row['modality'] != 'anat':
            df.at[i, 'qc_scandate'] = '20121101'
            df.at[i, 'session'] = '1'
            df.at[i, 'site'] = 'UC'
        elif row['session'] == '1' and row['site'] == 'MG' and row['modality'] != 'anat':
            df.at[i, 'qc_scandate'] = '20140814'
            df.at[i, 'session'] = '2'
            df.at[i, 'site'] = 'MG'
            
# Save the updated DataFrame
df.to_csv("code/dates_scraping/2_file_list_scandatesBIDSfile.csv", index=False, sep=',')

### it works!!! now apply to other sources to get dates
