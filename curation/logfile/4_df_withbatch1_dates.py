# can exclude rirc and bnk because batch1 doesnt have them. so only mg, uc and bnk txt filez

### BUT - you should add as many x as the matches it found. tbd tmr
### ALSO to add column in log where you indicate those that you dropped cause of duplicates or not useful (e.g. localizer)

import os
import pandas as pd
import numpy as np

df_main = pd.read_csv("code/dates_scraping/3_file_list_rushdates.csv")
df_main = df_main.copy()

df_main_1 = df_main[df_main['batch_number'] == 'batch_1']

df_main_1['protocol'] = df_main_1['protocol'].astype(str).str.replace('.0', '', regex=False)
df_main_1['protocol'] = df_main_1['protocol'].str.replace('-', '', regex=False)
df_main_1['sub_id'] = df_main_1['sub_id'].astype(str).str.replace('.0', '', regex=False)
df_main_1['sub_id'] = df_main_1['sub_id'].str.lstrip('0')
df_main_1['session'] = df_main_1['session'].replace('0.0', '0')
df_main_1['session'] = df_main_1['session'].astype(str).str.replace('.0', '', regex=False)
df_main_1['raw_scandate'] = df_main_1['raw_scandate'].astype(str).replace('-', '', regex=False)
df_main_1['scandate_log210326'] = df_main_1['scandate_log210326'].astype(str).str.replace('.0', '', regex=False)
df_main_1['qc_scandate'] = df_main_1['qc_scandate'].astype(str).str.replace('.0', '', regex=False)

# no rirc because batch_1 doesnt have them
df_bnk = pd.read_csv("code/dates_scraping/bannockburn_df.csv", sep=',', dtype=str, header=0)
df_uc = pd.read_csv("code/dates_scraping/uc_df.csv", sep=',', dtype=str, header=0)
df_mg = pd.read_csv("code/dates_scraping/mg_df.csv", sep=',', dtype=str, header=0)

def process_chunk_txt(chunk, df_txt):
    row = chunk.iloc[0]  # all rows in chunk refer to same session
    batch = str(row['batch_number'])
    subdir = str(row['sub_id'])
    protocol = str(row['protocol'])
    qc_date = str(row['qc_scandate'])
    site = str(row['site'])

    match = df_txt[
        (df_txt['sub'] == subdir) &
        (df_txt['ScanDate'] == qc_date)
    ]

    if not match.empty:
        # Update df_txt once per folder/session
        if not pd.isna(df_txt.loc[match.index, 'match_1']).all():
            df_txt.loc[match.index, 'match_1'] = df_txt.loc[match.index, 'match_1'].astype(str) + ', x'
        else:
            df_txt.loc[match.index, 'match_1'] = 'x'

        scandate = match.iloc[0]['ScanDate']
        visit_num = match.iloc[0]['visit_number']
        protocol_rush = match.iloc[0]['protocol']
        for i in chunk.index:
            df_main.at[i, 'rush_scandate'] = scandate
            df_main.at[i, 'rush_folder'] = scandate + '_' + visit_num + '_' + subdir
            df_main.at[i, 'rush_visit'] = visit_num
            df_main.at[i, 'rush_site'] = site
            df_main.at[i, 'rush_protocol'] = protocol_rush

    return df_main, df_txt
    
df_bnk = pd.read_csv("code/dates_scraping/bannockburn_df_updated.csv", sep=',', dtype=str, header=0)
df_uc = pd.read_csv("code/dates_scraping/uc_df_updated.csv", sep=',', dtype=str, header=0)
df_mg = pd.read_csv("code/dates_scraping/mg_df_updated.csv", sep=',', dtype=str, header=0)
df_rirc = pd.read_csv("code/dates_scraping/rirc_df_updated.csv", sep=',', dtype=str, header=0)

df_bnk['match_1'] = pd.NA
df_uc['match_1'] = pd.NA
df_mg['match_1'] = pd.NA
df_rirc['match_1'] = pd.NA

grouped = df_main_1.groupby(['session', 'sub_id'])

for _, group in grouped:
    site = group.iloc[0]['site']
    
    if site == 'BNK':
        df_main, df_bnk = process_chunk_txt(group, df_bnk)
    elif site == 'UC':
        df_main, df_uc = process_chunk_txt(group, df_uc)
    elif site == 'MG':
        df_main, df_mg = process_chunk_txt(group, df_mg)
        
df_bnk.to_csv("code/dates_scraping/bannockburn_df_updated.csv", index=False)
df_uc.to_csv("code/dates_scraping/uc_df_updated.csv", index=False)
df_mg.to_csv("code/dates_scraping/mg_df_updated.csv", index=False)

df_main = df_main.sort_values(by=df_main.columns[0], ascending=True)
df_main.to_csv("code/dates_scraping/4_dfscandates.csv", index=False)