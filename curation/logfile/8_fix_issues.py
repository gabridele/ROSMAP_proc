import os
import pandas as pd
import numpy as np

# Load main dataframe
df_main = pd.read_csv("code/dates_scraping/7_spot_duplicates.csv", header=0, dtype=str)
df_main = df_main.copy()

# ----------------------------------------------------------
# Fix known issues
# ----------------------------------------------------------

# Drop problematic scan due to known site mixup. UC has complete data in other batch, visit 8.
df_main.loc[
    df_main.iloc[:, 0] == 'sub-59497970_ses-0_task-rest_acq-UC20120221EPI_bold.',
    'suffix'
] = np.nan

# Manual correction for a specific entry in relative_path
pattern = 'batch_1/sub-73146926/ses-0/'
mask = df_main.iloc[:, 0].astype(str).str.contains(pattern)

df_main.loc[mask, 'rush_site'] = 'UC'
df_main.loc[mask, 'rush_protocol'] = '120221'
df_main.loc[mask, 'rush_scandate'] = '111207'
df_main.loc[mask, 'rush_visit'] = '1'

# ----------------------------------------------------------
# Load external site-specific dataframes
# ----------------------------------------------------------

df_bnk = pd.read_csv("code/dates_scraping/bannockburn_df_updated.csv", dtype=str)
df_uc = pd.read_csv("code/dates_scraping/uc_df_updated.csv", dtype=str)
df_mg = pd.read_csv("code/dates_scraping/mg_df_updated.csv", dtype=str)
df_rirc = pd.read_csv("code/dates_scraping/rirc_df_updated.csv", dtype=str)

# ----------------------------------------------------------
# Export version with suffix only, dropping scans that are superfluous
# ----------------------------------------------------------

df_new = df_main[df_main['suffix'].notna()]


mask = (
    (df_new['sub_id'] == '59497970') &
    (df_new['session'] == 1) &
    (df_new['batch_number'] == 'batch_1') &
    (df_new['modality'].isin(['fmap', 'func']))
)

df_new.loc[mask, 'session'] = 2
df_new.loc[mask, 'rush_scandate'] = '20140814'
df_new.loc[mask, 'rush_visit'] = '10'
df_new.loc[mask, 'site'] = 'UC'


df_new.to_csv("code/dates_scraping/8_updated_df.csv", index=False)
# ----------------------------------------------------------
# Filter out rows for quarantine and protocol issues
# ----------------------------------------------------------

df_new['sub_id'] = df_new['sub_id'].astype(str).str.replace('.0', '', regex=False)

df_new['scandate_log210326'] = df_new['scandate_log210326'].astype(str).str.replace('.0', '', regex=False)
df_filter = df_new[
    ~(
        (df_new['quarantine_2025'] == 'yes') |
        (df_new['issue'].str.contains('site, protocol', na=False)) |
        (df_new['issue'].str.contains('protocol', na=False))
    )
]

# ----------------------------------------------------------
# Fill missing rush data from external site dataframes
# ----------------------------------------------------------

df_nan = df_filter[df_filter['rush_scandate'].isna()]

for _, row in df_nan.iterrows():
    batch = str(row['batch_number'])
    subdir = str(row['sub_id'])
    scandate_log = str(row['scandate_log210326'])
    
    for df in [df_bnk, df_uc, df_mg, df_rirc]:
        #print('Processing:', df)
        match = df[
            (df['sub'] == subdir) & 
            (df['ScanDate'] == scandate_log)
            ]
        #print('Match:', match)
        if not match.empty:
            scandate = match.iloc[0]['ScanDate']
            visit_num = match.iloc[0]['visit_number']
            protocol_rush = match.iloc[0]['protocol']
            site = row['site']

            # Update only rows matching this subject and scan date
            update_mask = (
                (df_new['sub_id'] == subdir) &
                (df_new['scandate_log210326'] == scandate_log)
            )

            df_new.loc[update_mask, 'rush_scandate'] = scandate
            df_new.loc[update_mask, 'rush_folder'] = f"{scandate}_{visit_num}_{subdir}"
            df_new.loc[update_mask, 'rush_visit'] = visit_num
            df_new.loc[update_mask, 'rush_site'] = site
            df_new.loc[update_mask, 'rush_protocol'] = protocol_rush

# ----------------------------------------------------------
# Convert rush_visit to integer if possible
# ----------------------------------------------------------
df_new['rush_visit'] = df_new['rush_visit'].apply(
    lambda x: int(float(x)) if pd.notna(x) else x
)

# ----------------------------------------------------------
# Define new session number based on sub_id and rush_visit
# ----------------------------------------------------------
df_new = df_new[~df_new['duplicates'].str.contains('duplicate to', na=False)]

def assign_session_numbers(group):
    group = group.copy()
    group['rush_visit'] = pd.to_numeric(group['rush_visit'], errors='coerce')
    
    # Only rank non-NaN rush_visits
    if group['rush_visit'].isna().all():
        group['new_ses'] = np.nan
    else:
        group = group.sort_values('rush_visit')
        group['new_ses'] = group['rush_visit'].rank(method='dense').astype('Int64') - 1

    return group

df_new = (
    df_new
    .groupby(['sub_id'], group_keys=False)
    .apply(assign_session_numbers)
)

## ----------------------------------------------------------
# get rid of batch_1 bnk duplicates
## ----------------------------------------------------------

def remove_bnk_duplicates(group):
    is_anat = group['modality'] == 'anat'
    is_fmap = group['modality'] == 'fmap'
    
    anat_count = is_anat.sum()
    fmap_count = is_fmap.sum()

    # If condition is met: 8 anat and 4 fmap entries
    if anat_count == 8 and fmap_count == 4:
        # Set suffix to NaN only for batch 2 entries
        mask = group['batch_number'] == 'batch_2'
        group.loc[mask, 'suffix'] = np.nan
        group.loc[mask, 'duplicates'] = 'duplicate to batch_1'

    return group

df_neww = (
    df_new
    .groupby(['sub_id', 'new_ses'], group_keys=False)
    .apply(lambda g: remove_bnk_duplicates(g) if g['rush_site'].iloc[0] == 'BNK' else g)
)
df_new.update(df_neww)

## ----------------------------------------------------------
## MG has got no duplicates, i checked. because they either from batch_1 or batch_4 (the latter being only dwi)
## ----------------------------------------------------------

## ----------------------------------------------------------
## UC, some duplicates found in T1w and bold with batch_5
## ----------------------------------------------------------

def remove_uc_duplicates(group):
    is_anat = group['suffix'] == 'T1w'
    is_func = group['suffix'] == 'bold'
    is_fmap = group['modality'] == 'fmap'

    anat_count = is_anat.sum()
    func_count = is_func.sum()
    fmap_count = is_fmap.sum()
    # If condition is met: 8 anat and 4 func entries
    if anat_count == 4 and func_count == 4:
        # Set suffix to NaN only for batch 2 entries
        mask = (group['batch_number'] == 'batch_5') & (
            group['suffix'].isin(['T1w', 'bold'])
        )
        group.loc[mask, 'suffix'] = np.nan
        group.loc[mask, 'duplicates'] = 'duplicate to batch_1'
    if fmap_count == 8:
        # Set suffix to NaN only for batch 2 entries
        mask = (group['batch_number'] == 'batch_5') & (
            group['suffix'].isin(['fmap'])
        )
        group.loc[mask, 'suffix'] = np.nan
        group.loc[mask, 'duplicates'] = 'duplicate to batch_1'
    return group

df_neww = (
    df_new
    .groupby(['sub_id', 'new_ses'], group_keys=False)
    .apply(lambda g: remove_uc_duplicates(g) if g['rush_site'].iloc[0] == 'UC' else g)
)
df_new.update(df_neww)

# keeping this duplicate because i will drop the same bold scans from batch_1 as they have issues, and keep these
folder_name_to_exclude = "121101_08_59497970"
df_new.loc[df_new['folder_name'].str.contains(folder_name_to_exclude, na=False), 'duplicates'] = np.nan
## and set those from batch1 to quarantine
df_new.loc[
    (df_new['sub_id'] == '59497970') & 
    (df_new['rush_visit'] == 8) & 
    (df_new['batch_number'] == 'batch_1'),
    'quarantine_2025'
] = 'yes'

df_new = df_new[~df_new['duplicates'].str.contains('duplicate to', na=False)]
df_neww = (
    df_new
    .groupby(['sub_id'], group_keys=False)
    .apply(assign_session_numbers)
)
df_new.update(df_neww)
## ----------------------------------------------------------
## no more duplicates apparently
## ----------------------------------------------------------
mask = (
    (df_new['rush_protocol'] == '120221') &
    (df_new['sub_id'] == '73146926')
)

df_new.loc[mask, 'rush_protocol'] = '20120221'
df_new.loc[mask, 'rush_scandate'] = '20111207'

mask = (
    (df_new['rush_scandate'] == '20110727') &
    (df_new['sub_id'] == '37865636')
)

df_new.loc[mask, 'rush_site'] = 'BNK'

df_new.loc[
    df_new['issue'].str.contains('site, protocol', na=False),
    'quarantine_2025'
] = 'yes'


## get cleaned up version of df, without quarantined
df_new = df_new[~df_new['quarantine_2025'].str.contains('yes', na=False)]

df_new.to_csv("code/dates_scraping/8_cleanedup_df.csv", index=False)

## looks promising, but theres still duplicates from batch_1 and the rest of the batches. gotta find a way to spot them
## for now we skipped site, protocol issue 
## discard SWI as they are not necessary for fmriprep, nay theyre not used at all for the pipelines

## site, protocol issue still persists but they dont appear in latest csv saving. 
## protocol only issue gets solved easily

## putting 'site, protocol' issues in quarantine because BNK I got them already in batch2, so theyre duplicates, and the MGs I requested Chris


