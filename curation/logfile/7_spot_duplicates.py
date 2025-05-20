import os
import pandas as pd
import numpy as np
## 167 scans are without rush_scandate, not counting those with issues (otherwise it'd be 807). ao basically those are excluded from the duplicates check
## how to deal with those?
## they're missing the qc date because most of them are flagged as quarantined. but those who are not, idk why they god logBIDS scandate but not qc date....??
## should spot another way to check for duplicates
df_main = pd.read_csv("code/dates_scraping/6_highlightissues_batch1.csv", header=0, dtype=str)
df_main = df_main.copy()

# Drop rows if the entry of the first row contains specific values
values_to_check = ['backup_', '.PAR', '.REC', 'PARREC.zip', 'ParRec.zip', 'parrec.zip', 'Parrec.zip', '.pdf', '.txt', 'nifti.zip', 'niftii.zip', 'dicomdata', 'dicomdata.zip', 'extra', 'PCASL', 'BAT.', 'Localizer', 'Perfusion_Weighted', 'relCBF', 'WIP', '.7']
print('Listing subjects where theres inconsistency in protocol and site, within each session')
df_main = df_main[~df_main.iloc[:, 0].apply(lambda x: any(value in str(x) for value in values_to_check))]

# Reset the index after dropping
df_main = df_main.reset_index(drop=True)
filtered_df = df_main[df_main['suffix'].notna()]
filtered_df = filtered_df[filtered_df['issue'].isna() | (filtered_df['issue'] == '')].copy()

#### since batch_1 doesnt have visit number, for now i'll just look at duplicates in the rest of the batches
grouped = filtered_df.groupby(['sub_id', 'rush_scandate'])

def check_consistency(group):
    if not group.empty:
        for col in ['protocol', 'site']:
            if group[col].nunique() > 1:
                print(f"Inconsistent '{col}' in sub_id={group['sub_id'].iloc[0]}, rush_visit={group['rush_visit'].iloc[0]}")
                if col == 'site':
                    print(group['site'])
    return group


def flag_duplicates(group):
    group = group.copy()
    filenames = group['relative_path'].str.split('/').str[-1]
    duplicated_mask = filenames.duplicated(keep='first')
    
    # Mark duplicates
    group['duplicates'] = duplicated_mask.map({True: 'duplicate', False: ''})
    
    # Mark the "mother" of the duplicates with 'x' in the duplicates column
    group.loc[~duplicated_mask & filenames.duplicated(keep=False), 'duplicates'] = 'x'
    
    return group

def annotate_duplicate_batches(group):
    filenames = group['relative_path'].str.split('/').str[-1]
    first_batches = {}

    for i, fname in enumerate(filenames):
        if fname not in first_batches:
            first_batches[fname] = group.iloc[i]['batch_number']

    group = group.copy()

    group['duplicates'] = [
        f'duplicate to {first_batches[fname]}' if dup == 'duplicate' else dup
        for fname, dup in zip(filenames, group['duplicates'])
    ]

    return group
df_checked = grouped.apply(check_consistency)
df_checked = df_checked.reset_index(drop=True)

df_intermediate = df_checked.groupby(['sub_id', 'rush_visit']).apply(flag_duplicates).reset_index(drop=True)
df_final = df_intermediate.groupby(['sub_id', 'rush_visit']).apply(annotate_duplicate_batches).reset_index(drop=True)

# Merge to include batch_1

df_main = df_main.merge(
    df_final[['sub_id', 'rush_visit', 'relative_path', 'duplicates']],
    on=['sub_id', 'rush_visit', 'relative_path'],
    how='left'
)

print('----------------------')
print("Number of duplicates found:", df_intermediate['duplicates'].value_counts().get('duplicate', 0))

df_main.to_csv("code/dates_scraping/7_spot_duplicates.csv", index=False)