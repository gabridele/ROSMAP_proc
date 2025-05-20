import pandas as pd
import numpy as np
import re
df_main = pd.read_csv("code/dates_scraping/8_cleanedup_df.csv", header=0, dtype=str)
df_main = df_main.copy()

df_main['protocol'] = df_main['protocol'].replace('0.0', '0')
df_main['protocol'] = df_main['protocol'].astype(str).str.replace('.0', '', regex=False)
df_main['visit_number'] = df_main['visit_number'].replace('0.0', '0')
df_main['visit_number'] = df_main['visit_number'].astype(str).str.replace('.0', '', regex=False)
df_main['rush_scandate'] = df_main['rush_scandate'].replace('0.0', '0')
df_main['rush_scandate'] = df_main['rush_scandate'].astype(str).str.replace('.0', '', regex=False)
df_main['rush_visit'] = df_main['rush_visit'].replace('0.0', '0')
df_main['rush_visit'] = df_main['rush_visit'].astype(str).str.replace('.0', '', regex=False)
df_main['rush_protocol'] = df_main['rush_protocol'].replace('0.0', '0')
df_main['rush_protocol'] = df_main['rush_protocol'].astype(str).str.replace('.0', '', regex=False)

## initialize new_acq string and new_path
df_main['new_acq'] = np.nan
df_main['new_path'] = np.nan
# new_acq
# new path
# creating the folders

## ----
## creating acq- string
## ----

mask = df_main['batch_number'] != 'batch_1'
last_part = df_main.loc[mask, 'relative_path'].str.split('/').str[-1]

# Remove extension (get only the filename without .nii.gz or similar)
basename = last_part.str.split('.').str[0]
basename = basename.str.replace(r'_ADC', '', regex=False)
basename = basename.str.replace(r'FLAIR', '', regex=False)
basename = basename.str.replace(r'T2', '', regex=False)
basename = basename.str.replace(r'T1', '', regex=False)
basename = basename.str.replace(r'qalas', '', regex=False)
basename = basename.str.replace(r'DWI', '', regex=False)
basename = basename.str.replace(r'-AP', '', regex=False)
basename = basename.str.replace(r'-PA', '', regex=False)
basename = basename.str.replace(r'_e1', '', regex=False)
basename = basename.str.replace(r'_echoes_e\d+', 'echoes', regex=True)
basename = basename.str.replace(r'VARIANT[^_]*_', '', regex=True)

# Remove '-', '(', ')' from the basename
basename_clean = basename.str.replace(r'[-()_]+', '', regex=True)

df_main.loc[mask, 'new_acq'] = (
    'acq-' +
    df_main.loc[mask, 'rush_site'].astype(str) +
    df_main.loc[mask, 'rush_protocol'].astype(str) +
    basename_clean
)

mask1 = df_main['batch_number'] == 'batch_1'
df_main.loc[mask1, 'new_acq'] = (
    'acq-' +
    df_main.loc[mask1, 'rush_site'].astype(str) +
    df_main.loc[mask1, 'rush_protocol'].astype(str) +
    df_main.loc[mask1, 'relative_path']
        .str.replace(r'VARIANT[^_]*_', '_', regex=True)
        .str.extract(r'acq-[A-Za-z]*\d{8}([A-Za-z0-9]+)_', expand=False)
        .fillna('')
)

## ------
## creating new path
## -----

df_main['file_ext'] = df_main['relative_path'].str.extract(r'([a-zA-Z0-9]+(?:\.gz)?)$')

## still gotta handle multi echo
## task-rest
## dir-
## other suffixes

def build_new_path(row):
    rel_path = str(row['relative_path'])
    
    # Handle echo
    echo_match = re.search(r'_echoes_e(\d+)', rel_path)
    echo2_match = re.search(r'_mapping_e(\d+)', rel_path)
    if echo_match:
        echo_part = f'_echo-{echo_match.group(1)}'
    elif echo2_match: 
        echo_part = f'_echo-{echo2_match.group(1)}'
    else:
        echo_part = ''

    # Direction tag
    if '-PA' in rel_path:
        direction = '_dir-PA'
    elif '-AP' in rel_path:
        direction = '_dir-AP'
    else:
        direction = ''
    
    # task-rest logic
    task_rest = '_task-rest' if row['suffix'] == 'bold' else ''

    # File extension
    file_ext = row['file_ext']
    
    
    # Build the full new_path
    return (
        '../raw/' +
        'sub-' + str(row['sub_id']).zfill(8) +
        '/ses-' + str(row['new_ses']) +
        '/' + str(row['modality']) +
        '/sub-' + str(row['sub_id']).zfill(8) +
        '_ses-' + str(row['new_ses']) +
        task_rest +
        '_' + str(row['new_acq']) +
        echo_part +
        direction +
        '_' + str(row['suffix']) +
        '.' + file_ext
    )

df_main['new_path'] = df_main.apply(build_new_path, axis=1)
print(df_main[['new_path']])  # Check if 'new_path' is correct
df_main.to_csv("code/dates_scraping/9_logfile.csv", index=False)
