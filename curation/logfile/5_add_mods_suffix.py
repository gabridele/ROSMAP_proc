import os
import pandas as pd
import numpy as np

## script to add modalities of non-batch_1 scans (because batch_1 had them already)
## suffixes taken from https://github.com/bids-standard/bids-specification/blob/master/src/schema/objects/suffixes.yaml
df_main = pd.read_csv("code/dates_scraping/4_dfscandates.csv")
df_main = df_main.copy()
df_main.insert(9, 'suffix', np.nan)
df_main.insert(10, 'acq', np.nan)
#df_main['acq'] = df_main['relative_path'].str.extract(r'(acq-[^_\/]+)(?=(_|VARIANT))', expand=False)
mask1 = df_main['batch_number'] == 'batch_1'
df_main.loc[mask1, 'acq'] = df_main.loc[mask1, 'relative_path'].str.extract(r'(acq-[^_/]+)(?=_|VARIANT|$)', expand=False)

# Mask for everything EXCEPT batch_1
mask = df_main['batch_number'] != 'batch_1'

# Ensure relative_path is string
df_main['relative_path'] = df_main['relative_path'].astype(str)

# Extract the last part of the relative path (filename)
last_part = df_main.loc[mask, 'relative_path'].str.split('/').str[-1]

# Remove extension (get only the filename without .nii.gz or similar)
basename = last_part.str.split('.').str[0]
basename = basename.str.replace(r'_ADC', '', regex=False)
basename = basename.str.replace(r'FLAIR', '', regex=False)
basename = basename.str.replace(r'T2', '', regex=False)
basename = basename.str.replace(r'T1', '', regex=False)
basename = basename.str.replace(r'qalas', '', regex=False)
basename = basename.str.replace(r'DWI', '', regex=False)
basename = basename.str.replace(r'-AP', '-', regex=False)
basename = basename.str.replace(r'-PA', '-', regex=False)

# Remove '-', '(', ')' from the basename
basename_clean = basename.str.replace(r'[-()_]+', '', regex=True)

# Build the 'acq' field
df_main.loc[mask, 'acq'] = (
    'acq-' +
    df_main.loc[mask, 'site'].astype(str) +
    df_main.loc[mask, 'protocol'].astype(str).str.replace('.0', '', regex=False) +
    basename_clean
)

#df_main['acq'] = acq_extracted[0]
#print(acq_extracted.shape, df_main.shape)

def define_modality(row):
    if row['batch_number'] == 'batch_1':
        return row['modality']
    elif row['batch_number'] == 'batch_2':
        if '3D_MPRAGE' in row['relative_path']:
            return 'anat'
        elif 'phase_map' in row['relative_path']:
            return 'fmap'
        elif 'rfMRI.zip' in row['relative_path']:
            return 'func'
    elif row['batch_number'] == 'batch_3':
        if 'DWI' in row['relative_path']:
            return 'dwi'
        elif any(x in row['relative_path'] for x in ['FLAIR', 'MEMPRAGE', 'T2']):
            return 'anat'
    elif row['batch_number'] == 'batch_4':
        return 'dwi'
    elif row['batch_number'] == 'batch_5':
        if row['site'] == 'UC':
            if '.PAR' in row['relative_path'] or '.REC' in row['relative_path'] or 'PARREC.zip' in row['relative_path'] or 'ParRec.zip' in row['relative_path'] or 'parrec.zip' in row['relative_path'] or 'Parrec.zip' in row['relative_path']:
                return np.nan
            elif 'extra' in row['relative_path'] or '.zip' in row['relative_path']:
                return np.nan
            elif 'DWI' in row['relative_path']:
                return 'dwi'
            elif 'HARDI' in row['relative_path']:
                return 'dwi'
            elif 'EPI.' in row['relative_path']:
                return 'func'
            elif 'EPI_' in row['relative_path']:
                return 'fmap'
            elif 'FieldMap' in row['relative_path']:
                return 'fmap'
            elif 'FLAIR' in row['relative_path']:
                return 'anat'
            elif 'MPRAGE' in row['relative_path']:
                return 'anat'
            elif 'mprage-rms' in row['relative_path']:
                return 'anat'
            elif 'ME-GRE' in row['relative_path']:
                return 'anat'
            elif 'T2' in row['relative_path']:
                return 'anat'
            elif 'SWI' in row['relative_path']:
                return 'anat'
        elif row['site'] == 'MG':
            if 'extra' in row['relative_path'] or '.zip' in row['relative_path']:
                return np.nan
            elif 'nii' not in row['relative_path']:
                return np.nan
            elif 'DTI' in row['relative_path']:
                return 'dwi'
            elif 'bold' in row['relative_path']:
                return 'func'
            elif 'FLAIR' in row['relative_path']:
                return 'anat'
            elif 'field_mapping' in row['relative_path']:
                return 'fmap'
            elif 't1' in row['relative_path']:
                return 'anat'
            elif 't2' in row['relative_path']:
                return 'anat'
            elif 'T2' in row['relative_path']:
                return 'anat'
        elif row['site'] == 'BNK':
            if '3D_MPRAGE' in row['relative_path']:
                return 'anat'
            elif 'phase_map' in row['relative_path']:
                return 'fmap'
            elif 'DTI' in row['relative_path']:
                return 'dwi'
            elif 'T2' in row['relative_path']:
                return 'anat'
            
        elif row['site'] == 'RIRC':
            if 'DWI' in row['relative_path']:
                return 'dwi'
            elif 'rsfMRI' in row['relative_path']:
                return 'func'
            elif any(x in row['relative_path'] for x in ['FLAIR', 'ME-GRE', 'MEMPRAGE', 'T2', 'qalas', '(Align_with_QALAS)']):
                return 'anat'

def define_suffix(row):
    if '.PAR' in row['relative_path'] or '.REC' in row['relative_path'] or 'PARREC.zip' in row['relative_path'] or 'ParRec.zip' in row['relative_path'] or 'parrec.zip' in row['relative_path'] or 'Parrec.zip' in row['relative_path'] or '_ph.' in row['relative_path']:
                return np.nan
    if row['batch_number'] == 'batch_1':
        if 'NoSWI' in row['relative_path']:
            return np.nan
        return row['relative_path'].split('_')[-1].split('.')[0]

    elif row['batch_number'] == 'batch_2':
        if '3D_MPRAGE' in row['relative_path']:
            return 'T1w'
        elif '_e1.' in row['relative_path']:
            return 'e1'
        elif '(phase_map)' in row['relative_path']:
            return 'phasemap'
        elif 'rfMRI.zip' in row['relative_path']:
            return 'bold'
    elif row['batch_number'] == 'batch_3':
        if '_ADC' in row['relative_path']:
            return 'ADC'
        elif 'Ax-DWI-PA-40-6.' in row['relative_path']:
            return 'dwi'
        elif 'FLAIR' in row['relative_path']:
            return 'FLAIR'
        elif 'MEMPRAGE' in row['relative_path']: # use T1w for now, but it should be split into multiple scans, given the multi-echo nature 
            return 'T1w'
        elif 'T2' in row['relative_path']:
            return 'T2w'
        
    elif row['batch_number'] == 'batch_4':
        if '_ADC' in row['relative_path']:
            return 'ADC'
        elif 'Ax-DWI-PA-40-' in row['relative_path']: 
            return 'dwi'
        elif 'HARDI' in row['relative_path']:
            return 'dwi'
        elif 'DTI' in row['relative_path']:
            return 'dwi'
        elif 'directions.' in row['relative_path']:
            return 'dwi'
        elif '_ColFA' in row['relative_path']:
            return 'ColFA'
        elif '_FA' in row['relative_path']:
            return 'FA'
        elif '_TRACEW' in row['relative_path']:
            return 'trace'
        
    elif row['batch_number'] == 'batch_5':
        if '_ADC' in row['relative_path']:
            return 'ADC'
        elif any(x in row['relative_path'] for x in ['Ax-DWI-AP-6-2.', 'Ax-DWI-PA-40-']):
            return 'dwi'
        elif 'rsfMRI' in row['relative_path']:
            return 'bold'
        elif 'bold' in row['relative_path']:
            return 'bold'
        elif 'FLAIR' in row['relative_path']:
            return 'FLAIR'
        elif 'ME-GRE' in row['relative_path']:
            return 'T1w'
        elif 'MEMPRAGE' in row['relative_path']: # use T1w for now, but it should be split into multiple scans, given the multi-echo nature 
            return 'T1w'
        elif 'MPRAGE' in row['relative_path']:
            return 'T1w'
        elif 'mprage' in row['relative_path']:
            return 'T1w'
        elif 'T2' in row['relative_path']:
            return 'T2w'
        elif 'qalas' in row['relative_path']: # use QALAS for now, but it should be split into multiple maps (T1map, T2map, PDmap)
            return 'QALAS'
        elif 'Align_with_QALAS' in row['relative_path']:
            return 'TB1map'
        elif 'HARDI' in row['relative_path']:
            return 'dwi'
        elif 'DTI' in row['relative_path']:
            return 'dwi'
        elif '(phase_map)' in row['relative_path']:
            return 'phasemap'
        elif 'field_mapping' in row['relative_path']:
            return 'fieldmap'
        elif 't1_mpr' in row['relative_path']:
            return 'T1w'
        elif 'EPI.' in row['relative_path']:
            return 'bold'
        elif 'EPI_' in row['relative_path']:
            return 'epi'
        elif 'FieldMap' in row['relative_path']:
            return 'fieldmap'
        elif '3D_MPRAGE' in row['relative_path']:
            return 'T1w'
        elif 'SWI' in row['relative_path']:
            return np.nan
        elif '_ph.' in row['relative_path']:
            return np.nan

    return np.nan


df_main['modality'] = df_main.apply(define_modality, axis=1)
df_main['suffix'] = df_main.apply(define_suffix, axis=1)

def check_quarantine_2025(row):
    if 'qalas' in row['relative_path']:
        return 'yes'
    elif 'Align_with_QALAS' in row['relative_path']:
        return 'yes'
    elif pd.isna(row['suffix']):
        return 'yes'
    
    return 'no'

df_main['quarantine_2025'] = df_main.apply(check_quarantine_2025, axis=1)

output_path = "code/dates_scraping/5_mods_suffix.csv"
df_main.to_csv(output_path, index=False)
