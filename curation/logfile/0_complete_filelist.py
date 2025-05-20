import os
import csv
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

print("Starting file processing...")

# Precompiled patterns
folder_pattern = re.compile(r'(\d{6})_(\d{2})_(\d+)')
sub_id_pattern = re.compile(r'sub-(\d+)')
ses_id_pattern = re.compile(r'ses-(\d+)')

def is_hidden_or_ignored(path_parts, ignore_dirs):
    return any(part.startswith('.') or part in ignore_dirs for part in path_parts)

### ADD ses- for batch_1 because u gonna need it for matching with qc and bids log

def extract_info_from_path(path):
    # Assuming the path follows a structure like:
    # - batch_5/rirc/230711/230713_07_88463949/88463949_07_nii/Ax-DWI-AP-6-2.bval
    # - batch_2/110810_02_95638231/3D_MPRAGEa.json
    # - batch_3/sub-59500816/241118_09_59500816/Sag-3D-T2-SPACE.nii.gz
    # - batch_4/radc-dti/dti-2025-03/MG/120501/31373325_05/DIFF_DTI_45_directions.bvec
    # - batch_1/sub-00228190/ses-0/anat/sub-00228190_ses-0_acq-20150706MPRAGE_T1w.nii.gz
    
    parts = path.split(os.sep)
    
    batch_number = parts[0]
    site = parts[1]
    session_number = 'NA'
    scan_date = 'NA'
    visit_number = 'NA'
    sub_id = 'NA'
    modality = 'NA'

    # Handle batch_5
    # Example: batch_5/rirc/230711/230713_07_88463949/88463949_07_nii/Ax-DWI-AP-6-2.bval
    if batch_number == 'batch_5':
        if site == 'mg' or site == 'bnk':
            match = folder_pattern.match(parts[2])
        else:
            match = folder_pattern.match(parts[3])  # 110810_02_95638231
        if match:
            scan_date, visit_number, sub_id = match.groups()
            
    # Handle batch_4: skip scan_date
    elif batch_number == 'batch_4':
        sub_id = parts[-2].split('_')[0]
        modality = 'dwi'        

    # Handle batch_3: sub-59500816/241118_09_59500816, extract scan_date, visit, sub_id
    elif batch_number == 'batch_3':
        match = folder_pattern.match(parts[2])
        if match:
            scan_date, visit_number, sub_id = match.groups()

    # Handle batch_2: scan_date_visitor_sub_id
    elif batch_number == 'batch_2':
        match = folder_pattern.match(parts[1])
        if match:
            scan_date, visit_number, sub_id = match.groups()

    # Handle batch_1: sub-<sub_id>/ses-<session_id>/anat/<filename>
    elif batch_number == 'batch_1':
        match = sub_id_pattern.match(parts[1])
        if match:
            sub_id = match.group(1)
            ses_match = ses_id_pattern.match(parts[2])
            if ses_match:
                session_number = ses_match.group(1)
            if 'anat' in parts[3]:
                modality = 'anat'
            elif 'fmap' in parts[3]:
                modality = 'fmap'
            elif 'func' in parts[3]:
                modality = 'func'

    visit_number = visit_number.replace('.0','')
    return batch_number, session_number, scan_date, visit_number, sub_id, modality

def get_all_file_paths(base_dir, ignore_dirs=None):
    if ignore_dirs is None:
        ignore_dirs = set()

    all_paths = []
    for root, dirs, files in os.walk(base_dir):
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ignore_dirs]
        for file in files:
            if not file.startswith('.'):
                rel_path = os.path.relpath(os.path.join(root, file), base_dir)
                if not is_hidden_or_ignored(rel_path.split(os.sep), ignore_dirs):
                    all_paths.append(rel_path)
    return all_paths

def extract_all_info_parallel(paths, max_workers=48):
    file_list = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(extract_info_from_path, path): path for path in paths}
        for future in as_completed(futures):
            path = futures[future]
            try:
                result = future.result()
                file_list.append([path] + list(result))
            except Exception as e:
                print(f"Error with {path}: {e}")
    return file_list

def save_to_tsv(file_list, output_file):
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['relative_path', 'batch_number', 'session', 'scan_date', 'visit_number', 'sub_id', 'modality'])
        writer.writerows(file_list)

if __name__ == "__main__":
    base_directory = "."
    output_tsv = "code/dates_scraping/0_file_list.tsv"
    ignore = {'batch_2old', 'code'}
    print('getting the file paths')
    paths = get_all_file_paths(base_directory, ignore_dirs=ignore)
    print(f"Found {len(paths)} files.")

    files = extract_all_info_parallel(paths)
    print(f"Processed {len(files)} files.")

    # check_for_duplicates(files)
    save_to_tsv(files, output_tsv)

    print(f"Saved to: {output_tsv}")
