import pandas as pd
import numpy as np
import re

def format_date(date_str):
    if pd.isna(date_str) or date_str == 'NA':
        return 'NA'
    s = date_str.strip()
    if len(s) == 6 and s.isdigit():            # YYMMDD
        return f"20{s[:2]}{s[2:4]}{s[4:6]}"
    if len(s) == 8 and s.isdigit():            # YYYYMMDD
        return f"{s[:4]}{s[4:6]}{s[6:8]}"
    if len(s) == 10 and re.match(r'\d{4}\d{2}\d{2}', s):
        return s                               # already YYYY‑MM‑DD
    return 'NA'

def pad_sub_id(sub_id):
    if pd.isna(sub_id) or sub_id == 'NA':
        return 'NA'
    digits = ''.join(filter(str.isdigit, str(sub_id)))
    return digits.zfill(8) if digits else 'NA'

def clean_visit_number(visit):
    s = str(visit).strip().upper()
    if s in ('NA', ''):
        return None
    s = s.lstrip('0') or '0'       # interpret '00' as '0'
    s = s.replace('.0', '')        
    return int(s) if s.isdigit() else None

def extract_visit_from_folder(folder_name):
    """
    Given a folder name like '50300372_15', return the visit number (15).
    Returns None if no trailing '_<number>'.
    """
    m = re.search(r'_(\d+)$', folder_name)

    return int(m.group(1))

SITE_CODES = {'BNK', 'MG', 'UC', 'RIRC'}

def extract_site_and_protocol(filename):
    if 'acq-' not in filename:
        return 'NA', 'NA'
    m = re.search(r'acq-((?:[A-Za-z]+)?(\d{6,8}))', filename)
    if m:
        site, prot = m.group(1), m.group(2)
        site = site if site in SITE_CODES else 'NA'
        return site, prot
    return 'NA', 'NA'

BATCH2 = {'site': 'BNK', 'protocol': '20090211'}
BATCH3 = {'site': 'RIRC'}

def update_batch3_protocol(sub_id):
    return '20230726' if sub_id in {'13044513', '70876731'} else '20230803'

PROTOCOL_TO_SITE = {
    '20090211': 'BNK','20140922': 'UC','20150706': 'UC','20151120': 'UC',
    '20160125': 'UC','20120221': 'UC','20220113': 'UC','20220419': 'UC',
    '20120501': 'MG','20150715': 'MG','20160621': 'MG','20160627': 'MG'
}

def format_and_fill(input_tsv, output_tsv):
    df = pd.read_csv(input_tsv, sep='\t', dtype=str)
    # Drop rows where the first column contains 'batch_1/dataset_description.json'
    df = df[~df.iloc[:, 0].str.contains('batch_1/dataset_description.json', na=False)]

    # core transforms
    df['raw_scandate'] = df['scan_date'].apply(format_date)
    df['sub_id']       = df['sub_id'].apply(pad_sub_id)
    df['visit_number'] = df['visit_number'].apply(clean_visit_number)
    df = df.drop(columns=['scan_date'])

    # insert new cols
    idx = df.columns.get_loc('batch_number') + 1
    for c in ('site','protocol','folder_name'):
        df.insert(idx, c, 'NA')
        idx += 1

    # batch‑wise logic
    for i, row in df.iterrows():
        path, batch = row[df.columns[0]], row['batch_number']
        parts = path.split('/')

        if batch == 'batch_1':
            site, prot = extract_site_and_protocol(path)
            df.at[i, 'site']     = site
            df.at[i, 'protocol'] = prot

        elif batch == 'batch_2':
            df.at[i, 'site']        = BATCH2['site']
            df.at[i, 'protocol']    = BATCH2['protocol']
            df.at[i, 'folder_name'] = parts[1] if len(parts) > 1 else 'NA'

        elif batch == 'batch_3':
            df.at[i, 'site']        = BATCH3['site']
            df.at[i, 'protocol']    = update_batch3_protocol(row['sub_id'])
            df.at[i, 'folder_name'] = parts[2] if len(parts) > 2 else 'NA'

        elif batch == 'batch_4':
            if len(parts) >= 6:
                site_candidate = parts[3]
                df.at[i, 'site']     = site_candidate if site_candidate in SITE_CODES else 'NA'
                df.at[i, 'protocol'] = format_date(parts[4])
                folder               = parts[5]
                df.at[i, 'folder_name'] = folder
                # new: extract visit number from folder if missing
                if pd.isna(df.at[i, 'visit_number']):
                    df.at[i, 'visit_number'] = extract_visit_from_folder(folder)

        elif batch == 'batch_5':
            site_candidate = parts[1].upper()
            df.at[i, 'site']     = site_candidate if site_candidate in SITE_CODES else 'NA'

            if site_candidate == 'MG':
                df.at[i, 'protocol'] = '20120501'
                df.at[i, 'folder_name'] = parts[2]
            elif site_candidate == 'BNK': 
                df.at[i, 'protocol'] = '20090211'
                df.at[i, 'folder_name'] = parts[2]
            else:
                df.at[i, 'protocol'] = format_date(parts[2])
                df.at[i, 'folder_name'] = parts[3]

    # back‑fill site from protocol where still NA
    mask = df['site'] == 'NA'
    df.loc[mask, 'site'] = df.loc[mask, 'protocol'].map(PROTOCOL_TO_SITE).fillna('NA')

    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved to: {output_tsv}")

if __name__ == "__main__":
    format_and_fill(
        input_tsv="code/dates_scraping/0_file_list.tsv",
        output_tsv="code/dates_scraping/1_file_list.tsv"
    )
