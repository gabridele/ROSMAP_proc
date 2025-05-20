import os
import re
import datalad.api as dl
from pathlib import Path
import concurrent.futures
import glob


def rename_bids_dwi(filename):
    """Renames BIDS files by replacing underscores with dashe in dwi acquisition blocks."""
    match = re.search(r'(acq-)([^_]+)(_dir-PA)(_.*?)(_(dwi)\..*)$', filename)
    if match:
        acq_block = match.group(2) + match.group(4).replace('_', '').replace('-', '') + match.group(3)
        #print(f"Renaming: {filename} -> {match.group(1)}{acq_block}{match.group(5)}")
        return filename.replace(match.group(0), f"{match.group(1)}{acq_block}{match.group(5)}")

    return None  # No match found

def process_file(path_):
    """Handles renaming logic while avoiding unwanted modifications.
    
    - Skips hidden directories (paths containing "/.").
    - Skips files that already follow the correct BIDS format.
    - Renames dwi files if necessary.
    """
    # Ignore hidden files/directories
    if os.path.basename(str(path_)).startswith("."):
        return None
    # Skip already correctly formatted anatomical/functional files
    if re.search(r'acq-[A-Z]{1,4}\d{6,8}_(T1w|T2w|FLAIR|bold|e1|magnitude\d?|phase\d?)\..*$', str(path_)):
        return None  # Skip already correct filenames
    
    # Rename DWI files if they match the pattern
    if re.search(r'_acq-[A-Z]{1,4}\d{6}_dir-PA_([^_]*)(_(dwi)\..*)?$', str(path_)):
        return path_, rename_bids_dwi(str(path_))

    return None  # No changes needed

def main():

    paths_to_process = [
        path for path in sorted(glob.glob("sub-*/ses-*/dwi/*"))
        if not path.startswith("sub-41635233")
    ]

    # Step 1: Parallelize Renaming Decisions
    renamed_files = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        results = executor.map(process_file, paths_to_process)
        for result in results:
   
            if result and result[1]:  # If renaming needed
                renamed_files.append(result)
                
    # Step 2: Unlock and rename Files (Sequential, to avoid conflicts)
    for old_path, new_filename in renamed_files:
        sub = old_path.split("/")[0]
        print(f"Unlocking files in dataset: {sub}")
        dl.unlock(dataset=sub, path=old_path)
        print(f"Renaming: {old_path} -> {new_filename}")
        os.rename(old_path, new_filename)

    
        dl.save(dataset=sub, message="Renamed dwi acquisition block in BIDS filenames")

if __name__ == "__main__":
    main()
