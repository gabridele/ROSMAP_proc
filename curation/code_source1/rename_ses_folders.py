import os
import re
import datalad.api as dl
### script not used
def get_session_number(filename):
    """Extracts the ses- number from a filename."""
    match = re.search(r'ses-(\d+)', filename)
    return match.group(1) if match else None

def rename_mismatched_folders(base_path):
    """Navigates sub-* directories and verifies ses-* folder names."""
    for subdir in os.listdir(base_path):
        subdir_path = os.path.join(base_path, subdir)
        
        if os.path.isdir(subdir_path) and subdir.startswith("sub-"):
            for ses_folder in os.listdir(subdir_path):
                ses_path = os.path.join(subdir_path, ses_folder)
                
                if os.path.isdir(ses_path) and ses_folder.startswith("ses-"):
                    for inner_folder in os.listdir(ses_path):
                        inner_path = os.path.join(ses_path, inner_folder)
                        
                        if os.path.isdir(inner_path):
                            for file in os.listdir(inner_path):
                                ses_number = get_session_number(file)
                                if ses_number and ses_folder != f"ses-{ses_number}":
                                    new_ses_path = os.path.join(subdir_path, f"ses-{ses_number}")
                                    print(f"Found mismatched folder: {ses_path}")
                                    print(f"Renaming {ses_path} -> {new_ses_path}")
                                    # raise Exception("This is a dry run. Remove the raise statement to rename folders.")
                                    dl.unlock(dataset=subdir, path=ses_path)
                                    os.rename(ses_path, new_ses_path)
                                    dl.save(dataset=subdir, message=f"Rename {ses_folder} to ses-{ses_number}")
                                    break  # Rename once per folder, assuming consistency

if __name__ == "__main__":
    base_directory = os.getcwd()  # Current directory
    rename_mismatched_folders(base_directory)