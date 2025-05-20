import os
import re
import datalad.api as dl
from multiprocessing import Pool
import subprocess

def find_files_in_dir(args):
    """ Helper function to find files with 'rfMRI' in a given directory """
    directory, files = args
    return [os.path.join(directory, file) for file in files if 'rfMRI' in file]

def find_rfmri_files(directory):
    """ Parallelized function to find all files containing 'rfMRI' in their name """
    file_paths = []
    
    with Pool() as pool:
        results = pool.map(find_files_in_dir, [(root, files) for root, _, files in os.walk(directory)])
    
    for result in results:
        file_paths.extend(result)
    
    return file_paths
    
def unlock_files(file_paths):
    """ Unlock files using datalad """
    for file_path in file_paths:
        subprocess.run(["datalad", "unlock", file_path], check=True)

def rename_rfMRI_file(old_path):
    """ Rename acq- string section by replacing 'rfMRI_' """
    dir_name, file_name = os.path.split(old_path)
    
    # Replace 'rfMRI_' with 'EPI
    new_name = file_name.replace('rfMRI', 'EPI')
    new_path = os.path.join(dir_name, new_name)
    
    # Rename the file
    os.rename(old_path, new_path)
    print(f'Renamed: {old_path} -> {new_path}')
    
def main():
    directory = "."
    file_paths = find_rfmri_files(directory)
    print('proceeding with unlocking of paths')
    unlock_files(file_paths)
    
    print('now stripping rfMRI from filename')
    # Rename files containing 'rfMRI'
    with Pool() as pool:
        pool.map(rename_rfMRI_file, file_paths)

if __name__ == "__main__":
    main()
