import os
import re
import datalad.api as dl
from multiprocessing import Pool
import subprocess
from concurrent.futures import ThreadPoolExecutor

def find_files_in_dir(args):
    """ Helper function to find files with 'ADC' in a given directory """
    directory, files = args
    return [os.path.join(directory, file) for file in files if 'ADC' in file]

def find_ADC_files(directory):
    """ Parallelized function to find all files containing 'ADC' in their name """
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

def rename_ADC_file(old_path):
    """ 
    Rename files by replacing 'dwi' with 'ADC' suffix (_ADC.nii.gz)
    
    Args:
        old_path (str): Path to the file to be renamed.
    """
    dir_name, file_name = os.path.split(old_path)
    
    new_name = file_name.replace('ADC', '').replace('dwi', 'ADC')
    new_path = os.path.join(dir_name, new_name)
    
    # Rename the file
    os.rename(old_path, new_path)
    print(f'Renamed: {old_path} -> {new_path}')
    
def main():
    
    directory = "."  # Change this to your target directory
    file_paths = find_ADC_files(directory)
    print('proceeding with unlocking of paths')
    unlock_files(file_paths)
    
    print('now stripping ADC from filename')
    # Rename files containing 'ADC'

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        executor.map(rename_ADC_file, file_paths)

if __name__ == "__main__":
    main()