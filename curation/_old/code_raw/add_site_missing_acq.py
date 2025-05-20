import os
import glob
import re
import concurrent.futures
import datalad.api as dl

def modify_filenames(directory):
    """
    Add site info to acq- sections missing it
    """

    #print(f'Processing {directory}')
    # Read the string from ../anat
    # anat_path = os.path.join(os.path.dirname(directory), 'anat')
    # anat_string = ''
    # for file in os.listdir(anat_path):
    #     if file.endswith('.json'):
    #         match = re.search(r'acq-([a-zA-Z]+)', file)
    #         if match:
    #             anat_string = match.group(1)
    #             anat_string = ''.join(re.findall(r'[a-zA-Z]+', anat_string))
    #         break
    
    # Iterate through files in the directory
    for filename in os.listdir(directory):
        if 'ProtUnk' in filename:
            continue
        if 'RIRC' in filename:
            continue
        # Check if the filename contains no site in 'acq-' section
        match = re.search(r'acq-(\d+)', filename)
        #print(match)
        if match:
            substring = match.group(1)
            start = match.start(1)
            end = match.end(1)
            print('substring:', substring)
            print('start:', start)
            print('end:', end)
            if '20090211' in substring:
                new_filename = filename[:start] + 'BNK' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '20120501' in substring:
                new_filename = filename[:start] + 'MG' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '201606' in substring:
                new_filename = filename[:start] + 'MG' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '20120221' in substring:
                new_filename = filename[:start] + 'UC' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '20140922' in substring:
                new_filename = filename[:start] + 'UC' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '20150706' in substring:
                new_filename = filename[:start] + 'UC' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '20151120' in substring:
                new_filename = filename[:start] + 'UC' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)
            if '20160125' in substring:
                new_filename = filename[:start] + 'UC' + substring + filename[end:]
            # Extract the substring matched by the regex
                print('new_filename:', new_filename)

            # If modified, rename the file
            filename_path = os.path.join(directory, filename)
            dl.unlock(filename_path)
            print(f'Renaming {filename} to {new_filename}')
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))

def rename_RIRC(directory):
    print(f'Processing {directory}')
    
    # Iterate through files in the directory
    for filename in os.listdir(directory):
        if 'RIRC' not in filename:
            continue
        if 'sub-13044513_ses-0_acq-RIRC' in filename:
            new_filename = filename.replace('RIRC', 'RIRC20')
        if 'sub-70876731_ses-0_acq-RIRC' in filename:
            new_filename = filename.replace('RIRC', 'RIRC20')
        else:
            match = re.search(r'acq-([^-_]+)', filename)
            if match:
                # Extract the substring matched by the regex
                substring = match.group(1)
                if 'RIRC23' in substring:
                    new_substring = substring.replace('RIRC23', 'RIRC2023')[:8] + '0803' + substring[12:]
                else:
                    new_substring = substring[:8] + '0803' + substring[12:]
                start = match.start(1)
                end = match.end(1)
                # Rename the file
                new_filename = filename[:start] + new_substring + filename[end:]
                
        filename_path = os.path.join(directory, filename)
        dl.unlock(filename_path)
        print(f'Renaming {filename} to {new_filename}')
        os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))


def main():
    # Get the base directory
    base_directory = os.getcwd()

    # Collect all sub-*/ses-*/anat directories
    anat_dirs = sorted(glob.glob(os.path.join(base_directory, 'sub-*', 'ses-*', 'anat')))
    rirc_dirs = sorted(glob.glob(os.path.join(base_directory, 'sub-*', 'ses-*', '*')))

    # Run the function in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(modify_filenames, anat_dirs)
    
    # Wait for the first parallel execution to complete
    executor.shutdown(wait=True)

    # Run the second function in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(rename_RIRC, rirc_dirs)

if __name__ == '__main__':
    main()