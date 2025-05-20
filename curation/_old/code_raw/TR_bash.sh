#!/bin/bash

# Check current working directory
pwd

# Check if the file 'code/funcs_TR.txt' does not exist
if [[ ! -f code/funcs_TR.txt ]]; then
    # Find all symbolic links to files with names ending in '*bold.nii.gz', sort them, and save the list to 'code/funcs_TR.txt'
    find . -type l -name '*bold.nii.gz' | sort > code/funcs_TR.txt
    
    datalad save -m 'created list of all func files to make sure TR matches in json and header' code/funcs_TR.txt
fi

# Run the Python script 'change_TR.py' with 'code/funcs_TR.txt' as an argument
python code/change_TR.py "code/funcs_TR.txt"
