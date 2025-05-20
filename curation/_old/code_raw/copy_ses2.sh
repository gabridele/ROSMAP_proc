#!/bin/bash

# Set the base directory (change this to the target directory)
SOURCE="../sourcedata/batch_1"

# Loop through all directories starting with "sub" in the base directory
for sub_dir in "$SOURCE"/sub*; do
    # Ensure it's a directory
    if [[ -d "$sub_dir" ]]; then
        sub=$(basename "$sub_dir")
        # Loop through subdirectories inside the "sub-" directory
        rsync -av --copy-links --exclude=".*" "$sub_dir"/* "$sub"/
        datalad save -m 'copied sessions from sourcedata batch_1' -d $sub_dir
    fi
done