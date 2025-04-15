#!/bin/bash

# Set the base directory (change this to the target directory)
BASE_DIR="."

# Loop through all directories starting with "sub" in the base directory
for sub_dir in "$BASE_DIR"/sub*; do
    # Ensure it's a directory
    if [[ -d "$sub_dir" ]]; then
        # Loop through subdirectories inside the "sub-" directory
        for ses_dir in "$sub_dir"/ses-*; do
            # Ensure it's a directory before removing
            if [[ -d "$ses_dir" ]]; then
                echo "Removing: $ses_dir"
                rm -rf "$ses_dir"
            fi
        done
    fi
done