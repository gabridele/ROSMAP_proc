#!/bin/bash

# Define the source directory containing the data to be processed
SOURCE_DIR="../sourcedata/batch_2"

# Define the destination directory where the processed files will be copied
DEST_DIR="."

# Loop through each subject directory in the source directory
for sub in "$SOURCE_DIR"/sub-*; do
    # Skip if the item is not a directory
    [ -d "$sub" ] || continue

    # Loop through each session directory within the subject directory
    for ses_dir in "$sub"/*; do
        # Skip if the item is not a directory
        [ -d "$ses_dir" ] || continue
        
        # Get a list of all files in the session directory, excluding hidden files
        files=($(ls "$ses_dir" | grep -v "^\."))
        
        # Loop through each file in the session directory
        for file in "${files[@]}"; do
            # Check if the file name contains a session identifier (_ses-<number>_)
            if [[ "$file" =~ _ses-([0-9]+)_ ]]; then
                # Extract the session number and construct the destination session directory
                ses_num="ses-${BASH_REMATCH[1]}"
                dest_ses="$DEST_DIR/$(basename "$sub")/$ses_num"
                mkdir -p "$dest_ses"  # Create the destination session directory if it doesn't exist
                
                # Determine the subfolder based on the file type
                case "$file" in
                    *T1w.*) subfolder="anat" ;;       # Files with "T1w" go to the "anat" subfolder
                    *phase_map*) subfolder="fmap" ;;  # Files with "phase_map" go to the "fmap" subfolder
                    *dwi.*) subfolder="dwi" ;;        # Files with "dwi" go to the "dwi" subfolder
                    *task-rest*) subfolder="func" ;;  # Files with "task-rest" go to the "func" subfolder
                    *) subfolder="" ;;                # Other files are ignored
                esac
                
                # If a subfolder is determined, copy the file to the appropriate location
                if [[ -n "$subfolder" ]]; then
                    mkdir -p "$dest_ses/$subfolder"  # Create the subfolder if it doesn't exist
                    rsync -av --copy-links --exclude=".*" "$ses_dir/$file" "$dest_ses/$subfolder/"  # Copy the file
                fi
            fi
        done

    done
done