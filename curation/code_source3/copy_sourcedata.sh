#!/bin/bash
# This script copies sourcedata from the source directory to the target directory
# for each subject listed in the subs_batch3.txt file.
# It uses datalad to create a new dataset for each subject and saves the changes after copying files into it.
N=1
source_base="/home/gabridele/backup/ROSMAP/batch_3/RIRC"
target_base="/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_3"

copy_source_to_target() {
    IFS=',' read -r target_id middle_digits <<< "$1"

    # Construct absolute paths for source and target
    source_folder=$(find "$source_base" -maxdepth 2 -type d -name "*_${middle_digits}_${target_id#sub-}" -print -quit)
    target_folder="$target_base/$target_id"

    if [ -d "$source_folder" ] && [ ! -e "$target_folder" ]; then
    
    	datalad create -d "$target_folder"
    	
        cp -r "$source_folder" "$target_folder/"
        echo "Copied $source_folder to $target_folder/"

        datalad save -m "populated subject directory with sourcedata" -d "$target_folder"
    else
        echo "Skipping: $source_folder or $target_folder does not exist"
    fi
}
i=0
(
while IFS= read -r line; do
    ((i=i%N)); ((i++==0)) && wait
    copy_source_to_target "$line" &
done < "code/subs_batch3.txt"
wait
)
