#!/bin/bash
N=1
# dont parallelize because datalad doesnt support it
source_base="/home/gabridele/backup/ROSMAP/batch_2/BNK/090211"
target_base="/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2"

copy_source_to_target() {
    IFS=',' read -r target_id middle_digits <<< "$1"

    # Construct absolute paths for source and target
    source_folder=$(find "$source_base" -maxdepth 1 -type d -name "*_${middle_digits}_${target_id#sub-}" -print -quit)
    target_folder="$target_base/$target_id"

    if [ -d "$source_folder" ] && [ -d "$target_folder" ]; then
        cp -r "$source_folder" "$target_folder/"
        echo "Copied $source_folder to $target_folder/"
        cd "$target_folder"
        datalad save -m "populated subject directory with sourcedata"
	cd ..
    else
        echo "Skipping: $source_folder or $target_folder does not exist"
    fi
}
i=0
(
while IFS= read -r line; do
    ((i=i%N)); ((i++==0)) && wait
    copy_source_to_target "$line" &
done < "code/subs_batch2.txt"
wait
)
