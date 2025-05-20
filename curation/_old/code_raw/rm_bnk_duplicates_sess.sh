#!/bin/bash
# code to unlock and remove the duplicates of the listed subjects
# code to be run in the directory where the subjects are located
set -e
cwd=$(pwd)

folders=("sub-73177635" "sub-15938020" "sub-65952361" "sub-52052115" "sub-83216408" "sub-78010047")
for folder in "${folders[@]}"; do
    if [ ! -d "$folder" ]; then
        echo "Directory $folder does not exist. Please check the path."
        exit 1
    fi

    datalad unlock ${folder}/ses-0/anat/${folder}_ses-0_acq-20090211*
    rm ${folder}/ses-0/anat/${folder}_ses-0_acq-20090211*

done

folders=("sub-24644776" "sub-48480640" "sub-18455382" "sub-44842532" "sub-49094578")
for folder in "${folders[@]}"; do
    if [ ! -d "$folder" ]; then
        echo "Directory $folder does not exist. Please check the path."
        exit 1
    fi
    datalad unlock ${folder}/ses-1/anat/${folder}_ses-1_acq-20090211*
    rm ${folder}/ses-1/anat/${folder}_ses-1_acq-20090211*
done