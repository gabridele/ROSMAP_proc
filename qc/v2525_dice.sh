#!/bin/bash
set -euo pipefail

path_fmriprep="/Users/ga0034de/Desktop/output_priority4_freesurfnobbr"
tmp_dir="$HOME/Desktop/tmp_dir"
dice_results="$HOME/Desktop/2704v2525freesurfnobbr_dice_results.txt"

echo "sub_id, ses_id, dice_val" > "$dice_results"
# ===============================
# Dice (FSL)
# ===============================

compute_dice() {
    local sub_id="$1"
    local ses_id="$2"
    local anat_mask="$3"
    local mask_func="$4"

    intersection="$tmp_dir/${sub_id}_${ses_id}_mask_intersection.nii.gz"

    anat=$(fslstats "$anat_mask" -V | awk '{print $1}')
    func=$(fslstats "$mask_func" -V | awk '{print $1}')
    echo "Anat voxels: $anat, Func voxels: $func"
    echo "Computing intersection for $sub_id $ses_id"
    fslmaths "$anat_mask" -mul "$mask_func" "$intersection"
    echo "created intersection file: $intersection"
    I=$(fslstats "$intersection" -V | awk '{print $1}')

    dice_val=$(python3 -c "print(round(2*$I/($anat+$func),3))")
    echo "$sub_id, $ses_id, $dice_val" >> "$dice_results"
}

export -f compute_dice

for sub_ses in $(cat "$HOME/Desktop/priority_v2511_freesurfernobbr_coreg_boldref_files.txt"); do
    sub_id=$(grep -o 'sub-[0-9]\+' <<< "$sub_ses" | head -1)
    ses_id=$(grep -o 'ses-[0-9]\+' <<< "$sub_ses" | head -1)
    find "$path_fmriprep" -type f -name "anat/${sub_id}_${ses_id}.*space-MNI.*desc-brain_mask.nii.gz"

    anat_mask=$(find "$path_fmriprep" -type f -name "${sub_id}_${ses_id}*space-MNI*desc-brain_mask.nii.gz" ! -name "*task-rest*")
    mask_func=$(find "$path_fmriprep" -type f -name "${sub_id}_${ses_id}_task-rest*space-MNI*desc-brain_mask.nii.gz")
    
    echo "Computing Dice for $sub_id $ses_id with anat mask: $anat_mask and func mask: $mask_func"
    compute_dice $sub_id $ses_id $anat_mask $mask_func
done
