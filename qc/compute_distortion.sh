#!/usr/bin/env bash
set -euo pipefail

path_fmriprep="/Users/ga0034de/Desktop/output_priority4_freesurfnobbr"
tmp_dir="$HOME/Desktop/tmp_dir"
distortion_results="$HOME/Desktop/2704freesurfnobbr_distortion_results.txt"

distortion_ratio() {
    local sub_id="$1"
    local ses_id="$2"
    local anat_mask="$3"
    local mask_func="$4"

    # fslmaths "$mask_func" -bin "$func_bin"
    # fslmaths "$anat_mask" -bin "$anat_bin"
    # echo "Binary masks created."
    # Functional voxels outside anatomical brain
    func_outside_anat="$tmp_dir/${sub_id}_${ses_id}_func_outside_anat.nii.gz"
    fslmaths "$mask_func" -sub "$anat_mask" -thr 0 -bin "$func_outside_anat"

    # Volumes
    vol_outside=$(fslstats "$func_outside_anat" -V | awk '{print $2}')
    vol_func=$(fslstats "$mask_func" -V | awk '{print $2}')

    # Ratio
    distortion_val=$(python3 -c "print(round($vol_outside/$vol_func,4))")
    echo "$sub_id, $ses_id, $distortion_val" >> "$distortion_results"
}


echo "sub_id, ses_id, distortion_val" > "$distortion_results"

export -f distortion_ratio

for sub_ses in $(cat "$HOME/Desktop/priority_v2511_freesurfernobbr_coreg_boldref_files.txt"); do
    sub_id=$(grep -o 'sub-[0-9]\+' <<< "$sub_ses" | head -1)
    ses_id=$(grep -o 'ses-[0-9]\+' <<< "$sub_ses" | head -1)
    find "$path_fmriprep" -type f -name "anat/${sub_id}_${ses_id}.*space-MNI.*desc-brain_mask.nii.gz"

    anat_mask=$(find "$path_fmriprep" -type f -name "${sub_id}_${ses_id}*space-MNI*desc-brain_mask.nii.gz" ! -name "*task-rest*")
    mask_func=$(find "$path_fmriprep" -type f -name "${sub_id}_${ses_id}_task-rest*space-MNI*desc-brain_mask.nii.gz")
    
    echo "Computing distortion for $sub_id $ses_id with anat mask: $anat_mask and func mask: $mask_func"
    distortion_ratio $sub_id $ses_id $anat_mask $mask_func
done