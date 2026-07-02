#!/bin/bash
set -euo pipefail

# module load c3d/1.4.4
# module load fsl/6.0
# source ${FSLDIR}/etc/fslconf/fsl.sh

pathroot="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP"

list_sid="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP/raw/code/ID_list.csv"
output_file="/Users/ga0034de/Desktop/dropout10_old.txt"
log_file="/Users/ga0034de/Desktop/dropout10_old_errors.log"

echo "sid session volume_gm nvox_gm intensity_gm volume_dropout nvox_dropout intensity_dropout" > "$output_file"
: > "$log_file"

wdir="/Users/ga0034de/Desktop/workdir"

# helper: find first file matching pattern
find_first_file() {
  local base_dir="$1"; shift
  local pattern="$*"
  find "$base_dir" -type f -name "$pattern" -print -quit 2>/dev/null || true
}

# helper: find first directory matching a -path pattern
find_first_dir() {
  local base_dir="$1"; shift
  local path_pattern="$*"
  find "$base_dir" -type d -path "$path_pattern" -print -quit 2>/dev/null || true
}

# Read CSV and iterate rows
while IFS=',' read -r col1 col2 rest; do
    # Skip header / empty lines
    if [[ "$col1" == "sub_id" || -z "${col1:-}" ]]; then
        continue
    fi

    sid="$col1"
    session="$col2"
    echo "Processing subject: $sid, session: $session"

    # 1) Try the exact expected path (old style you used)
    fmriprep_dir="${pathroot}/derivatives/fmriprep/unzipped/${sid}_${session}_fmriprep-25-1-1/fmriprep/${sid}/${session}"

    # 3) If still missing, try to locate any directory under derivatives that looks like the subject/session
    if [[ ! -d "$fmriprep_dir" ]]; then
        fmriprep_dir="${pathroot}/derivatives/fmriprep/unzipped/${sid}_${session}_fmriprep-25-1-1/${sid}/${session}"
    fi

    if [[ ! -d "${fmriprep_dir}" ]]; then
        echo "$(date) - WARNING: fmriprep dir not found for ${sid} ${session}" | tee -a "$log_file"
        # log and skip this subject
        echo "$sid $session missing_fmriprep" >> "$output_file"
        continue
    fi

    echo "Using fmriprep dir: $fmriprep_dir"

    # find required files (take first match only)
    mask_func=$(find_first_file "${fmriprep_dir}/func" "*space-MNI152NLin6Asym_res-2_desc-brain_mask.nii.gz")
    mask_anat=$(find_first_file "${fmriprep_dir}/anat" "*space-MNI152NLin6Asym_res-2_desc-brain_mask.nii.gz")
    boldref=$(find_first_file "${fmriprep_dir}/func" "*space-MNI152NLin6Asym_res-2_boldref.nii.gz")
    bold=$(find_first_file "${fmriprep_dir}/func" "*space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz")
    gm_seg=$(find_first_file "${fmriprep_dir}/anat" "*space-MNI152NLin6Asym_res-2_label-GM_probseg.nii.gz")

    echo "Found files:"
    echo " mask_func: ${mask_func:-NONE}"
    echo " mask_anat: ${mask_anat:-NONE}"
    echo " boldref:   ${boldref:-NONE}"
    echo " bold:      ${bold:-NONE}"
    echo " gm_seg:    ${gm_seg:-NONE}"

    # check required files exist before processing
    missing=0
    for f in mask_func mask_anat boldref gm_seg; do
        val="${!f:-}"
        if [[ -z "$val" ]]; then
            echo "$(date) - ERROR: required file '$f' missing for ${sid} ${session}" | tee -a "$log_file"
            missing=1
        fi
    done

    if [[ $missing -eq 1 ]]; then
        echo "$sid $session missing_files" >> "$output_file"
        continue
    fi

    # prepare working dirs
    mkdir -p "$wdir/${sid}_${session}/func"
    mkdir -p "$wdir/${sid}_${session}/anat"

    boldref_masked="$wdir/${sid}_${session}/func/${sid}_${session}_task-rest_space-MNI152NLin6Asym_res-2_boldref_masked.nii.gz"
    mask_gm_thr="$wdir/${sid}_${session}/anat/${sid}_${session}_space-MNI152NLin6Asym_res-2_label-GM_mask-03thrbin.nii.gz"
    new_mask_func="$wdir/${sid}_${session}/func/${sid}_${session}_space-MNI152NLin6Asym_res-2_boldref_new-mask.nii.gz"
    mask_dropout="$wdir/${sid}_${session}/func/${sid}_${session}_space-MNI152NLin6Asym_res-2_dropout10_mask.nii.gz"
    mask_gm_thr_clean="$wdir/${sid}_${session}/anat/${sid}_${session}_space-MNI152NLin6Asym_res-2_GM_wo-dropout_temp.nii.gz"
    mask_merged="$wdir/${sid}_${session}/func/${sid}_${session}_space-MNI152NLin6Asym_res-2_anat-func-masks-merged.nii.gz"

    mask_dropout_old="$wdir/${sid}_${session}/func/${sid}_${session}_space-MNI152NLin6Asym_res-2_desc-dropoutsmaskold.nii.gz"

    echo "Creating GM binary mask (thr=0.3) at: $mask_gm_thr"
    fslmaths "$gm_seg" -thr 0.3 -bin "$mask_gm_thr"

    fslmaths "$mask_gm_thr" -sub "$mask_func" "$mask_dropout_old"
    fslmaths "$mask_dropout_old" -thr 0 -bin "$mask_dropout_old"
    
    # stats
    vol_dropout=$(c3d "$mask_dropout_old" -dup -lstat | awk 'NR==3 {print $7}')
    nvox_dropout=$(c3d "$mask_dropout_old" -dup -lstat | awk 'NR==3 {print $6}')

    vol_gm=$(c3d "$mask_gm_thr" -dup -lstat | awk 'NR==3 {print $7}')
    nvox_gm=$(c3d "$mask_gm_thr" -dup -lstat | awk 'NR==3 {print $6}')

    fslmaths "$mask_gm_thr" -sub "$mask_dropout_old" "$mask_gm_thr_clean"
    intensity_gm=$(fslstats "$boldref" -k "$mask_gm_thr_clean" -M)
    intensity_dropout=$(fslstats "$boldref" -k "$mask_dropout_old" -M)

    echo "$sid $session | GM: $vol_gm $intensity_gm | dropouts: $vol_dropout $intensity_dropout"
    echo "$sid $session $vol_gm $nvox_gm $intensity_gm $vol_dropout $nvox_dropout $intensity_dropout" >> "$output_file"

    # cleanup temp files (don't fail if missing)
    rm -f "$mask_gm_thr_clean" "$mask_merged" "$boldref_masked"

done < "$list_sid"