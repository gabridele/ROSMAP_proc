#!/bin/bash
set -euo pipefail

# module load c3d/1.4.4
# module load fsl/6.0
# source ${FSLDIR}/etc/fslconf/fsl.sh

pathroot="/Users/ga0034de/Desktop/output_priority4_freesurfnobbr"

list_sid="/Users/ga0034de/Documents/R_projs/priority_rosmap/ID_list.csv"
output_file="/Users/ga0034de/Desktop/2704nobbr_dropout10_new.txt"
log_file="/Users/ga0034de/Desktop/dropout10_errors.log"

echo "sid, session, volume_gm, nvox_gm, intensity_gm, volume_dropout, nvox_dropout, intensity_dropout" > "$output_file"
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
    fmriprep_dir="${pathroot}/${sid}_${session}_fmriprep-25-2-5/fmriprep/${sid}/${session}"

    # 3) If still missing, try to locate any directory under derivatives that looks like the subject/session
    if [[ ! -d "$fmriprep_dir" ]]; then
        fmriprep_dir="${pathroot}/${sid}_${session}_fmriprep-25-2-5/${sid}/${session}"
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

    echo "Creating GM binary mask (thr=0.3) at: $mask_gm_thr"
    fslmaths "$gm_seg" -thr 0.3 -bin "$mask_gm_thr"

    # merge anat+func masks (c3d add, replace) — use both masks
    c3d "$mask_anat" "$mask_func" -add -replace 2 1 -o "$mask_merged"
	echo "Merged anat+func mask created at: $mask_merged"
    # mask the boldref with the merged mask
    fslmaths "$boldref" -mul "$mask_merged" "$boldref_masked"

    # compute threshold from masked boldref and create new mask
    thresh=$(fslstats "$boldref_masked" -l 0.001 -P 10 2>/dev/null | awk '{print $1}')
	echo "DEBUG thresh raw: '$(fslstats "$boldref_masked" -l 0.001 -P 10)'"
	echo "DEBUG thresh var: '$thresh'"
    fslmaths "$boldref_masked" -thr "$thresh" -bin "$new_mask_func"

    # compute dropout: GM mask minus new func mask (inverting new_mask_func before add by -scale -1)
    c3d "$mask_gm_thr" "$new_mask_func" -scale -1 -add -o "$mask_dropout"
    # convert >1 values -> 0 and 1 stays 1
    c3d "$mask_dropout" -replace 1 1 0 0 -1 0 -o "$mask_dropout"

    # stats
    vol_dropout=$(c3d "$mask_dropout" -dup -lstat | awk 'NR==3 {print $7}')
    nvox_dropout=$(c3d "$mask_dropout" -dup -lstat | awk 'NR==3 {print $6}')

    vol_gm=$(c3d "$mask_gm_thr" -dup -lstat | awk 'NR==3 {print $7}')
    nvox_gm=$(c3d "$mask_gm_thr" -dup -lstat | awk 'NR==3 {print $6}')

    fslmaths "$mask_gm_thr" -sub "$mask_dropout" "$mask_gm_thr_clean"
    intensity_gm=$(fslstats "$boldref" -k "$mask_gm_thr_clean" -M)
    intensity_dropout=$(fslstats "$boldref" -k "$mask_dropout" -M)

    echo "$sid $session | GM: $vol_gm $intensity_gm | dropouts: $vol_dropout $intensity_dropout"
    echo "$sid, $session, $vol_gm, $nvox_gm, $intensity_gm, $vol_dropout, $nvox_dropout, $intensity_dropout" >> "$output_file"

    # cleanup temp files (don't fail if missing)
    rm -f "$mask_gm_thr_clean" "$mask_merged" "$boldref_masked"

done < "$list_sid"