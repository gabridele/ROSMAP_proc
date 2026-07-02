#!/bin/bash
set -euo pipefail

# https://github.com/ANTsX/ANTs/discussions/1706

path_fmriprep="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP/derivatives/fmriprep/unzipped"

list_sid="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP/raw/code/ID_list.csv"
output_file="/Users/ga0034de/Desktop/nmi_masked.txt"
log_file="/Users/ga0034de/Desktop/nmi_masked_errors.log"

nmi_file="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP/derivatives/code/nmi_masked.txt"

mni="/Users/ga0034de/Library/Caches/templateflow/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_desc-brain_T1w.nii.gz"
mni_mask="/Users/ga0034de/Library/Caches/templateflow/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz"

entropy_mni=$(ImageIntensityStatistics 3 "$mni" "$mni_mask" | awk 'NR==2 {print $6}')

echo "sid session mattes_t1_bold mattes_wt1_mni mattes_wbold_mni entropy_t1 entropy_bold entropy_wt1 entropy_wbold entropy_mni" > "$nmi_file"

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
    fmriprep_dir="${path_fmriprep}/${sid}_${session}_fmriprep-25-1-1/fmriprep/${sid}/${session}"

    # 3) If still missing, try to locate any directory under derivatives that looks like the subject/session
    if [[ ! -d "$fmriprep_dir" ]]; then
        fmriprep_dir="${path_fmriprep}/${sid}_${session}_fmriprep-25-1-1/${sid}/${session}"
    fi

    if [[ ! -d "${fmriprep_dir}" ]]; then
        echo "$(date) - WARNING: fmriprep dir not found for ${sid} ${session}" | tee -a "$log_file"
        # log and skip this subject
        echo "$sid $session missing_fmriprep" >> "$output_file"
        continue
    fi

    echo "Using fmriprep dir: $fmriprep_dir"


    boldref_mni=$(find_first_file "${fmriprep_dir}/func" "*space-MNI152NLin6Asym_res-2_boldref.nii.gz")
    t1=$(find "${fmriprep_dir}/anat" -type f -name "*_desc-preproc_T1w.nii.gz" | grep -v "space-MNI152NLin6Asym_res-2")
    t1_mask=$(find "${fmriprep_dir}/anat" -type f -name "*_desc-brain_mask.nii.gz" | grep -v "space-MNI152NLin6Asym_res-2")
    boldref=$(find "${fmriprep_dir}/func" -type f -name "*desc-coreg_boldref.nii.gz" | grep -v "space-MNI152NLin6Asym_res-2")
    echo "Found files: $t1_mask"
    t1_mask_bold_space="${wdir}/mask_${sid}_${session}_bold_space.nii.gz"
    t1_resampled="${wdir}/t1_${sid}_${session}_bold_space.nii.gz"

    wt1=$(find_first_file "${fmriprep_dir}/anat" "*space-MNI152NLin6Asym_res-2_desc-preproc_T1w.nii.gz")

    # need masks to be in the same resolution for calculating entropy + then resample t1 to make it homogeneous
    antsApplyTransforms -d 3 -i "$t1_mask" -r "$boldref" -o "$t1_mask_bold_space" -n NearestNeighbor
    antsApplyTransforms -d 3 -i "$t1" -r "$boldref" -o "$t1_resampled" -n BSpline

    # Mattes
    mattes_t1_bold=$(MeasureImageSimilarity -d 3 -m Mattes["$t1_resampled","$boldref",1,64] -x "$t1_mask_bold_space")
    mattes_wt1_mni=$(MeasureImageSimilarity -d 3 -m Mattes["$wt1","$mni",1,64] -x "$mni_mask")
    mattes_wbold_mni=$(MeasureImageSimilarity -d 3 -m Mattes["$boldref_mni","$mni",1,64] -x "$mni_mask")
    echo "$wt1 $mni $mni_mask"
    # Entropies
    entropy_t1=$(ImageIntensityStatistics 3 "$t1_resampled" "$t1_mask_bold_space" | awk 'NR==2 {print $6}')
    entropy_bold=$(ImageIntensityStatistics 3 "$boldref" "$t1_mask_bold_space" | awk 'NR==2 {print $6}')
    entropy_wt1=$(ImageIntensityStatistics 3 "$wt1" "$mni_mask" | awk 'NR==2 {print $6}')
    entropy_wbold=$(ImageIntensityStatistics 3 "$boldref_mni" "$mni_mask" | awk 'NR==2 {print $6}')

    echo "$sid $session $mattes_t1_bold $mattes_wt1_mni $mattes_wbold_mni"
    echo "$sid $session $mattes_t1_bold $mattes_wt1_mni $mattes_wbold_mni $entropy_t1 $entropy_bold $entropy_wt1 $entropy_wbold $entropy_mni" >> "$nmi_file"
    echo sid session mattes_t1_bold mattes_wt1_mni mattes_wbold_mni entropy_t1 entropy_bold entropy_wt1 entropy_wbold entropy_mni
done < "$list_sid"