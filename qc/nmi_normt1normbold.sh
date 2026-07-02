#!/bin/bash
set -euo pipefail

path_fmriprep="/Users/ga0034de/Desktop/output_priority4_freesurfnobbr"

find "$path_fmriprep" -type f -name "*_desc-coreg_boldref.nii.gz" ! -name "._*" > ~/Desktop/priority_v2511_freesurfernobbr_coreg_boldref_files.txt

mni_priorityyes_nmi_file="$HOME/Desktop/2704priority_freesurfnobbr_v2525_nmi_norm.txt"
wdir="$HOME/Desktop/workdir"
outdir="$HOME/Desktop/freesurfnobbr_bold_t1space"

mkdir -p "$wdir" "$outdir"

mni="/Users/ga0034de/Library/Caches/templateflow/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_desc-brain_T1w.nii.gz"
mni_mask="/Users/ga0034de/Library/Caches/templateflow/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz"

entropy_mni=$(ImageIntensityStatistics 3 "$mni" "$mni_mask" | awk 'NR==2 {print $6}')
echo "Entropy of MNI template: $entropy_mni"
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

echo "sid, session, mattes_wt1_wbold, entropy_wt1, entropy_wbold" > "$mni_priorityyes_nmi_file"

while IFS= read -r file; do
    echo "Processing file: $file"
# if file starts with ./, change it to /Volume/research/LU26D1023-DemonLab/DemonLab/ROSMAP/derivatives/fmriprep/unzipped/
    if [[ "$file" == ./* ]]; then
        file="${path_fmriprep}/${file:2}"
    fi
    
    jobdir=$(echo "$file" | cut -d'/' -f6)
    sid=$(echo "$jobdir" | cut -d'_' -f1)
    session=$(echo "$jobdir" | cut -d'_' -f2)
    echo "$jobdir, Subject: $sid, Session: $session"

    # 1) Try the exact expected path
    fmriprep_dir="${path_fmriprep}/${jobdir}/fmriprep/${sid}/${session}"
    echo "Trying fmriprep dir: $fmriprep_dir"
    # 3) If still missing, remove fmriprep from the path and try again
    if [[ ! -d "$fmriprep_dir" ]]; then
        fmriprep_dir="${path_fmriprep}/${jobdir}/${sid}/${session}"
    fi

    if [[ ! -d "${fmriprep_dir}" ]]; then
        echo "$(date) - WARNING: fmriprep dir not found for ${sid} ${session}"
        exit 0
    fi

    echo "Using fmriprep dir: $fmriprep_dir"

    boldref="$file"
    boldref_mni=$(find_first_file "${fmriprep_dir}/func" "*space-MNI152NLin6Asym_res-2_boldref.nii.gz")
    wt1=$(find_first_file "${fmriprep_dir}/anat" "*space-MNI152NLin6Asym_res-2_desc-preproc_T1w.nii.gz")

    # matrix=$(find "${fmriprep_dir}/func" -type f -name "*from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt")
    # t1=$(find "${fmriprep_dir}/anat" -type f -name "*_desc-preproc_T1w.nii.gz" | grep -v "space-MNI152NLin6Asym_res-2")
    # t1_mask=$(find "${fmriprep_dir}/anat" -type f -name "*_desc-brain_mask.nii.gz" | grep -v "space-MNI152NLin6Asym_res-2")

    # bold_t1space="${outdir}/${sid}_${session}_space-T1w_desc-coreg_boldref.nii.gz"
    # t1_mask_bold_space="${wdir}/mask_${sid}_${session}_bold_space.nii.gz"

    echo "Found files:"
    echo " "

    # antsApplyTransforms -d 3 -i "$boldref" -r "$t1" -o "$bold_t1space" -t "$matrix" --interpolation LanczosWindowedSinc


    # antsApplyTransforms -d 3 -i "$t1_mask" -r "$bold_t1space" -o "$t1_mask_bold_space" -n NearestNeighbor

    # unmasked 
    # in mni space
    mattes_wt1_wbold=$(MeasureImageSimilarity -d 3 -m Mattes["$wt1","$boldref_mni",1,64])
    entropy_wt1=$(ImageIntensityStatistics 3 "$wt1" "$mni_mask" | awk 'NR==2 {print $6}')
    entropy_wbold=$(ImageIntensityStatistics 3 "$boldref_mni" "$mni_mask" | awk 'NR==2 {print $6}')

    # masked 
    # in mni space
    #mask_mattes_wt1_wbold=$(MeasureImageSimilarity -d 3 -m Mattes["$wt1","$boldref_mni",1,64] -x "$mni_mask")
    
    echo "$sid $session done"
    echo "$sid, $session, $mattes_wt1_wbold, $entropy_wt1, $entropy_wbold" | tee -a "$mni_priorityyes_nmi_file"

done < "$HOME/Desktop/priority_v2511_freesurfernobbr_coreg_boldref_files.txt"