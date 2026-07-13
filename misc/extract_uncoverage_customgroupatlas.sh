#!/bin/bash
set -euo pipefail

proj_dir=$1
tmpdir=$2
outdir=$3

ATLAS="/Users/ga0034de/Desktop/atlas-Schaefer400TianS2Cereb-space-MNI152NLin6Asym/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz"

LABELS="/Users/ga0034de/Desktop/atlas-Schaefer400TianS2Cereb-space-MNI152NLin6Asym/atlas-Schaefer400TianS2Cereb_dseg.tsv"

extract_ts() {
    local sub_ses=$1
    local denoised_bold
    local output

    #base=$(basename "$dir")
    # strip sub_ses of ""
    sub_ses=${sub_ses#\"}
    sub_ses=${sub_ses%\"}
    echo "##########################"

    echo "Subject ID: $sub_ses"

    if [[ ! -f "$LABELS" ]]; then
        echo "ERROR: Missing labels file: $LABELS" >&2
        return 1
    fi

    denoised_bold=$(find "${sub_ses}_xcpd-0-11-1/xcp_d_nifti/${sub_ses%_ses*}/${sub_ses##*_}/func" -type f -name "${sub_ses}_task-rest_acq-*_space-MNI152NLin6Asym_res-2_desc-denoised_bold.nii.gz")

    if [[ -z "$denoised_bold" ]]; then
        echo "ERROR: Could not find denoised BOLD file" >&2
        return 1
    fi

    echo "Found BOLD file: $denoised_bold"

    mkdir -p "$outdir"
    output="$outdir/${sub_ses}_task-rest_space-MNI152NLin6Asym_seg-Schaefer400TianS2Cereb17Networks_timeseries_desc-nocoverage.tsv"

    #individual_atlas="/Users/ga0034de/Desktop/individual_atlases/${sid}_${session}_400atlases/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz"

    python /Users/ga0034de/Desktop/code_tosort/extract_uncoverage.py "$ATLAS" "$LABELS" "$denoised_bold" "$output"

    echo "Time series extracted for $sub_ses"
    echo "Output: $output"
    echo "##########################"
}

export -f extract_ts
export LABELS ATLAS outdir
# sort dirs by subject and session
cat "/Users/ga0034de/Desktop/rosmap_IDlist_subses.csv" | parallel -j 2 extract_ts {}