#!/bin/bash
set -euo pipefail

proj_dir=$1
tmpdir=$2
outdir=$3

#ATLAS="/Users/ga0034de/Library/Caches/templateflow/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_atlas-Schaefer2018_desc-400Parcels17Networks_dseg.nii.gz"

LABELS="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP/derivatives/atlases/atlas-Schaefer400TianS2Cereb/MNINLin6Asym/atlas-Schaefer400TianS2Cereb_dseg.tsv"

extract_ts() {
    local base=$1
    local base
    local sid
    local session
    local denoised_bold
    local output

    #base=$(basename "$dir")

    echo "##########################"
    echo "Processing dir: $dir"

    sid=$(grep -oE 'sub-[^_]+' <<< "$base" || true)
    session=$(grep -oE 'ses-[^_]+' <<< "$base" || true)

    if [[ -z "$sid" || -z "$session" ]]; then
        echo "ERROR: Could not determine subject/session from directory name: $base" >&2
        return 1
    fi

    echo "Subject ID: $sid"
    echo "Session: $session"

    if [[ ! -f "$LABELS" ]]; then
        echo "ERROR: Missing labels file: $LABELS" >&2
        return 1
    fi

    denoised_bold=$1
    #(
    #    find "${sid}_${session}_xcpd-0-11-1/xcp_d_nifti/${sid}/${session}/func" -type f \
     #       -name "${sid}_${session}_task-rest_acq-*_space-MNI152NLin6Asym_res-2_desc-denoised_bold.nii.gz"
    #)

    if [[ -z "$denoised_bold" ]]; then
        echo "ERROR: Could not find denoised BOLD file in $dir" >&2
        return 1
    fi

    echo "Found BOLD file: $denoised_bold"

    mkdir -p "$outdir"
    output="$outdir/${sid}_${session}_task-rest_space-MNI152NLin6Asym_seg-456CerebParcels17Networks_timeseries_desc-nocoverage.tsv"

    individual_atlas="/Users/ga0034de/Desktop/individual_atlases/${sid}_${session}_400atlases/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz"

    python /Users/ga0034de/Desktop/code_tosort/extract_uncoverage.py "$individual_atlas" "$LABELS" "$denoised_bold" "$output"

    echo "Time series extracted for $sid $session"
    echo "Output: $output"
    echo "##########################"
}

export -f extract_ts
export LABELS outdir
# sort dirs by subject and session
find "$proj_dir" -maxdepth 1 -type f -name '*.nii.gz' | sort | parallel -j 10 extract_ts {}