#!/bin/bash
set -euo pipefail

#module load fsl/6.0

sub_ses=$1
# get sub_ses from filename up to second _
sub_ses=$(echo "$sub_ses" | cut -d'_' -f1,2)

bids_dir="/Users/ga0034de/Desktop/denoised_bnks/v1.3/bnks/probseg"
atlases_path="/Volumes/research/LU26D1023-DemonLab/DemonLab/ROSMAP/derivatives/atlases"

atlas400="${atlases_path}/atlas-Schaefer400TianS2Cereb/MNINLin6Asym/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz"


fmriprep_dir="${bids_dir}/${sub_ses}_generated"

mask_gm="/Users/ga0034de/Desktop/individual_atlases/mask_gm/${sub_ses}_space-MNI152NLin6Asym_res-2_label-GM_mask-03.nii.gz"
# skip if existing
#if [[ -f "$mask_gm" ]]; then
#  echo "GM mask already exists: $mask_gm"
#  exit 0
#fi

probseg_gm=$(find "." \
  -type f \
  -path "*_space-MNI152NLin6Asym_res-2_label-GM_probseg.nii.gz" \
  -print -quit)


# ------ create a GM mask

fslmaths "$probseg_gm" -thr 0.3 -bin "$mask_gm"
echo $mask_gm

# ------ synchronize necessary files for xcp_d
sid_atlas_dir="/Users/ga0034de/Desktop/individual_atlases/${sub_ses}_400atlases"
mkdir -p "$sid_atlas_dir"
#dir200="${sid_atlas_dir}/schaefer_supplemented/atlas-Schaefer200TianS2Cereb"
dir400="${sid_atlas_dir}"

# mkdir -p "$sid_atlas_dir"
# mkdir -p "$dir200"
# mkdir -p "$dir400"

# rsync -aP "${atlases_path}/dataset_description.json" "${sid_atlas_dir}/dataset_description.json" 

# rsync -aP "${atlases_path}/schaefer_supplemented/atlas-Schaefer200TianS2Cereb/atlas-Schaefer200TianS2Cereb_dseg.tsv" "${dir200}/atlas-Schaefer200TianS2Cereb_dseg.tsv" 
# rsync -aP "${atlases_path}/schaefer_supplemented/atlas-Schaefer200TianS2Cereb/atlas-Schaefer200TianS2Cereb_space-MNI152NLin2006cAsym_res-2_dseg.json" "${dir200}/atlas-Schaefer200TianS2Cereb_space-MNI152NLin2006cAsym_res-2_dseg.json" 

# rsync -aP "${atlases_path}/schaefer_supplemented/atlas-Schaefer400TianS2Cereb/atlas-Schaefer400TianS2Cereb_dseg.tsv" "${dir400}/atlas-Schaefer400TianS2Cereb_dseg.tsv" 
# rsync -aP "${atlases_path}/schaefer_supplemented/atlas-Schaefer400TianS2Cereb/atlas-Schaefer400TianS2Cereb_space-MNI152NLin2006cAsym_res-2_dseg.json" "${dir400}/atlas-Schaefer400TianS2Cereb_space-MNI152NLin2006cAsym_res-2_dseg.json" 

echo $dir400 
# ------ create an individualized parcellation
#fslmaths "$atlas200" -mul "$mask_gm" "${dir200}/atlas-Schaefer200TianS2Cereb_space-MNI152NLin2006cAsym_res-2_dseg.nii.gz" 
#skip if existing
#if [[ -f "${dir400}/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz" ]]; then
#  echo "Individualized atlas already exists: ${dir400}/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz"
#  exit 0
#fi
fslmaths "$atlas400" -mul "$mask_gm" "${dir400}/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz" 

#echo "${dir200}/atlas-Schaefer200TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz was created"
echo "${dir400}/atlas-Schaefer400TianS2Cereb_space-MNI152NLin6Asym_res-2_dseg.nii.gz was created"
