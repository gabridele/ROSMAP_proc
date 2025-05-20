#!/bin/bash

TARGET_DIR="/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2/sub-63698192/110119_04_63698192"

cp "/home/gabridele/backup/ROSMAP/batch_2/BNK/090211/110119_04_63698192/3D_MPRAGEa.json" "$TARGET_DIR"
cp "/home/gabridele/backup/ROSMAP/batch_2/BNK/090211/110119_04_63698192/3D_MPRAGE.json" "$TARGET_DIR"
cp "/home/gabridele/backup/ROSMAP/batch_2/BNK/090211/110119_04_63698192/3D_MPRAGE.nii.gz" "$TARGET_DIR"
cp "/home/gabridele/backup/ROSMAP/batch_2/BNK/090211/110119_04_63698192/Obl_2-echo_GRE_(phase_map)_e1.json" "$TARGET_DIR"

cd "$TARGET_DIR"
mv "3D_MPRAGEa.json" "sub-63698192_ses-#_acq-BNK20090211_3D_MPRAGEa_T1w.json"
mv "3D_MPRAGE.json" "sub-63698192_ses-#_acq-BNK20090211_3D_MPRAGE_T1w.json"
mv "3D_MPRAGE.nii.gz" "sub-63698192_ses-#_acq-BNK20090211_3D_MPRAGE_T1w.nii.gz"
mv "Obl_2-echo_GRE_(phase_map)_e1.json" "sub-63698192_ses-#_acq-BNK20090211_Obl_2-echo_GRE_(phase_map)_e1.json"

datalad save -m "Copied missing files from $SOURCE_DIR to $TARGET_DIR" "$TARGET_DIR"
############################################################################################################
# Define another source and target directories

ANOTHER_TARGET_DIR="/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2/sub-81569322/120315_02_81569322"

cp "/home/gabridele/backup/ROSMAP/batch_2/BNK/090211/120315_02_81569322/3D_MPRAGE.nii.gz" "$ANOTHER_TARGET_DIR"
cp "/home/gabridele/backup/ROSMAP/batch_2/BNK/090211/120315_02_81569322/P18432.zip" "$ANOTHER_TARGET_DIR"

cd "$ANOTHER_TARGET_DIR"
mv "3D_MPRAGE.nii.gz" "/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2/sub-81569322/120315_02_81569322/sub-81569322_ses-#_acq-BNK20090211_3D_MPRAGE_T1w.nii.gz"

datalad save -m "Copied missing files from $ANOTHER_SOURCE_DIR to $ANOTHER_TARGET_DIR" "$ANOTHER_TARGET_DIR"