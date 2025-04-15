#!/bin/bash

SOURCE_DIR="/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2/sub-99911705/100106_00_99911705"
DEST_DIR="."

for file in "$SOURCE_DIR"/*3D*; do
    rsync -av --copy-links "$file" "$DEST_DIR"
done