#!/bin/bash

for dir in sub-*; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        cd "$dir"
        for subdir in *; do
            cd "$subdir"
            if [ -d "dwi" ]; then 
                rm -rf "dwi"
            else
                echo "no folder found"
            fi
            cd ..
        cd ..
        done
    fi
done