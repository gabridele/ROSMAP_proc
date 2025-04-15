#!/bin/bash
echo "Removing brackets from filenames..."
pwd

for dir in sub-*; do
    if [[ -d "$dir" ]]; then
        cd "$dir" || exit  # Ensure we don't proceed if cd fails
        for subdir in */; do
            if [[ -d "$subdir" ]]; then
                cd "$subdir" || exit
                for file in *phase_map*; do
                    if [[ -f "$file" ]]; then
                        echo "Removing brackets from $file..."
                        # Remove round and square brackets from the filename
                        new_file=$(echo "$file" | sed 's/[()]//g')
                        if [[ "$file" != "$new_file" ]]; then  # Avoid unnecessary renaming
                            mv "$file" "$new_file"
                        fi
                    fi
                done
                
                cd .. || exit
            fi
        done
        datalad save -m "Removed brackets from phasemap filenames" .
        cd .. || exit
    fi
done
echo "Done."