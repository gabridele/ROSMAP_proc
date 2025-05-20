#!/bin/bash

# Output CSV file
output_file="code/batch1_subs.csv"

# Write header
echo "sub-ID,ses-number,protocol" > "$output_file.tmp"

# Iterate through subject directories
for sub_dir in sub-*; do
    sub_id=$(basename "$sub_dir")  # Extract sub-ID

    # Iterate through session directories
    for ses_dir in "$sub_dir"/ses-*; do
        [ -d "$ses_dir" ] || continue  # Skip if not a directory
        ses_number=$(basename "$ses_dir")  # Extract session number

        # Look for JSON files in func directory
        for json_file in "$ses_dir"/func/*.json; do
            [ -f "$json_file" ] || continue  # Skip if not a file

            # Extract protocol
            if [[ "$json_file" =~ acq-([A-Za-z]+[0-9]{8}) ]]; then
                protocol="${BASH_REMATCH[1]}"
                echo "$sub_id,$ses_number,$protocol" >> "$output_file.tmp"
            fi
        done
    done
done
# Sort the CSV file by subject, session, and protocol
(head -n 1 "$output_file.tmp" && tail -n +2 "$output_file.tmp" | sort -t',' -k1,1 -k2,2 -k3,3) > "$output_file"

rm "$output_file.tmp"
echo "CSV file generated: $output_file"