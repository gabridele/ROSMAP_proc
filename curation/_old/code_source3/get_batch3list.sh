#!/bin/bash

# Output CSV file
output_file="code/batch3_subs.csv"

# Write header
echo "sub-ID,visit,date,protocol" > "$output_file.tmp"

# Iterate through subject directories
for sub_dir in sub-*; do
    sub_id=$(basename "$sub_dir")  # Extract sub-ID from parent folder name (e.g., sub-12345678)
    
    # Iterate through child directories inside the sub-# folder
    for child_dir in "$sub_dir"/*/; do
        [ -d "$child_dir" ] || continue  # Skip if not a directory
        
        # Extract the date (first part before the underscore in the child folder name)
        date=$(basename "$child_dir" | awk -F'_' '{print $1}')
        
        # Extract the visit (second part between underscores in the child folder name)
        visit=$(basename "$child_dir" | awk -F'_' '{print $2}')
        
        # Protocol is fixed (shouldnt be hardcoded, but I fixed it later on)
        protocol="RIRC230726"
        
        # Write to the output CSV
        echo "$sub_id,$visit,$date,$protocol" >> "$output_file.tmp"
    done
done

# Optional: Sort the output file by sub-ID, visit, and date if needed
sort -t',' -k1,1 -k2,2 -k3,3 "$output_file.tmp" > "$output_file"
rm "$output_file.tmp"