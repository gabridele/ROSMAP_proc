import os
import re
import datalad.api as dl


# Set the path to your main folder
main_folder_path = '/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_2'

# Iterate over each sub-#### folder
for sub_folder in os.listdir(main_folder_path):
    if not sub_folder.startswith('sub-'):
        continue

    sub_folder_path = os.path.join(main_folder_path, sub_folder)
    if not os.path.isdir(sub_folder_path):
        continue

    # Extract digits from sub-####
    sub_digits = re.search(r'sub-(\d+)', sub_folder).group(1)

    # Iterate over each subfolder inside sub-#### (e.g., ###_##_###)
    for inner_folder in os.listdir(sub_folder_path):
        inner_folder_path = os.path.join(sub_folder_path, inner_folder)
        if not os.path.isdir(inner_folder_path):
            continue
        # Iterate over files in the inner folder
        for filename in os.listdir(inner_folder_path):
            file_path = os.path.join(inner_folder_path, filename)
            if os.path.isdir(file_path):
                continue

            # Only process files containing 'Obl_2' in the name
            if 'Obl_2' in filename:

                # Extract file extension
                file_name, file_ext = os.path.splitext(filename)
                if file_name.endswith('.nii'):
                    file_name = file_name[:-4]
                if file_ext == '.gz' and '.nii' not in file_ext:
                    file_ext = '.nii' + file_ext
                dl.unlock(dataset=sub_folder_path, path=file_path)
                # Construct the new filename
                new_filename = f"sub-{sub_digits}_ses-#_acq-BNK20090211_{file_name}{file_ext}"
                new_file_path = os.path.join(inner_folder_path, new_filename)

                # Rename the file
                os.rename(file_path, new_file_path)
                print(f"Renamed {filename} to {new_filename}")
    #outside of loop so that we get only one git log per subject dataset
    dl.save(dataset=sub_folder_path, message=f"Renamed Obl_2-echo_... files in {sub_folder}")

print("Renaming process completed.")
