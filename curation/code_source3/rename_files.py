import os
import re
import datalad.api as dl


# Set the path to your main folder
main_folder_path = '/home/gabridele/backup/ROSMAP_proc/sourcedata/batch_3'

# Iterate over each sub-#### folder
for sub_folder in os.listdir(main_folder_path):
    if not sub_folder.startswith('sub-'):
        continue

    sub_folder_path = os.path.join(main_folder_path, sub_folder)
    if not os.path.isdir(sub_folder_path):
        continue

    # Extract digits from sub-####
    sub_digits = re.search(r'sub-(\d+)', sub_folder).group(1)

    # Iterate over each subfolder inside sub-#### (e.g., ###_##_###) sub_folder_path = path to sub-### folder
    for inner_folder in os.listdir(sub_folder_path):
        inner_folder_path = os.path.join(sub_folder_path, inner_folder)
        if inner_folder.startswith("."):
            continue
        dwi_folder = os.path.join(inner_folder_path, 'dwi')
        os.makedirs(dwi_folder, exist_ok=True)
        if not os.path.isdir(inner_folder_path):
            continue
        # Iterate over files in the inner folder
        for filename in os.listdir(inner_folder_path):
            file_path = os.path.join(inner_folder_path, filename)
            if os.path.isdir(file_path):
                continue

            # Extract file extension
            file_name, file_ext = os.path.splitext(filename)
            if file_name.endswith('.nii'):
                file_name = file_name[:-4]
            if file_ext == '.gz' and '.nii' not in file_ext:
                file_ext = '.nii' + file_ext
            dl.unlock(dataset=sub_folder_path, path=file_path)

            if sub_digits == "13044513" or "70876731": ###!!! make sure tyhese subs dont have acqs in the other protocol too
                protocol = "230726"
            else:
                protocol = "230803"

            # Only process files containing 'x' in the name
            if 'DWI' in filename:

                # Construct the new filename
                if 'ADC' in filename:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_dir-PA_Ax-40-6-ADC_dwi{file_ext}"
                else:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_dir-PA_Ax-40-6_dwi{file_ext}"
                
            if 'FLAIR' in filename: 
                if 'ND' in filename:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_Sag-3D-ND_FLAIR{file_ext}"
                else:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_Sag-3D_FLAIR{file_ext}"
            if 'T1' in filename:
                if 'RMS' in filename:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_Sag-3D_MEMPRAGE_RMS_T1w{file_ext}"
                else:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_Sag-3D_MEMPRAGE_T1w{file_ext}"
            if 'T2' in filename:
                if 'ND' in filename:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_Sag-3D-SPACE-ND_T2w{file_ext}"
                else:
                    new_filename = f"sub-{sub_digits}_ses-#_acq-RIRC{protocol}_Sag-3D-SPACE_T2w{file_ext}"

            new_file_path = os.path.join(inner_folder_path, new_filename)
            os.rename(file_path, new_file_path)
            print(f"Renamed {filename} to {new_filename}")


        #outside of loop so that we get only one git log per subject dataset
        dl.save(dataset=sub_folder_path, message=f"Renamed files in {sub_folder} according to BIDS format")
        # dl.save(path='*FLAIR*', dataset=sub_folder_path, message=f"Renamed FLAIR files in {sub_folder}")
        # dl.save(path='*T1w*', dataset=sub_folder_path, message=f"Renamed T1w files in {sub_folder}")
        # dl.save(path='*T2w*', dataset=sub_folder_path, message=f"Renamed dwi files in {sub_folder}")

        #maybe add statement that makes sure each subdset is clean and everything got saved
print("Renaming process completed.")
