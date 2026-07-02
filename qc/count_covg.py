import pandas as pd
import glob
import os

# new pattern
pattern = "/Users/ga0034de/Desktop/freesurfernobbr_xcpd_output/sub-*xcp_d*/xcp_d_nifti/sub-*/ses-*/func/*seg-4S156Parcels_stat-coverage_bold.tsv"

# old 
#pattern="/Users/ga0034de/Desktop/old_proc_priorityIDs/sub-*xcpd*/xcp_d_nifti/sub-*/ses-*/func/*seg-4S156Parcels_stat-coverage_bold.tsv"
files = glob.glob(pattern)

print(f"Found {len(files)} files matching the pattern.")
results = []

for f in files:
    print(f"Processing file: {f}")
    df = pd.read_csv(f, sep="\t")
    
    
    # second column (index 1)
    col = df.iloc[:, 1]
    
    zero_count = (col < 0.5).sum()  # count where value < 0.5
    sub_session=os.path.basename(f).split("_task-rest")[0]  # extract sub- and ses- info
    results.append((sub_session, zero_count))

# print results
for fname, val in results:
    print(f"{fname}: {val:.2f}%")
    # save results to new tsv file

results_df = pd.DataFrame(results, columns=["sub_ses", "zero_count"])
results_df.to_csv("freesurfnobbr_count_non-coverage_156parcels_priorityl.tsv", sep="\t", index=False)