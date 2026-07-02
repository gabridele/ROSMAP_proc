import pandas as pd
import glob
import os

# pattern to match all your files
pattern = "/Volumes/GabrieleSSD/v2525_copy_ssd/sub-*xcpd*/xcp_d_nifti/sub-*/ses-*/func/*seg-4S156Parcels_stat-coverage_bold.tsv"

files = glob.glob(pattern)
print(f"Found {len(files)} files matching the pattern.")
results = []

for f in files:
    print(f"Processing file: {f}")
    df = pd.read_csv(f, sep="\t")
    
    # second column (index 1)
    col = df.iloc[:, 1]
    
    mean_val = col.mean() * 100  # convert to percentage
    sub_session=os.path.basename(f).split("_task-rest")[0]  # extract sub- and ses- info
    results.append((sub_session, mean_val))

# print results
for fname, val in results:
    print(f"{fname}: {val:.2f}%")
    # save results to new tsv file

results_df = pd.DataFrame(results, columns=["sub_ses", "mean_coverage_percent"])
results_df.to_csv("v2525_mean_coverage_156parcels_priorityl.tsv", sep="\t", index=False)