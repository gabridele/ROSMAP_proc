import pandas as pd
import glob
import os

# pattern to match all your files
pattern = "/Users/ga0034de/Desktop/BNK_BBRnoBBR_plot/covg_bbr/*seg-4S456Parcels_stat-coverage_bold.tsv"

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
results_df.to_csv("bbr_mean_coverage_456parcels_priorityl.tsv", sep="\t", index=False)