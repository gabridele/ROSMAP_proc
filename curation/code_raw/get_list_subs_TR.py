import pandas as pd

# Load the TSV file, output of bids validation
df = pd.read_csv("code/CuBIDS/v0_validation.tsv", sep="\t")
print(df.head())

# Filter by 'code' column for the value 'REPETITION_TIME_MISMATCH'
filtered_df = df[df['code'] == 'REPETITION_TIME_MISMATCH']

# Select only the 'location' column
location_series = filtered_df['location']

# Save to a text file
location_series.to_csv("code/TR_mismatch_raw.txt", index=False, header=False)