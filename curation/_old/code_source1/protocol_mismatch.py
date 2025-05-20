import pandas as pd
import re

data = pd.read_csv("code/qc.csv", delimiter=",")

df = pd.DataFrame(data)

# Function to extract the date from 'acq'
def extract_date(acq_entry):
    match = re.search(r'^[A-Za-z]+(\d{8})', acq_entry)
    return str(match.group(1)) if match else None

# Apply function to extract date
df["extracted_date"] = df["acq"].apply(extract_date)

# Now filter columns to retain only the desired ones, including 'extracted_date' temporarily
df_filtered = df[["sub", "ses", "task", "acq", "ScanDate", "Visit", "projid", "Protocol", "ScannerGroup", "extracted_date"]]

# Compare extracted date with Protocol and create a match column
df_filtered["match"] = df_filtered["extracted_date"] == df_filtered["Protocol"]

# Drop 'extracted_date' if it's not needed in the final output
df_filtered.drop(columns=["extracted_date"], inplace=True)

# Display result
# print(df_filtered) 

# Filter only rows where match is False
df_false = df_filtered[df_filtered["match"] == False]

# Drop 'extracted_date' if it's not needed in the final output
#df_false.drop(columns=["extracted_date", "match"], inplace=True)

# Display only the rows where the dates do not match
df_false.to_csv("code/filtered_false_entries.csv", index=False)