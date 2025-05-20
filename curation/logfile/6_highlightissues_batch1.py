import os
import pandas as pd
import numpy as np

df_main = pd.read_csv("code/dates_scraping/5_mods_suffix.csv")
df_main = df_main.copy()

df_main['issue'] = np.nan
df_main['sub_id'] = df_main['sub_id'].astype(str).str.replace('.0', '', regex=False)
df_main['sub_id'] = df_main['sub_id'].str.lstrip('0')
df_main['session'] = df_main['session'].replace('0.0', '0')
df_main['session'] = df_main['session'].astype(str).str.replace('.0', '', regex=False)

filtered_df = df_main[df_main['batch_number'] == 'batch_1']
filtered_df = filtered_df[filtered_df['suffix'].notna()]
grouped = filtered_df.groupby(['sub_id', 'session'])

groups_with_site_only = 0
groups_with_protocol_only = 0
groups_with_both_issues = 0

# Iterate through each group and identify issues
for (sub_id, session), group in grouped:
    issues = []
    
    # Check for multiple sites within the group
    has_site_issue = group['site'].nunique() > 1
    if has_site_issue:
        issues.append('site')
    # Check for multiple protocols within the group
    has_protocol_issue = group['protocol'].nunique() > 1
    if has_protocol_issue:
        issues.append('protocol')
    
    # Update the counters based on the issues found
    if has_site_issue and has_protocol_issue:
        groups_with_both_issues += 1
    elif has_site_issue:
        groups_with_site_only += 1
    elif has_protocol_issue:
        groups_with_protocol_only += 1
    
    # Update the 'issue' column for all rows in the group
    if issues:
        issue_str = ', '.join(issues)
        df_main.loc[group.index, 'issue'] = issue_str

# Save the updated DataFrame to a CSV file
df_main.to_csv("code/dates_scraping/6_highlightissues_batch1.csv", index=False)

# Print the counts for each type of issue
print(f"Number of groups with only 'site' issues: {groups_with_site_only}")
print(f"Number of groups with only 'protocol' issues: {groups_with_protocol_only}")
print(f"Number of groups with both 'site' and 'protocol' issues: {groups_with_both_issues}")