import pandas as pd
import re

# Read the filtered cohort
input_file = 'filtered_cohort.csv'
df = pd.read_csv(input_file)

# Define trigger words
idc_triggers = ['no special type', 'nst', 'duct', 'ductal', 'ductular']
ilc_triggers = ['lobular']

# Define the histology column
histo_col = 'Unnamed: 15'  # Column containing histology information
gh_id_col = 'Can delete duplicates of GH_ID'

# Function to check if any trigger word is in the text
def contains_trigger(text, triggers):
    if pd.isna(text):
        return False
    text = str(text).lower()
    return any(trigger in text for trigger in triggers)

# Classify each record
idc_mask = df[histo_col].apply(lambda x: contains_trigger(x, idc_triggers))
ilc_mask = df[histo_col].apply(lambda x: contains_trigger(x, ilc_triggers))

# Create separate dataframes
idc_df = df[idc_mask].copy()
ilc_df = df[ilc_mask].copy()

# Check for unclassified or double-classified records
total_classified = len(idc_df) + len(ilc_df)
unclassified = df[~(idc_mask | ilc_mask)]
double_classified = df[idc_mask & ilc_mask]

# Print summary
print("\nClassification Summary:")
print(f"Total records in filtered cohort: {len(df)}")
print(f"IDC records: {len(idc_df)}")
print(f"ILC records: {len(ilc_df)}")
print(f"Unclassified records: {len(unclassified)}")
print(f"Double-classified records: {len(double_classified)}")

if len(unclassified) > 0:
    print("\nUnclassified records histology descriptions:")
    for idx, row in unclassified.iterrows():
        print(f"GH_ID: {row[gh_id_col]}, Histology: {row[histo_col]}")

if len(double_classified) > 0:
    print("\nDouble-classified records histology descriptions:")
    for idx, row in double_classified.iterrows():
        print(f"GH_ID: {row[gh_id_col]}, Histology: {row[histo_col]}")

# Save the classified datasets
idc_output = 'idc_cohort.csv'
ilc_output = 'ilc_cohort.csv'
idc_df.to_csv(idc_output, index=False)
ilc_df.to_csv(ilc_output, index=False)

print(f"\nSaved IDC cohort to {idc_output}")
print(f"Saved ILC cohort to {ilc_output}") 