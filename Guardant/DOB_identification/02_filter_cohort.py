import pandas as pd

# Read the original data file
original_file = 'Guardant Project - Histology Data(IDC + ILC ONLY).csv'

# First read the third row to get proper column names
header_df = pd.read_csv(original_file, encoding='latin-1', nrows=3)
column_names = header_df.iloc[2]  # Get the third row (index 2) as column names

# Now read the data using these column names, skipping the first 3 rows
df = pd.read_csv(original_file, encoding='latin-1', skiprows=3, names=column_names)

# Read the list of patients to exclude
exclude_file = 'different_histologies_definitive.csv'
exclude_df = pd.read_csv(exclude_file)
gh_ids_to_exclude = exclude_df['GH_ID'].tolist()

# Remove the patients to exclude
df_filtered = df[~df['GH_ID'].isin(gh_ids_to_exclude)]

# Keep only the first occurrence for duplicates
df_final = df_filtered.drop_duplicates(subset=['GH_ID'], keep='first')

# Create histology mapping based on the description
df_final['Histology'] = df_final['Histo/Behavior ICD-O-3-Desc'].apply(
    lambda x: 'IDC' if 'ductal' in str(x).lower() else ('ILC' if 'lobular' in str(x).lower() else 'Unknown')
)

# Select relevant columns for output
columns_to_keep = [
    'GH_ID', 'Date of Birth', 'Age at DX', 'Sex', 'Primary Site-Desc',
    'Laterality-Desc', 'Histo/Behavior ICD-O-3-Desc', 'Grade/Differentiation-Desc',
    'Regional Nodes Positive', 'Regional Nodes Exam', 'Histology'
]

df_final = df_final[columns_to_keep]

# Split into IDC and ILC cohorts
idc_df = df_final[df_final['Histology'] == 'IDC']
ilc_df = df_final[df_final['Histology'] == 'ILC']

# Save the filtered datasets
filtered_output = 'filtered_cohort.csv'
idc_output = 'idc_cohort.csv'
ilc_output = 'ilc_cohort.csv'

df_final.to_csv(filtered_output, index=False)
idc_df.to_csv(idc_output, index=False)
ilc_df.to_csv(ilc_output, index=False)

# Print summary
print(f"\nOriginal number of unique patients: {df['GH_ID'].nunique()}")
print(f"Number of patients excluded: {len(gh_ids_to_exclude)}")
print(f"Final number of unique patients: {df_final['GH_ID'].nunique()}")
print(f"\nSaved filtered cohort to {filtered_output}")
print(f"IDC patients: {len(idc_df)}")
print(f"ILC patients: {len(ilc_df)}")
print(f"\nSaved IDC cohort to {idc_output}")
print(f"Saved ILC cohort to {ilc_output}") 