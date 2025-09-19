import pandas as pd

# Read the CSV file with latin-1 encoding
file_path = 'Guardant Project - Histology Data(IDC + ILC ONLY).csv'
df = pd.read_csv(file_path, encoding='latin-1')

# Skip the first two rows (empty row and header row) and reset index
df = df.iloc[2:].reset_index(drop=True)

# Define column names
gh_id_col = 'Can delete duplicates of GH_ID'
histo_col = 'Unnamed: 15'  # Column containing histology information

# Count total unique GH_IDs
total_unique_ids = df[gh_id_col].nunique()
print(f"\nTotal number of unique GH_IDs: {total_unique_ids}")

# Count total number of duplicated GH_IDs
duplicated_ids = df[df[gh_id_col].duplicated(keep=False)]
unique_duplicated_ids = duplicated_ids[gh_id_col].unique()

print(f"Total number of records with duplicated GH_IDs: {len(duplicated_ids)}")
print(f"Number of unique GH_IDs that have duplicates: {len(unique_duplicated_ids)}")

# Analyze each duplicated GH_ID and store different histologies
same_desc = 0
diff_desc = 0
different_histologies = []

print("\nAnalyzing duplicates...")
for gh_id in unique_duplicated_ids:
    # Get all records for this GH_ID
    id_records = df[df[gh_id_col] == gh_id]
    # Check if all Histo/Behavior descriptions are the same
    unique_desc = id_records[histo_col].nunique()
    if unique_desc == 1:
        same_desc += 1
    else:
        diff_desc += 1
        # Store the different histologies
        histologies = id_records[histo_col].unique()
        different_histologies.append({
            'GH_ID': gh_id,
            'Number_of_Records': len(id_records),
            'Histology_1': histologies[0],
            'Histology_2': histologies[1] if len(histologies) > 1 else None,
            'Histology_3': histologies[2] if len(histologies) > 2 else None
        })

print(f"\nAmong duplicated GH_IDs:")
print(f"Number of IDs with same histology description: {same_desc}")
print(f"Number of IDs with different histology description: {diff_desc}")

# Save the different histologies to a CSV file
if different_histologies:
    diff_df = pd.DataFrame(different_histologies)
    output_file = 'different_histologies.csv'
    diff_df.to_csv(output_file, index=False)
    print(f"\nSaved {diff_desc} records with different histologies to {output_file}")

# Display some examples of differences
print("\nExample of records with different histology descriptions:")
for gh_id in unique_duplicated_ids:
    id_records = df[df[gh_id_col] == gh_id]
    if id_records[histo_col].nunique() > 1:
        print(f"\nGH_ID: {gh_id}")
        print("Records:")
        for idx, row in id_records.iterrows():
            print(f"  Histology: {row[histo_col]}")
        break  # Just show one example 