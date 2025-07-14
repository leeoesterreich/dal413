import pandas as pd

# Define trigger words
idc_triggers = ['no special type', 'nst', 'duct', 'ductal', 'ductular']
ilc_triggers = ['lobular']

# Function to check if any trigger word is in the text
def contains_trigger(text, triggers):
    if pd.isna(text):
        return False
    text = str(text).lower()
    return any(trigger in text for trigger in triggers)

# Read the genomic data
genomic_file = 'different_histologies_definitive.csv'
genomic_df = pd.read_csv(genomic_file)

# Classify each record based on Histology_1 (first histology description)
idc_mask = genomic_df['Histology_1'].apply(lambda x: contains_trigger(x, idc_triggers))
ilc_mask = genomic_df['Histology_1'].apply(lambda x: contains_trigger(x, ilc_triggers))

# Create separate dataframes
idc_genomic = genomic_df[idc_mask].copy()
ilc_genomic = genomic_df[ilc_mask].copy()

# Check for unclassified or double-classified records
unclassified = genomic_df[~(idc_mask | ilc_mask)]
double_classified = genomic_df[idc_mask & ilc_mask]

# Save the filtered genomic data
idc_output = 'IDC_genomic.csv'
ilc_output = 'ILC_genomic.csv'
idc_genomic.to_csv(idc_output, index=False)
ilc_genomic.to_csv(ilc_output, index=False)

# Print summary
print("\nGenomic Data Classification Summary:")
print(f"Total records in genomic data: {len(genomic_df)}")
print(f"IDC genomic records: {len(idc_genomic)}")
print(f"ILC genomic records: {len(ilc_genomic)}")
print(f"Unclassified records: {len(unclassified)}")
print(f"Double-classified records: {len(double_classified)}")

if len(unclassified) > 0:
    print("\nUnclassified genomic records:")
    print(unclassified[['GH_ID', 'Histology_1']].to_string())

if len(double_classified) > 0:
    print("\nDouble-classified genomic records:")
    print(double_classified[['GH_ID', 'Histology_1']].to_string())

print(f"\nSaved IDC genomic data to {idc_output}")
print(f"Saved ILC genomic data to {ilc_output}")

# Print samples of the classified data
print("\nSample of IDC genomic records:")
if len(idc_genomic) > 0:
    print(idc_genomic[['GH_ID', 'Histology_1', 'Histology_2', 'Histology_3']].head().to_string())

print("\nSample of ILC genomic records:")
if len(ilc_genomic) > 0:
    print(ilc_genomic[['GH_ID', 'Histology_1', 'Histology_2', 'Histology_3']].head().to_string()) 