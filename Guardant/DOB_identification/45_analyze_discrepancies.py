import pandas as pd
import numpy as np

# Read the updated histology map
df = pd.read_csv('updated_histology_map.csv')

# Function to check if description matches subtype
def check_description_match(row):
    if pd.isna(row['Histo/Behavior ICD-O-3-Desc']):
        return 'No Description'
    
    desc = row['Histo/Behavior ICD-O-3-Desc'].lower()
    subtype = row['Histological Subtype']
    
    if subtype == 'ILC':
        if 'lobular' in desc:
            return 'Match'
        else:
            return 'Mismatch'
    elif subtype == 'IDC':
        if any(term in desc for term in ['duct', 'ductal', 'no special type']):
            return 'Match'
        else:
            return 'Mismatch'
    return 'Unknown'

# Add match status
df['Match_Status'] = df.apply(check_description_match, axis=1)

# Print overall statistics
print("\nOverall Statistics:")
print("-" * 50)
match_counts = df['Match_Status'].value_counts()
for status, count in match_counts.items():
    percentage = (count / len(df)) * 100
    print(f"{status}: {count} ({percentage:.1f}%)")

# Analyze mismatches
print("\nDetailed Mismatch Analysis:")
print("-" * 50)
mismatches = df[df['Match_Status'] == 'Mismatch']
print("\nSample of Mismatches:")
print(mismatches[['Effective Patient ID', 'Histological Subtype', 'Histo/Behavior ICD-O-3-Desc']].head())

# Analyze description patterns
print("\nUnique ICD-O-3 Descriptions:")
print("-" * 50)
desc_counts = df[df['Histo/Behavior ICD-O-3-Desc'].notna()]['Histo/Behavior ICD-O-3-Desc'].value_counts()
print(desc_counts.head())

# Save mismatches to file for further analysis
mismatches.to_csv('histology_mismatches.csv', index=False) 