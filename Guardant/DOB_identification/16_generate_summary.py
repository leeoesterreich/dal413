import pandas as pd
import numpy as np
import sys
from io import StringIO
import datetime
import os

# Create output directory if it doesn't exist
os.makedirs('Vus_filtering', exist_ok=True)

# Create a string buffer to capture all output
output_buffer = StringIO()
original_stdout = sys.stdout
sys.stdout = output_buffer

def analyze_alterations(df, cohort_name):
    print(f"\n{cohort_name} Cohort Analysis (Filtered for Vus=FALSE):")
    print("-" * (len(cohort_name) + 35))
    
    # Count unique alterations
    unique_alterations = df['Alteration'].dropna().unique()
    print(f"\nTotal unique alterations: {len(unique_alterations)}")
    
    # Count patients with alterations
    patients_with_alterations = df['Effective Patient ID'].nunique()
    print(f"Total patients with alterations: {patients_with_alterations}")
    
    # Most common alterations
    print("\nTop 10 most common alterations:")
    alteration_counts = df['Alteration'].value_counts().head(10)
    for alt, count in alteration_counts.items():
        print(f"- {alt}: {count} occurrences")
    
    # Analyze alteration types
    print("\nAlteration types distribution:")
    type_counts = df['Type'].value_counts()
    for type_name, count in type_counts.items():
        print(f"- {type_name}: {count}")
    
    # Analyze mutation types
    print("\nMutation types distribution:")
    mutation_type_counts = df['Mutation Type'].value_counts()
    for mut_type, count in mutation_type_counts.items():
        print(f"- {mut_type}: {count}")

# Add timestamp
print(f"Analysis performed on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Read the data files
print("\nReading data files...")
try:
    idc_genomic = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/IDC_genomic.csv', 
                             encoding='latin-1', 
                             low_memory=False)
    ilc_genomic = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/ILC_genomic.csv', 
                             encoding='latin-1', 
                             low_memory=False)
except FileNotFoundError as e:
    print(f"Error: An input file was not found. {e}")
    sys.exit(1)

# Print column names for debugging
print("\nAvailable columns in IDC genomic data:")
for col in idc_genomic.columns:
    print(f"'{col}'")

print("\nAvailable columns in ILC genomic data:")
for col in ilc_genomic.columns:
    print(f"'{col}'")

# Filter for Vus = FALSE
idc_genomic = idc_genomic[idc_genomic['Vus'] == "FALSE"]
ilc_genomic = ilc_genomic[ilc_genomic['Vus'] == "FALSE"]

# Analyze both cohorts
analyze_alterations(idc_genomic, "IDC")
analyze_alterations(ilc_genomic, "ILC")

# Generate detailed alteration report
print("\nGenerating detailed alteration report...")

def generate_alteration_details(df, cohort):
    # Get alteration counts
    alt_counts = df['Alteration'].value_counts().reset_index()
    alt_counts.columns = ['Alteration', 'Total_Occurrences']
    alt_counts['Cohort'] = cohort
    alt_counts['Unique_Patients'] = df.groupby('Alteration')['Effective Patient ID'].nunique().reindex(alt_counts['Alteration']).values
    
    # Add type information
    alt_counts['Type'] = df.groupby('Alteration')['Type'].first().reindex(alt_counts['Alteration']).values
    alt_counts['Mutation_Type'] = df.groupby('Alteration')['Mutation Type'].first().reindex(alt_counts['Alteration']).values
    
    return alt_counts

idc_details = generate_alteration_details(idc_genomic, "IDC")
ilc_details = generate_alteration_details(ilc_genomic, "ILC")

# Combine and save details
all_details = pd.concat([idc_details, ilc_details])
all_details = all_details.sort_values(['Cohort', 'Total_Occurrences'], ascending=[True, False])
all_details.to_csv('Vus_filtering/alteration_details.csv', index=False)
print("Detailed alteration report saved to: Vus_filtering/alteration_details.csv")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('Vus_filtering/alteration_analysis_report.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: Vus_filtering/alteration_analysis_report.txt") 