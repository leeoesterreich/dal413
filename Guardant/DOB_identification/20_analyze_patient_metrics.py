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

def analyze_molecular_markers(df, cohort_name):
    print(f"\n{cohort_name} Cohort Molecular Markers Analysis (Filtered for Vus=FALSE):")
    print("-" * (len(cohort_name) + 45))
    
    # Count unique patients
    unique_patients = df['Effective Patient ID'].nunique()
    print(f"\nTotal unique patients: {unique_patients}")
    
    # Analyze molecular response score
    print("\nMolecular Response Score Analysis:")
    response_scores = df.groupby('Effective Patient ID')['Molecular Response Score'].first().value_counts()
    for score, count in response_scores.items():
        if pd.notna(score):
            print(f"- {score}: {count} patients ({count/unique_patients*100:.1f}%)")
    
    # Analyze ctDNA level changes
    print("\nMolecular Response ctDNA Level Changes:")
    ctdna_changes = df.groupby('Effective Patient ID')['Molecular Response ctDNA Level Change'].first().value_counts()
    for change, count in ctdna_changes.items():
        if pd.notna(change):
            print(f"- {change}: {count} patients ({count/unique_patients*100:.1f}%)")
    
    # Analyze suspected germline variants
    print("\nSuspected Germline Variants:")
    germline_vars = df.groupby('Effective Patient ID')['Suspected Germline Variant'].first().value_counts()
    for variant, count in germline_vars.items():
        if pd.notna(variant):
            print(f"- {variant}: {count} patients ({count/unique_patients*100:.1f}%)")
    
    # Analyze promoter methylation
    print("\nPromoter Methylation Status:")
    methylation = df.groupby('Effective Patient ID')['Promoter Methylation'].first().value_counts()
    for status, count in methylation.items():
        if pd.notna(status):
            print(f"- {status}: {count} patients ({count/unique_patients*100:.1f}%)")

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
analyze_molecular_markers(idc_genomic, "IDC")
analyze_molecular_markers(ilc_genomic, "ILC")

# Generate detailed molecular markers report
print("\nGenerating detailed molecular markers report...")

def generate_molecular_markers_details(df, cohort):
    # Get patient-level molecular marker data
    marker_data = df.groupby('Effective Patient ID').agg({
        'Molecular Response Score': 'first',
        'Molecular Response ctDNA Level Change': 'first',
        'Suspected Germline Variant': 'first',
        'Promoter Methylation': 'first'
    }).reset_index()
    
    marker_data['Cohort'] = cohort
    
    return marker_data

idc_details = generate_molecular_markers_details(idc_genomic, "IDC")
ilc_details = generate_molecular_markers_details(ilc_genomic, "ILC")

# Combine and save details
all_details = pd.concat([idc_details, ilc_details])
all_details = all_details.sort_values('Cohort')
all_details.to_csv('Vus_filtering/molecular_markers_details.csv', index=False)
print("Detailed molecular markers report saved to: Vus_filtering/molecular_markers_details.csv")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('Vus_filtering/molecular_markers_analysis_report.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: Vus_filtering/molecular_markers_analysis_report.txt") 