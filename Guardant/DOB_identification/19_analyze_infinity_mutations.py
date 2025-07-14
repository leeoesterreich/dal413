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

def analyze_therapies(df, cohort_name):
    print(f"\n{cohort_name} Cohort Therapy Analysis (Filtered for Vus=FALSE):")
    print("-" * (len(cohort_name) + 35))
    
    # Count unique patients
    unique_patients = df['Effective Patient ID'].nunique()
    print(f"\nTotal unique patients: {unique_patients}")
    
    # Analyze approved therapies in disease
    print("\nApproved Therapies in Disease:")
    therapy_counts = df.groupby('Effective Patient ID')['Therapies Approved in Disease'].first().value_counts()
    for therapy, count in therapy_counts.items():
        if pd.notna(therapy):
            print(f"- {therapy}: {count} patients ({count/unique_patients*100:.1f}%)")
    
    # Analyze approved therapies in other diseases
    print("\nApproved Therapies in Other Diseases:")
    other_therapy_counts = df.groupby('Effective Patient ID')['Approved in Other Diseases'].first().value_counts()
    for therapy, count in other_therapy_counts.items():
        if pd.notna(therapy):
            print(f"- {therapy}: {count} patients ({count/unique_patients*100:.1f}%)")
    
    # Analyze resistance therapies
    print("\nTherapies Associated with Resistance:")
    resistance_counts = df.groupby('Effective Patient ID')['Therapies Associated Resistance'].first().value_counts()
    for therapy, count in resistance_counts.items():
        if pd.notna(therapy):
            print(f"- {therapy}: {count} patients ({count/unique_patients*100:.1f}%)")

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
analyze_therapies(idc_genomic, "IDC")
analyze_therapies(ilc_genomic, "ILC")

# Generate detailed therapy report
print("\nGenerating detailed therapy report...")

def generate_therapy_details(df, cohort):
    # Get patient-level therapy data
    therapy_data = df.groupby('Effective Patient ID').agg({
        'Therapies Approved in Disease': 'first',
        'Approved in Other Diseases': 'first',
        'Therapies Associated Resistance': 'first'
    }).reset_index()
    
    therapy_data['Cohort'] = cohort
    
    return therapy_data

idc_details = generate_therapy_details(idc_genomic, "IDC")
ilc_details = generate_therapy_details(ilc_genomic, "ILC")

# Combine and save details
all_details = pd.concat([idc_details, ilc_details])
all_details = all_details.sort_values('Cohort')
all_details.to_csv('Vus_filtering/therapy_details.csv', index=False)
print("Detailed therapy report saved to: Vus_filtering/therapy_details.csv")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('Vus_filtering/therapy_analysis_report.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: Vus_filtering/therapy_analysis_report.txt") 