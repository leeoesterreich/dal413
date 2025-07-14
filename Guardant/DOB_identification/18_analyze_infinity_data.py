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

def analyze_test_metrics(df, cohort_name):
    print(f"\n{cohort_name} Cohort Test Metrics Analysis (Filtered for Vus=FALSE):")
    print("-" * (len(cohort_name) + 45))
    
    # Count unique tests
    unique_tests = df['Effective Patient ID'].nunique()
    print(f"\nTotal unique tests: {unique_tests}")
    
    # Analyze tumor fraction
    print("\nTumor Fraction Analysis:")
    tumor_stats = df.groupby('Effective Patient ID')['Tumor Fraction'].first().describe()
    print(f"- Mean tumor fraction: {tumor_stats['mean']:.3f}")
    print(f"- Median tumor fraction: {tumor_stats['50%']:.3f}")
    print(f"- Range: {tumor_stats['min']:.3f} - {tumor_stats['max']:.3f}")
    
    # TMB Score analysis
    print("\nTMB Score Analysis:")
    tmb_stats = df.groupby('Effective Patient ID')['TMB Score'].first().describe()
    print(f"- Mean TMB score: {tmb_stats['mean']:.2f}")
    print(f"- Median TMB score: {tmb_stats['50%']:.2f}")
    print(f"- Range: {tmb_stats['min']:.2f} - {tmb_stats['max']:.2f}")
    
    # MSI Status
    print("\nMSI Status Distribution:")
    msi_dist = df.groupby('Effective Patient ID')['MSI High'].first().value_counts()
    for status, count in msi_dist.items():
        print(f"- {status}: {count} tests ({count/unique_tests*100:.1f}%)")
    
    # PDL1 Analysis
    print("\nPDL1 Score Analysis:")
    pdl1_stats = df.groupby('Effective Patient ID')['PDL1 Value'].first().describe()
    print(f"- Mean PDL1 value: {pdl1_stats['mean']:.2f}")
    print(f"- Median PDL1 value: {pdl1_stats['50%']:.2f}")
    print(f"- Range: {pdl1_stats['min']:.2f} - {pdl1_stats['max']:.2f}")

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
analyze_test_metrics(idc_genomic, "IDC")
analyze_test_metrics(ilc_genomic, "ILC")

# Generate detailed test metrics report
print("\nGenerating detailed test metrics report...")

def generate_test_metrics_details(df, cohort):
    # Get test-level metrics
    test_metrics = df.groupby('Effective Patient ID').agg({
        'Tumor Fraction': 'first',
        'TMB Score': 'first',
        'MSI High': 'first',
        'PDL1 Value': 'first',
        'PDL1 Score Type': 'first'
    }).reset_index()
    
    test_metrics['Cohort'] = cohort
    
    return test_metrics

idc_details = generate_test_metrics_details(idc_genomic, "IDC")
ilc_details = generate_test_metrics_details(ilc_genomic, "ILC")

# Combine and save details
all_details = pd.concat([idc_details, ilc_details])
all_details = all_details.sort_values(['Cohort', 'Tumor Fraction'], ascending=[True, False])
all_details.to_csv('Vus_filtering/test_metrics_details.csv', index=False)
print("Detailed test metrics report saved to: Vus_filtering/test_metrics_details.csv")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('Vus_filtering/test_metrics_analysis_report.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: Vus_filtering/test_metrics_analysis_report.txt") 