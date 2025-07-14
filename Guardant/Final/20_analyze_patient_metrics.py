import pandas as pd
import numpy as np
from io import StringIO
import sys

# Create string buffer to capture output
output_buffer = StringIO()
original_stdout = sys.stdout
sys.stdout = output_buffer

def clean_tumor_fraction(value):
    if pd.isna(value):
        return np.nan
    if isinstance(value, str):
        value = value.strip()
        if value == '<0.05' or value == '<0.05%':
            return 0.0
        if value == '>40' or value == '>40%':
            return 41.0
        # Remove '%' if present
        value = value.rstrip('%')
        # Handle any other '>' or '<' cases
        if value.startswith('>'):
            return float(value[1:])
        if value.startswith('<'):
            return float(value[1:])
        return float(value)
    return float(value)

def analyze_patient_metrics(df, cohort_name):
    print(f"\n{cohort_name} Patient Analysis:")
    print("-" * (len(cohort_name) + 18))
    
    # Get unique patients
    unique_patients = df.drop_duplicates(subset=['Effective Patient ID'])
    total_patients = len(unique_patients)
    print(f"\nTotal Unique Patients: {total_patients}")
    
    # 1. MSI Status Analysis
    print("\nMSI Status Distribution:")
    msi_status = unique_patients['MSI High'].fillna('NOT DETECTED')
    msi_counts = msi_status.map(lambda x: 'MSI-H' if x == 'DETECTED' else 'MSS').value_counts()
    for status, count in msi_counts.items():
        percentage = (count / total_patients) * 100
        print(f"  {status}: {count} patients ({percentage:.1f}%)")
    
    # 2. Tumor Fraction Analysis
    print("\nTumor Fraction Analysis:")
    # Convert percentage strings to float and handle special cases
    tumor_fractions = unique_patients['Tumor Fraction'].apply(clean_tumor_fraction)
    valid_tf = tumor_fractions.dropna()
    if len(valid_tf) > 0:
        median_tf = np.median(valid_tf)
        print(f"  Median Tumor Fraction: {median_tf:.2f}%")
        print(f"  Patients with valid Tumor Fraction: {len(valid_tf)}")
        print(f"  Patients with Tumor Fraction <0.05%: {len(valid_tf[valid_tf == 0])}")
        print(f"  Patients with Tumor Fraction >40%: {len(valid_tf[valid_tf == 41])}")
    else:
        print("  No valid Tumor Fraction values found")
    
    # Save tumor fraction data to CSV
    tumor_data = pd.DataFrame({
        'Patient_ID': unique_patients['Effective Patient ID'],
        'Original_Tumor_Fraction': unique_patients['Tumor Fraction'],
        'Processed_Tumor_Fraction': tumor_fractions
    }).sort_values('Patient_ID')
    tumor_data.to_csv(f'{cohort_name}_tumor_fractions_corrected.csv', index=False)
    print(f"  Tumor fraction data saved to: {cohort_name}_tumor_fractions_corrected.csv")
    
    # 3. TMB Score Analysis
    print("\nTMB Score Analysis:")
    # Convert to numeric and handle any non-numeric values
    tmb_scores = pd.to_numeric(unique_patients['TMB Score'], errors='coerce')
    valid_tmb = tmb_scores.dropna()
    if len(valid_tmb) > 0:
        median_tmb = np.median(valid_tmb)
        print(f"  Median TMB Score: {median_tmb:.2f}")
        print(f"  Patients with valid TMB Score: {len(valid_tmb)}")
    else:
        print("  No valid TMB Score values found")

# Read the Infinity data files
print("Reading Infinity data files...")
idc_infinity = pd.read_csv('IDC_genomic_infinity_corrected.csv')
ilc_infinity = pd.read_csv('ILC_genomic_infinity_corrected.csv')

# Analyze each cohort
analyze_patient_metrics(idc_infinity, "IDC")
analyze_patient_metrics(ilc_infinity, "ILC")

# Analyze combined cohort
print("\nCombined IDC and ILC Analysis:")
print("-" * 28)
combined_infinity = pd.concat([idc_infinity, ilc_infinity])
analyze_patient_metrics(combined_infinity, "Combined")

# Restore original stdout and save the output
sys.stdout = original_stdout
output_text = output_buffer.getvalue()

# Save to file
with open('patient_metrics_analysis_corrected.txt', 'w') as f:
    f.write(output_text)

print("Full analysis report saved to: patient_metrics_analysis_corrected.txt") 