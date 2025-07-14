import pandas as pd
import numpy as np

# Read the genomic data files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Find overlapping patients (those that should be ILC only)
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])

# Remove overlapping patients from IDC data (as they should be ILC)
idc_data = idc_data[~idc_data['Effective Patient ID'].isin(overlapping_ids)]

def analyze_testing_data(df, subtype):
    # Group by patient ID and get their testing info
    patient_data = df.groupby('Effective Patient ID').agg({
        'Total Test': 'max',  # Get the total number of tests for each patient
        'Serial Test': ['nunique', 'max'],  # Get number of unique serial tests and max serial test number
        'GH_ID': 'first',
        'Patient DOB': 'first'
    }).reset_index()
    
    # Rename columns for clarity
    patient_data.columns = ['Effective Patient ID', 'Total Tests', 'Unique Serial Tests', 'Max Serial Number', 'GH_ID', 'Patient DOB']
    
    # Count patients with multiple tests
    total_patients = len(patient_data)
    patients_with_multiple = len(patient_data[patient_data['Total Tests'] > 1])
    
    # Get distribution of total tests
    test_distribution = patient_data['Total Tests'].value_counts().sort_index()
    
    return patient_data, total_patients, patients_with_multiple, test_distribution

# Analyze IDC patients
print("\nAnalyzing IDC patients:")
print("-" * 50)
idc_patient_data, idc_total, idc_multiple, idc_dist = analyze_testing_data(idc_data, 'IDC')
print(f"Total unique IDC patients: {idc_total}")
print(f"IDC patients with multiple tests: {idc_multiple}")
print("\nDistribution of total tests for IDC:")
print(idc_dist)

# Analyze ILC patients
print("\nAnalyzing ILC patients:")
print("-" * 50)
ilc_patient_data, ilc_total, ilc_multiple, ilc_dist = analyze_testing_data(ilc_data, 'ILC')
print(f"Total unique ILC patients: {ilc_total}")
print(f"ILC patients with multiple tests: {ilc_multiple}")
print("\nDistribution of total tests for ILC:")
print(ilc_dist)

# Add type indicator
idc_patient_data['Type'] = 'IDC'
ilc_patient_data['Type'] = 'ILC'

# Combine all patient data
all_patient_data = pd.concat([idc_patient_data, ilc_patient_data])

# Overall statistics
total_patients = idc_total + ilc_total
total_multiple = idc_multiple + ilc_multiple

print("\nOverall Statistics:")
print("-" * 50)
print(f"Total unique patients (IDC + ILC): {total_patients}")
print(f"Total patients with multiple tests: {total_multiple}")
print(f"Percentage with multiple tests: {(total_multiple/total_patients)*100:.1f}%")

# Save detailed patient information
output_file = 'patient_testing_details_fixed.csv'
all_patient_data.to_csv(output_file, index=False)
print(f"\nDetailed patient testing information saved to: {output_file}")

# Show sample of patients with multiple tests
print("\nSample of patients with multiple tests:")
multi_test_patients = all_patient_data[all_patient_data['Total Tests'] > 1].sort_values('Total Tests', ascending=False)
print(multi_test_patients.head()) 