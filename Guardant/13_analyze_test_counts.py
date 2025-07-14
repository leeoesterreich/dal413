# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

# Read the data files
idc_genomic = pd.read_csv('IDC_genomic_with_histology.csv', encoding='latin-1')
ilc_genomic = pd.read_csv('ILC_genomic_with_histology.csv', encoding='latin-1')

# Print column names to verify
print("\nColumns in IDC genomic data:")
print(idc_genomic.columns.tolist())
print("\nColumns in ILC genomic data:")
print(ilc_genomic.columns.tolist())

def analyze_test_counts(df, cohort_name):
    # Group by patient ID and get their Total Test values
    gh_id_col = 'GH_ID'  # Updated column name
    test_counts = df.groupby(gh_id_col)['Total Test'].first()  # Taking first since it should be same for each patient
    
    # Get histology types
    histology_types = df.groupby(gh_id_col)['Histology Type'].first()
    
    # Basic statistics
    total_patients = len(test_counts)
    total_tests = test_counts.sum()
    avg_tests = total_tests / total_patients
    
    # Patients with multiple tests
    multiple_tests = test_counts[test_counts > 1]
    multiple_test_median = multiple_tests.median()
    
    print(f"\n{cohort_name} Cohort Analysis:")
    print("-" * (len(cohort_name) + 16))
    print(f"Total unique patients: {total_patients}")
    print(f"Total number of tests: {total_tests}")
    print(f"Average tests per patient: {avg_tests:.2f}")
    print(f"\nPatients with multiple tests (Total Test > 1):")
    print(f"Number of patients: {len(multiple_tests)} ({len(multiple_tests)/total_patients*100:.1f}%)")
    print(f"Median number of tests: {multiple_test_median:.1f}")
    
    # Distribution of test counts
    print("\nDistribution of Total Test values:")
    print(test_counts.value_counts().sort_index())
    
    return test_counts, histology_types

# Analyze both cohorts
print("Analyzing test counts in genomic data...")
idc_counts, idc_histology = analyze_test_counts(idc_genomic, "IDC")
ilc_counts, ilc_histology = analyze_test_counts(ilc_genomic, "ILC")

# Save detailed patient-level data
output_data = {
    'Patient_ID': [],
    'Cohort': [],
    'Total_Tests': [],
    'Histology_Type': []
}

for patient, tests in idc_counts.items():
    output_data['Patient_ID'].append(patient)
    output_data['Cohort'].append('IDC')
    output_data['Total_Tests'].append(tests)
    output_data['Histology_Type'].append(idc_histology[patient])

for patient, tests in ilc_counts.items():
    output_data['Patient_ID'].append(patient)
    output_data['Cohort'].append('ILC')
    output_data['Total_Tests'].append(tests)
    output_data['Histology_Type'].append(ilc_histology[patient])

# Create and save detailed output
output_df = pd.DataFrame(output_data)
output_df = output_df.sort_values(['Cohort', 'Total_Tests', 'Patient_ID'])
output_df.to_csv('test_counts_by_patient.csv', index=False)
print("\nDetailed patient-level data saved to: test_counts_by_patient.csv") 