# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys

# --- Configuration ---
FULL_GENOMIC_FILE = 'serial_genomic_with_patient_id.csv'
ILC_POSITIVE_COHORT_FILE = 'ILC_genomic_corrected.csv'
IDC_POSITIVE_COHORT_FILE = 'IDC_genomic_corrected.csv'
ORIGINAL_ILC_FILE = 'ilc_cohort.csv' # For histology description
ORIGINAL_IDC_FILE = 'idc_cohort.csv' # For histology description
PATIENT_ID_COLUMN = 'Effective Patient ID'
OUTPUT_FILE = 'test_counts_for_positive_cohort.csv'

# --- Data Loading and Preparation ---
try:
    print("Reading data files for the positive cohort analysis...")
    df_full_genomic = pd.read_csv(FULL_GENOMIC_FILE, encoding='latin-1')
    df_ilc_cohort = pd.read_csv(ILC_POSITIVE_COHORT_FILE, encoding='latin-1')
    df_idc_cohort = pd.read_csv(IDC_POSITIVE_COHORT_FILE, encoding='latin-1')
    # Load original files to get histology description
    df_ilc_orig = pd.read_csv(ORIGINAL_ILC_FILE, usecols=['GH_ID', 'Histo/Behavior ICD-O-3-Desc'])
    df_idc_orig = pd.read_csv(ORIGINAL_IDC_FILE, usecols=['GH_ID', 'Histo/Behavior ICD-O-3-Desc'])
    print("...done")
except FileNotFoundError as e:
    print(f"Error: A required file was not found - {e}", file=sys.stderr)
    sys.exit(1)

# Get the sets of patient IDs from the positive cohorts
ilc_ids = set(df_ilc_cohort[PATIENT_ID_COLUMN].unique())
idc_ids = set(df_idc_cohort[PATIENT_ID_COLUMN].unique())

# Filter the main genomic file to get ALL tests for ONLY the positive patients
ilc_genomic = df_full_genomic[df_full_genomic[PATIENT_ID_COLUMN].isin(ilc_ids)]
idc_genomic = df_full_genomic[df_full_genomic[PATIENT_ID_COLUMN].isin(idc_ids)]

# --- Merge Histology Information ---
# Create a mapping from GH_ID to histology for the original files
ilc_histo_map = df_ilc_orig.drop_duplicates('GH_ID').set_index('GH_ID')['Histo/Behavior ICD-O-3-Desc']
idc_histo_map = df_idc_orig.drop_duplicates('GH_ID').set_index('GH_ID')['Histo/Behavior ICD-O-3-Desc']

# Apply the mapping to our filtered dataframes
ilc_genomic = ilc_genomic.assign(**{'Histo/Behavior ICD-O-3-Desc': ilc_genomic['GH_ID'].map(ilc_histo_map)})
idc_genomic = idc_genomic.assign(**{'Histo/Behavior ICD-O-3-Desc': idc_genomic['GH_ID'].map(idc_histo_map)})

def analyze_test_counts(df, cohort_name):
    # Group by patient ID and get their pre-calculated Total Test value
    test_counts = df.groupby(PATIENT_ID_COLUMN)['Total Test'].first()
    
    # Get histology types
    histology_types = df.groupby(PATIENT_ID_COLUMN)['Histo/Behavior ICD-O-3-Desc'].first()
    
    # Basic statistics
    total_patients = len(test_counts)
    if total_patients == 0:
        print(f"\n{cohort_name} Cohort Analysis:")
        print("No patients found in this cohort.")
        return None, None

    # Correctly sum the 'Total Test' value for each unique patient
    total_tests = test_counts.sum()
    avg_tests = total_tests / total_patients
    
    # Patients with multiple tests
    multiple_tests = test_counts[test_counts > 1]
    multiple_test_median = multiple_tests.median() if not multiple_tests.empty else 0
    
    print(f"\n{cohort_name} Cohort Analysis:")
    print("-" * (len(cohort_name) + 16))
    print(f"Total unique patients: {total_patients}")
    print(f"Total number of tests (sum of 'Total Test' per patient): {total_tests}")
    print(f"Average tests per patient: {avg_tests:.2f}")
    print(f"\nPatients with multiple tests (Total Test > 1):")
    print(f"Number of patients: {len(multiple_tests)} ({len(multiple_tests)/total_patients*100:.1f}%)")
    print(f"Median number of tests for these patients: {multiple_test_median:.1f}")
    
    # Distribution of test counts
    print("\nDistribution of Total Test values:")
    print(test_counts.value_counts().sort_index())
    
    return test_counts, histology_types

# Analyze both cohorts
print("\nAnalyzing test counts for the cohort with detected alterations...")
idc_counts, idc_histology = analyze_test_counts(idc_genomic, "IDC")
ilc_counts, ilc_histology = analyze_test_counts(ilc_genomic, "ILC")

# Save detailed patient-level data
output_data = {
    'Patient_ID': [],
    'Cohort': [],
    'Total_Tests': [],
    'Histology_Description': []
}

if idc_counts is not None:
for patient, tests in idc_counts.items():
    output_data['Patient_ID'].append(patient)
    output_data['Cohort'].append('IDC')
    output_data['Total_Tests'].append(tests)
        output_data['Histology_Description'].append(idc_histology[patient])

if ilc_counts is not None:
for patient, tests in ilc_counts.items():
    output_data['Patient_ID'].append(patient)
    output_data['Cohort'].append('ILC')
    output_data['Total_Tests'].append(tests)
        output_data['Histology_Description'].append(ilc_histology[patient])

# Create and save detailed output
output_df = pd.DataFrame(output_data)
output_df = output_df.sort_values(['Cohort', 'Total_Tests', 'Patient_ID'])
output_df.to_csv(OUTPUT_FILE, index=False)
print(f"\nDetailed patient-level data saved to: {OUTPUT_FILE}") 