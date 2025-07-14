# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os

# --- Configuration ---
IDC_INPUT_PATH = 'final_genomic_cohorts/IDC_genomic.csv'
ILC_INPUT_PATH = 'final_genomic_cohorts/ILC_genomic.csv'
OUTPUT_DIR = 'final_genomic_cohorts'
DETAILED_OUTPUT_FILENAME = os.path.join(OUTPUT_DIR, 'test_counts_by_patient.csv')
SUMMARY_OUTPUT_FILENAME = os.path.join(OUTPUT_DIR, 'test_counts_summary.csv')

# --- Main Functions ---
def analyze_test_counts(df, cohort_name):
    """Analyzes test counts for a given dataframe and returns key statistics."""
    gh_id_col = 'GH_ID'
    
    # Check if 'Total Test' column exists. If not, calculate it.
    if 'Total Test' in df.columns and df['Total Test'].notna().any():
        print(f"Using existing 'Total Test' column for {cohort_name} cohort.")
        # Ensure we handle multiple rows per test correctly by taking the first value per patient
        test_counts = df.groupby(gh_id_col)['Total Test'].first()
    else:
        print(f"'Total Test' column not found or empty for {cohort_name} cohort. Calculating from unique dates.")
        # Ensure the date column is in datetime format to count unique dates properly
        if 'Sample Received Date' in df.columns:
            df['Sample Received Date'] = pd.to_datetime(df['Sample Received Date'], errors='coerce')
            test_counts = df.groupby(gh_id_col)['Sample Received Date'].nunique()
        else:
            print(f"Error: Neither 'Total Test' nor 'Sample Received Date' found for {cohort_name}.")
            return None, None # Return None if essential columns are missing

    if test_counts.empty:
        print(f"No test counts could be determined for {cohort_name}.")
        return None, None

    # Get histology types for detailed report
    histology_col = 'hist_Histo/Behavior ICD-O-3-Desc'
    histology_types = df.groupby(gh_id_col)[histology_col].first() if histology_col in df.columns else pd.Series(None, index=test_counts.index)
    
    # --- Calculate Statistics ---
    total_patients = len(test_counts)
    total_tests = int(test_counts.sum())
    avg_tests = total_tests / total_patients if total_patients > 0 else 0
    
    multiple_tests_series = test_counts[test_counts > 1]
    patients_with_multiple_tests = len(multiple_tests_series)
    percentage_multiple_tests = (patients_with_multiple_tests / total_patients * 100) if total_patients > 0 else 0
    median_num_tests = multiple_tests_series.median() if patients_with_multiple_tests > 0 else 0
    
    # --- Print to Console ---
    print(f"\n--- {cohort_name} Cohort Analysis ---")
    print(f"Total unique patients: {total_patients}")
    print(f"Total number of tests: {total_tests}")
    print(f"Average tests per patient: {avg_tests:.2f}")
    print(f"Patients with multiple (serial) tests: {patients_with_multiple_tests} ({percentage_multiple_tests:.1f}%)")
    if patients_with_multiple_tests > 0:
        print(f"Median number of tests for patients with serial testing: {median_num_tests:.1f}")
    
    # --- Prepare Summary Data for CSV ---
    summary_data = {
        'Cohort': cohort_name,
        'Total Unique Patients': total_patients,
        'Total Tests': total_tests,
        'Patients with Serial Testing (>1 test)': patients_with_multiple_tests,
        'Median # of Tests (for serial testing patients)': median_num_tests
    }
    
    return test_counts, histology_types, summary_data

def main():
    """Main function to run the analysis."""
    try:
        idc_genomic = pd.read_csv(IDC_INPUT_PATH, encoding='latin-1', low_memory=False)
        ilc_genomic = pd.read_csv(ILC_INPUT_PATH, encoding='latin-1', low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: Input file not found. Ensure '{e.filename}' exists.")
        return

    # Combine dataframes for total analysis
    combined_genomic = pd.concat([idc_genomic, ilc_genomic], ignore_index=True)

    # Analyze each cohort and the combined total
    all_summary_data = []
    
    idc_counts, idc_histology, idc_summary = analyze_test_counts(idc_genomic, "IDC")
    ilc_counts, ilc_histology, ilc_summary = analyze_test_counts(ilc_genomic, "ILC")
    _, _, combined_summary = analyze_test_counts(combined_genomic, "Combined Total")
    
    if idc_summary: all_summary_data.append(idc_summary)
    if ilc_summary: all_summary_data.append(ilc_summary)
    if combined_summary: all_summary_data.append(combined_summary)

    # --- Save Summary CSV ---
    if all_summary_data:
        summary_df = pd.DataFrame(all_summary_data)
        summary_df.to_csv(SUMMARY_OUTPUT_FILENAME, index=False)
        print(f"\nSummary statistics saved to: {SUMMARY_OUTPUT_FILENAME}")
    else:
        print("\nNo summary data was generated.")

    # --- Save Detailed Patient-Level CSV ---
    if idc_counts is not None and ilc_counts is not None:
        output_data = {'Patient_ID': [], 'Cohort': [], 'Total_Tests': [], 'Histology_Description': []}
        for patient, tests in idc_counts.items():
            output_data['Patient_ID'].append(patient)
            output_data['Cohort'].append('IDC')
            output_data['Total_Tests'].append(tests)
            output_data['Histology_Description'].append(idc_histology.get(patient))

        for patient, tests in ilc_counts.items():
            output_data['Patient_ID'].append(patient)
            output_data['Cohort'].append('ILC')
            output_data['Total_Tests'].append(tests)
            output_data['Histology_Description'].append(ilc_histology.get(patient))

        detailed_df = pd.DataFrame(output_data).sort_values(['Cohort', 'Total_Tests', 'Patient_ID'])
        detailed_df.to_csv(DETAILED_OUTPUT_FILENAME, index=False)
        print(f"\nDetailed patient-level data saved to: {DETAILED_OUTPUT_FILENAME}")
    else:
        print("\nCould not generate detailed patient-level data due to earlier errors.")

if __name__ == "__main__":
    main() 