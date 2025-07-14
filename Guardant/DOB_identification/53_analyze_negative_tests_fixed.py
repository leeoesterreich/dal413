import pandas as pd

# Read the genomic data files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Find overlapping patients (those that should be ILC only)
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])

# Remove overlapping patients from IDC data (as they should be ILC)
idc_data = idc_data[~idc_data['Effective Patient ID'].isin(overlapping_ids)]

def analyze_test_results(df, subtype):
    # First, let's look at the unique values in the Alteration Detected field
    print(f"\nUnique values in 'Aleration Detected?' field for {subtype}:")
    print(df['Aleration Detected?'].unique())
    
    # Group by patient and analyze their test results
    patient_results = df.groupby('Effective Patient ID').agg({
        'Aleration Detected?': ['count', lambda x: (x == 'False').all()],
        'Total Test': 'max'
    }).reset_index()
    
    # Rename columns
    patient_results.columns = ['Effective Patient ID', 'Total Tests', 'All Negative', 'Max Tests']
    
    # Count statistics
    total_patients = len(patient_results)
    all_negative = len(patient_results[patient_results['All Negative']])
    
    # Get breakdown by number of tests for negative patients
    negative_by_tests = patient_results[patient_results['All Negative']]['Max Tests'].value_counts().sort_index()
    
    return total_patients, all_negative, negative_by_tests, patient_results

# Analyze IDC patients
print("\nAnalyzing IDC patients:")
print("-" * 50)
idc_total, idc_negative, idc_dist, idc_results = analyze_test_results(idc_data, 'IDC')
print(f"Total IDC patients: {idc_total}")
print(f"IDC patients with all negative tests: {idc_negative}")
print("\nDistribution of tests for all-negative IDC patients:")
print(idc_dist)

# Analyze ILC patients
print("\nAnalyzing ILC patients:")
print("-" * 50)
ilc_total, ilc_negative, ilc_dist, ilc_results = analyze_test_results(ilc_data, 'ILC')
print(f"Total ILC patients: {ilc_total}")
print(f"ILC patients with all negative tests: {ilc_negative}")
print("\nDistribution of tests for all-negative ILC patients:")
print(ilc_dist)

# Overall statistics
total_patients = idc_total + ilc_total
total_negative = idc_negative + ilc_negative

print("\nOverall Statistics:")
print("-" * 50)
print(f"Total patients: {total_patients}")
print(f"Total patients with all negative tests: {total_negative}")
print(f"Percentage with all negative tests: {(total_negative/total_patients)*100:.1f}%")

# Add type indicator and combine results
idc_results['Type'] = 'IDC'
ilc_results['Type'] = 'ILC'
all_results = pd.concat([idc_results, ilc_results])

# Save detailed information
output_file = 'patient_test_results_fixed.csv'
all_results.to_csv(output_file, index=False)
print(f"\nDetailed patient test results saved to: {output_file}")

# Show sample of patients with all negative tests and multiple tests
print("\nSample of patients with all negative tests and multiple tests:")
negative_multi = all_results[(all_results['All Negative']) & (all_results['Max Tests'] > 1)]
negative_multi = negative_multi.sort_values('Max Tests', ascending=False)
print(negative_multi.head()) 