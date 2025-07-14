import pandas as pd

# Read the genomic data files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Find overlapping patients (those that should be ILC only)
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])

# Remove overlapping patients from IDC data (as they should be ILC)
idc_data = idc_data[~idc_data['Effective Patient ID'].isin(overlapping_ids)]

# Function to analyze patient testing data
def analyze_patient_data(df, subtype):
    # Count unique patients
    total_patients = df['Effective Patient ID'].nunique()
    
    # Get patients with serial testing
    serial_patients = df[df['Serial Test'] > 1]['Effective Patient ID'].nunique()
    
    # Get distribution of serial test numbers
    serial_dist = df.groupby('Effective Patient ID')['Serial Test'].max().value_counts().sort_index()
    
    return total_patients, serial_patients, serial_dist

# Analyze both datasets
print("\nAnalyzing IDC patients:")
print("-" * 50)
idc_total, idc_serial, idc_dist = analyze_patient_data(idc_data, 'IDC')
print(f"Total unique IDC patients: {idc_total}")
print(f"IDC patients with serial testing: {idc_serial}")
print("\nDistribution of serial tests for IDC:")
print(idc_dist)

print("\nAnalyzing ILC patients:")
print("-" * 50)
ilc_total, ilc_serial, ilc_dist = analyze_patient_data(ilc_data, 'ILC')
print(f"Total unique ILC patients: {ilc_total}")
print(f"ILC patients with serial testing: {ilc_serial}")
print("\nDistribution of serial tests for ILC:")
print(ilc_dist)

# Combined statistics
total_unique_patients = idc_total + ilc_total
total_serial_patients = idc_serial + ilc_serial

print("\nOverall Statistics:")
print("-" * 50)
print(f"Total unique patients (IDC + ILC): {total_unique_patients}")
print(f"Total patients with serial testing: {total_serial_patients}")
print(f"Percentage with serial testing: {(total_serial_patients/total_unique_patients)*100:.1f}%")

# Detailed analysis of serial testing patterns
print("\nDetailed Serial Testing Analysis:")
print("-" * 50)

def get_patient_testing_details(df):
    patient_tests = df.groupby('Effective Patient ID').agg({
        'Serial Test': ['max', 'count'],
        'GH_ID': 'first',
        'Patient DOB': 'first'
    }).reset_index()
    patient_tests.columns = ['Effective Patient ID', 'Max Serial Number', 'Total Tests', 'GH_ID', 'Patient DOB']
    return patient_tests

# Get testing details for both types
idc_details = get_patient_testing_details(idc_data)
ilc_details = get_patient_testing_details(ilc_data)

# Add type indicator
idc_details['Type'] = 'IDC'
ilc_details['Type'] = 'ILC'

# Combine all patient details
all_patient_details = pd.concat([idc_details, ilc_details])

# Save detailed patient testing information
output_file = 'patient_testing_details.csv'
all_patient_details.to_csv(output_file, index=False)
print(f"\nDetailed patient testing information saved to: {output_file}")

# Show sample of patients with multiple tests
print("\nSample of patients with multiple tests:")
multi_test_patients = all_patient_details[all_patient_details['Total Tests'] > 1].sort_values('Total Tests', ascending=False)
print(multi_test_patients.head()) 