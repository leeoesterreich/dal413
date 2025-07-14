import pandas as pd

# Read the genomic data files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

# Find overlapping patients (those that should be ILC only)
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])

# Remove overlapping patients from IDC data (as they should be ILC)
idc_data = idc_data[~idc_data['Effective Patient ID'].isin(overlapping_ids)]

def analyze_alterations(df, subtype):
    # Group by patient ID
    patient_groups = df.groupby('Effective Patient ID')
    
    # Initialize counters
    total_patients = len(patient_groups)
    patients_with_alterations = 0
    patients_without_alterations = 0
    
    # Store patient details
    patient_details = []
    
    for patient_id, patient_data in patient_groups:
        # Check if any alterations exist
        has_alterations = False
        total_tests = patient_data['Total Test'].max()
        
        # Check each test for alterations
        for _, test in patient_data.iterrows():
            if pd.notna(test['Alteration']):
                has_alterations = True
                break
        
        # Update counters
        if has_alterations:
            patients_with_alterations += 1
        else:
            patients_without_alterations += 1
        
        # Store patient details
        patient_details.append({
            'Effective Patient ID': patient_id,
            'Total Tests': total_tests,
            'Has Alterations': has_alterations,
            'Type': subtype
        })
    
    return total_patients, patients_with_alterations, patients_without_alterations, patient_details

# Analyze IDC patients
print("\nAnalyzing IDC patients:")
print("-" * 50)
idc_total, idc_pos, idc_neg, idc_details = analyze_alterations(idc_data, 'IDC')
print("Total IDC patients: {}".format(idc_total))
print("Patients with alterations: {}".format(idc_pos))
print("Patients without alterations: {}".format(idc_neg))

# Analyze ILC patients
print("\nAnalyzing ILC patients:")
print("-" * 50)
ilc_total, ilc_pos, ilc_neg, ilc_details = analyze_alterations(ilc_data, 'ILC')
print("Total ILC patients: {}".format(ilc_total))
print("Patients with alterations: {}".format(ilc_pos))
print("Patients without alterations: {}".format(ilc_neg))

# Overall statistics
total_patients = idc_total + ilc_total
total_negative = idc_neg + ilc_neg

print("\nOverall Statistics:")
print("-" * 50)
print("Total patients: {}".format(total_patients))
print("Total patients without alterations: {}".format(total_negative))
print("Percentage without alterations: {:.1f}%".format((total_negative/float(total_patients))*100))

# Combine and save patient details
all_details = pd.DataFrame(idc_details + ilc_details)
output_file = 'patient_alterations.csv'
all_details.to_csv(output_file, index=False)
print("\nDetailed patient alteration information saved to: {}".format(output_file))

# Show distribution of tests for patients without alterations
print("\nDistribution of tests for patients without alterations:")
negative_patients = all_details[~all_details['Has Alterations']]
test_dist = negative_patients['Total Tests'].value_counts().sort_index()
print("\nBy number of tests:")
print(test_dist)

# Show sample of patients without alterations and multiple tests
print("\nSample of patients without alterations and multiple tests:")
negative_multi = negative_patients[negative_patients['Total Tests'] > 1].sort_values('Total Tests', ascending=False)
print(negative_multi.head()) 