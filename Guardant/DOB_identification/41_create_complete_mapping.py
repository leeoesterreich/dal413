import pandas as pd

# Read the master list and cohort files
print("Reading files...")
master_df = pd.read_csv('patient_master_list.csv')
idc_df = pd.read_csv('idc_cohort.csv')
ilc_df = pd.read_csv('ilc_cohort.csv')

# Initialize the result dataframe with master list
result_df = master_df.copy()
result_df['Histology'] = None
result_df['Histo/Behavior ICD-O-3-Desc'] = None

# Process IDC cohort
print("Processing IDC cohort...")
idc_info = idc_df[['GH_ID', 'Histo/Behavior ICD-O-3-Desc']].drop_duplicates()
idc_patients = set(idc_df['GH_ID'].unique())

# Process ILC cohort
print("Processing ILC cohort...")
ilc_info = ilc_df[['GH_ID', 'Histo/Behavior ICD-O-3-Desc']].drop_duplicates()
ilc_patients = set(ilc_df['GH_ID'].unique())

# First mark all IDC patients
result_df.loc[result_df['GH_ID'].isin(idc_patients), 'Histology'] = 'IDC'
# Then override with ILC for those in ILC cohort (as per previous analysis)
result_df.loc[result_df['GH_ID'].isin(ilc_patients), 'Histology'] = 'ILC'

# Add descriptions
idc_desc_dict = dict(zip(idc_info['GH_ID'], idc_info['Histo/Behavior ICD-O-3-Desc']))
ilc_desc_dict = dict(zip(ilc_info['GH_ID'], ilc_info['Histo/Behavior ICD-O-3-Desc']))

# Add descriptions for IDC first
result_df.loc[result_df['GH_ID'].isin(idc_patients), 'Histo/Behavior ICD-O-3-Desc'] = \
    result_df.loc[result_df['GH_ID'].isin(idc_patients), 'GH_ID'].map(idc_desc_dict)

# Override descriptions for ILC patients
result_df.loc[result_df['GH_ID'].isin(ilc_patients), 'Histo/Behavior ICD-O-3-Desc'] = \
    result_df.loc[result_df['GH_ID'].isin(ilc_patients), 'GH_ID'].map(ilc_desc_dict)

# Print statistics
print("\nStatistics:")
print("Total patients: {}".format(len(result_df)))
print("IDC patients: {}".format(sum(result_df['Histology'] == 'IDC')))
print("ILC patients: {}".format(sum(result_df['Histology'] == 'ILC')))

# Save the result
print("\nSaving complete mapping file...")
result_df.to_csv('complete_patient_mapping.csv', index=False)
print("Done!") 