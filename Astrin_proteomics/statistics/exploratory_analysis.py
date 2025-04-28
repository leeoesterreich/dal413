import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/Projects/Astrin_proteomic/Project 2070_cancer.csv')

# Group by PatientID and check if each patient has any non-zero values in plasma/ds_plasma and tumor
# Also check if AvailVolume > 0.1 for plasma/ds_plasma entries
patient_summary = df.groupby('PatientID').agg({
    'plasma': lambda x: ((x != 0) & (df.loc[x.index, 'AvailVolume'] > 0.1)).any(),
    'ds_plasma': lambda x: ((x != 0) & (df.loc[x.index, 'AvailVolume'] > 0.1)).any(),
    'tumor': lambda x: (x != 0).any()
}).reset_index()

# Count patients with non-zero values in both plasma/ds_plasma and tumor
plasma_tumor_patients = patient_summary[
    ((patient_summary['plasma']) | (patient_summary['ds_plasma'])) & 
    (patient_summary['tumor'])
]['PatientID'].nunique()

# Count patients with non-zero values in either plasma or ds_plasma (AvailVolume > 0.1)
plasma_ds_plasma_patients = patient_summary[
    (patient_summary['plasma']) | (patient_summary['ds_plasma'])
]['PatientID'].nunique()

# Count patients with non-zero values in tumor
tumor_patients = patient_summary[patient_summary['tumor']]['PatientID'].nunique()

# Print summary
print(f"Total number of unique patients: {len(patient_summary)}")
print(f"Number of patients with non-zero values in either plasma or ds_plasma (AvailVolume > 0.1): {plasma_ds_plasma_patients}")
print(f"Number of patients with non-zero values in tumor: {tumor_patients}")
print(f"Number of patients with non-zero values in both plasma/ds_plasma (AvailVolume > 0.1) and tumor: {plasma_tumor_patients}")

# Read the aggregated file for stage analysis
df_aggregated = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/Projects/Astrin_proteomic/cancer_plasma>0.1_aggregated.csv')

# Replace blank stages with 'Unavailable'
df_aggregated['Stage'] = df_aggregated['Stage'].fillna('Unavailable')

# Count stages
stage_counts = df_aggregated['Stage'].value_counts()

# Create pie chart
plt.figure(figsize=(10, 8))
plt.pie(stage_counts, labels=stage_counts.index, autopct='%1.1f%%')
plt.title('Distribution of Patient Stages')
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle

# Save the plot
plt.savefig('patient_stages_pie.png')
plt.close()

# Print stage counts
print("\nStage Distribution:")
for stage, count in stage_counts.items():
    print(f"{stage}: {count} patients ({count/len(df_aggregated)*100:.1f}%)")

# Total patients in original dataset
total_patients = 2397

# Current patients (>100µL)
current_patients = len(df_aggregated)

# Patients with >300µL
patients_300ul = len(df_aggregated[df_aggregated['AvailVolume'] >= 0.3])

# Patients with paired tumor
patients_with_tumor = len(df_aggregated[df_aggregated['tumor'] > 0])

# Patients with both >300µL and paired tumor
patients_300ul_and_tumor = len(df_aggregated[(df_aggregated['AvailVolume'] >= 0.3) & (df_aggregated['tumor'] > 0)])

# Create figure
plt.figure(figsize=(12, 4))

# Create bars
bars = plt.barh(['Patient\nComposition'], [total_patients], color='lightgray', label='All Patients')
plt.barh(['Patient\nComposition'], [current_patients], color='#66b3ff', label='>100µL')

# Calculate positions for the overlapping bars
tumor_only = patients_with_tumor - patients_300ul_and_tumor
vol_300_only = patients_300ul - patients_300ul_and_tumor

# Add >300µL and tumor bars
plt.barh(['Patient\nComposition'], [patients_300ul], left=[current_patients-patients_300ul], 
         color='#99ff99', label='>300µL')
plt.barh(['Patient\nComposition'], [patients_with_tumor], left=[current_patients-patients_with_tumor], 
         color='#ff9999', label='With Tumor', alpha=0.7)

# Add text annotations
plt.text(total_patients/2, 0, f'Total: {total_patients}', 
         ha='center', va='center')
plt.text(current_patients-100, 0, f'>100µL: {current_patients}', 
         ha='right', va='center')
plt.text(current_patients+50, 0, 
         f'>300µL: {patients_300ul}\nWith Tumor: {patients_with_tumor}\nBoth: {patients_300ul_and_tumor}',
         ha='left', va='center')

# Customize the plot
plt.xlabel('Number of Patients')
plt.title('Patient Sample Composition')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig('patient_composition.png', bbox_inches='tight', dpi=300)
plt.close()

# Print statistics
print(f"Total patients in original dataset: {total_patients}")
print(f"Patients with >100µL samples: {current_patients}")
print(f"Patients with >300µL samples: {patients_300ul}")
print(f"Patients with paired tumor: {patients_with_tumor}")
print(f"Patients with both >300µL and paired tumor: {patients_300ul_and_tumor}")

# Calculate percentages
print("\nPercentages:")
print(f">100µL of total: {current_patients/total_patients*100:.1f}%")
print(f">300µL of >100µL: {patients_300ul/current_patients*100:.1f}%")
print(f"With tumor of >100µL: {patients_with_tumor/current_patients*100:.1f}%")
print(f"Both >300µL and tumor of >100µL: {patients_300ul_and_tumor/current_patients*100:.1f}%") 