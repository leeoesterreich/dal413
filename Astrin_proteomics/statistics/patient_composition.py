import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read both files
df_aggregated = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/Projects/Astrin_proteomic/cancer_plasma>0.1_aggregated.csv')
df_original = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/Projects/Astrin_proteomic/Project 2070_cancer.csv')

# Total patients in original dataset
total_patients = 2397

# Current patients (>100µL)
current_patients = len(df_aggregated)

# Patients with >300µL
patients_300ul = len(df_aggregated[df_aggregated['AvailVolume'] >= 0.3])

# Patients with paired tumor (tumor > 0 grams)
patients_with_tumor = len(df_original[df_original['tumor'] > 0])

# Patients with both >100µL and paired tumor
patients_overlap = len(df_aggregated.merge(
    df_original[df_original['tumor'] > 0], 
    on='PatientID', 
    how='inner'
))

# Calculate percentages
percentages = [
    patients_overlap/total_patients*100,
    patients_300ul/total_patients*100,
    current_patients/total_patients*100,
    100  # total patients is 100%
]

# Create figure
plt.figure(figsize=(10, 6))

# Create bars (reversed order)
data = [patients_overlap, patients_300ul, current_patients, total_patients]
labels = ['>100µL & Tumor', '>300µL', '>100µL', 'All Patients']
colors = ['#9b59b6', '#2ecc71', '#3498db', '#95a5a6']  # Purple, Green, Blue, Muted Gray

bars = plt.barh(range(len(data)), data, color=colors)

# Add value labels on the right of bars
for bar, percentage in zip(bars, percentages):
    width = bar.get_width()
    plt.text(width, bar.get_y() + bar.get_height()/2.,
             f'{int(width)} ({percentage:.1f}%)',
             ha='left', va='center')

# Customize the plot
plt.title('Patient Distribution')
plt.xlabel('Number of Patients')
plt.yticks(range(len(data)), labels)

# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig('patient_composition.png', bbox_inches='tight', dpi=300)
plt.close()

# Print statistics
print(f"Total patients in original dataset: {total_patients}")
print(f"Patients with >100µL samples: {current_patients}")
print(f"Patients with >300µL samples: {patients_300ul}")
print(f"Patients with paired tumor: {patients_with_tumor}")
print(f"Patients with both >100µL and paired tumor: {patients_overlap}")

# Calculate percentages
print("\nPercentages:")
print(f">100µL of total: {current_patients/total_patients*100:.1f}%")
print(f">300µL of total: {patients_300ul/total_patients*100:.1f}%")
print(f">300µL of >100µL: {patients_300ul/current_patients*100:.1f}%")
print(f"With tumor of total: {patients_with_tumor/total_patients*100:.1f}%")
print(f"Both >100µL and tumor of total: {patients_overlap/total_patients*100:.1f}%") 