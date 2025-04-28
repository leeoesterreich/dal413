import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the CSV file
df = pd.read_csv('Project 2070_normal.csv')

# Count samples with non-zero plasma
non_zero_plasma = df[df['plasma'] != 0]
print(f"\nTotal number of samples with non-zero plasma: {len(non_zero_plasma)}")

# Filter rows where plasma is not zero
plasma_df = df[df['plasma'] != 0]

# 1. Bar chart for AvailVolume and uom
plt.figure(figsize=(10, 6))
volume_uom_counts = plasma_df.groupby(['AvailVolume', 'uom']).size().unstack(fill_value=0)
ax = volume_uom_counts.plot(kind='bar', stacked=True)
plt.title('Distribution of Available Volume and Units')
plt.xlabel('Available Volume')
plt.ylabel('Count')
plt.legend(title='Unit of Measure')
plt.xticks(rotation=45)

# Add total count annotations on top of each bar
for i, row in enumerate(volume_uom_counts.itertuples()):
    total = sum(row[1:])  # Sum all values except the index
    ax.text(i, total, f'n={total}', ha='center', va='bottom')

plt.tight_layout()
plt.savefig('volume_distribution.png')
plt.close()

# 2. Pie chart for Diagnosis
# Standardize the Donor category
plasma_df['Diagnosis'] = plasma_df['Diagnosis'].replace('DONOR', 'Donor')
plt.figure(figsize=(12, 8))
diagnosis_counts = plasma_df['Diagnosis'].value_counts()
plt.pie(diagnosis_counts, labels=diagnosis_counts.index, autopct='%1.1f%%')
plt.title('Distribution of Diagnoses')
plt.tight_layout()
plt.savefig('diagnosis_distribution.png')
plt.close()

# Print some summary statistics
print("\nSummary of filtered data:")
print(f"Total number of samples with plasma: {len(plasma_df)}")
print("\nDiagnosis distribution:")
print(diagnosis_counts)
print("\nVolume and unit distribution:")
print(volume_uom_counts) 