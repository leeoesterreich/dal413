import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set the backend to non-interactive 'Agg'
import matplotlib.pyplot as plt

# Read the data file
df = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/Projects/Astrin_proteomic/cancer_plasma>0.1_aggregated.csv')

# Print total number of rows for verification
print(f"Total number of rows in dataset: {len(df)}")

# Replace blank values and NaN with 'Unknown'
df['diag_year'] = df['diag_year'].replace('', np.nan)
df['surg_year'] = df['surg_year'].replace('', np.nan)
df['diag_year'] = pd.to_numeric(df['diag_year'], errors='coerce')
df['surg_year'] = pd.to_numeric(df['surg_year'], errors='coerce')

# Convert years to integers before counting (NaN values will stay as NaN)
df['diag_year'] = df['diag_year'].apply(lambda x: int(x) if pd.notnull(x) else x)
df['surg_year'] = df['surg_year'].apply(lambda x: int(x) if pd.notnull(x) else x)

# Replace NaN with 'Unknown' after converting to integers
df['diag_year'] = df['diag_year'].fillna('Unknown')
df['surg_year'] = df['surg_year'].fillna('Unknown')

# Count the years
diagnosis_counts = df['diag_year'].value_counts()
surgery_counts = df['surg_year'].value_counts()

# Sort the counts (keeping Unknown at the end if it exists)
def sort_counts(counts):
    numeric_years = counts[counts.index != 'Unknown'].sort_index()
    unknown = counts[counts.index == 'Unknown']
    return pd.concat([numeric_years, unknown])

diagnosis_counts = sort_counts(diagnosis_counts)
surgery_counts = sort_counts(surgery_counts)

# Verify total counts
print(f"\nTotal diagnosis years count: {diagnosis_counts.sum()}")
print(f"Total surgery years count: {surgery_counts.sum()}")

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Function to create colors
def create_colors(counts):
    return ['#ff9999' if x == 'Unknown' else '#66b3ff' for x in counts.index]

# Function to format x-axis labels
def format_year_labels(years):
    return [str(int(y)) if isinstance(y, (int, float)) else y for y in years]

# Plot Year of Diagnosis
colors1 = create_colors(diagnosis_counts)
bars1 = ax1.bar(range(len(diagnosis_counts)), diagnosis_counts.values, color=colors1)
ax1.set_title('Year of Diagnosis (Plasma >100µL)')
ax1.set_xlabel('Year')
ax1.set_ylabel('Number of Patients')
ax1.set_xticks(range(len(diagnosis_counts)))
ax1.set_xticklabels(format_year_labels(diagnosis_counts.index), rotation=45, ha='right')

# Add value labels on top of bars
for bar in bars1:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}',
             ha='center', va='bottom')

# Plot Year of Surgery
colors2 = create_colors(surgery_counts)
bars2 = ax2.bar(range(len(surgery_counts)), surgery_counts.values, color=colors2)
ax2.set_title('Year of Surgery/Plasma Extraction (Plasma >100µL)')
ax2.set_xlabel('Year')
ax2.set_ylabel('Number of Patients')
ax2.set_xticks(range(len(surgery_counts)))
ax2.set_xticklabels(format_year_labels(surgery_counts.index), rotation=45, ha='right')

# Add value labels on top of bars
for bar in bars2:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height)}',
             ha='center', va='bottom')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig('year_distribution_bar.png', bbox_inches='tight', dpi=300)
plt.close()

# Print statistics
print("\nYear of Diagnosis Distribution:")
for year, count in diagnosis_counts.items():
    year_str = str(int(year)) if isinstance(year, (int, float)) else year
    print(f"{year_str}: {count} patients")

print("\nYear of Surgery Distribution:")
for year, count in surgery_counts.items():
    year_str = str(int(year)) if isinstance(year, (int, float)) else year
    print(f"{year_str}: {count} patients")

print("\nBar plots have been saved as 'year_distribution_bar.png'") 