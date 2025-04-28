import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the aggregated file
df = pd.read_csv('/bgfs/alee/LO_LAB/Personal/Daisong/Projects/Astrin_proteomic/cancer_plasma>0.1_aggregated.csv')

# Print total number of rows
total_rows = len(df)
print(f"Total number of patients: {total_rows}")

# Check for different types of missing values
print("\nMissing value types in stage column:")
print(df['stage'].unique())

# Replace all types of missing/blank values and '0' with '0(NA)'
df['stage'] = df['stage'].replace('', '0(NA)')  # Empty strings
df['stage'] = df['stage'].replace('0', '0(NA)')  # Stage 0
df['stage'] = df['stage'].fillna('0(NA)')      # NA/NaN values

# Clean up stage names (remove extra spaces)
df['stage'] = df['stage'].str.strip()

# Count stages
stage_counts = df['stage'].value_counts(dropna=False)  # Include NA values in counting

# Print stage distribution
print("\nStage Distribution:")
for stage, count in stage_counts.items():
    print(f"{stage}: {count} patients ({count/len(df)*100:.1f}%)")

# Create pie chart
plt.figure(figsize=(10, 8))
colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99', '#ff99cc', '#99ccff']

# Create custom autopct function to show both count and percentage
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return f'{val}\n({pct:.1f}%)'
    return my_autopct

plt.pie(stage_counts.values, labels=stage_counts.index, 
        autopct=make_autopct(stage_counts.values),
        colors=colors, startangle=90,
        labeldistance=1.05)  # Bring labels even closer to the pie chart
plt.title('Distribution of Patient Stages (Plasma >100µL)')
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle

# Save the plot
plt.savefig('stage_distribution_pie.png', bbox_inches='tight', dpi=300)
plt.close()

# Verify total
print(f"\nTotal (verification): {stage_counts.sum()} patients")
print("\nPie chart has been saved as 'stage_distribution_pie.png'") 