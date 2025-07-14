import matplotlib
matplotlib.use('Agg')  # Set the backend to Agg for headless environments
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Set global font sizes
plt.rcParams['font.size'] = 18  # Base font size
plt.rcParams['axes.titlesize'] = 24  # Title font size
plt.rcParams['axes.labelsize'] = 20  # Axis label font size
plt.rcParams['xtick.labelsize'] = 18  # X-axis tick label font size
plt.rcParams['ytick.labelsize'] = 24  # Y-axis tick label font size (increased from 18)

# ============ Data Preparation ============
file_path = "Concentration_summary_IDC.csv"

# Read CSV file
df = pd.read_csv(file_path)

# Clean sample names by removing parentheses and numbers
df['Samples'] = df['Sample_name'].apply(lambda x: re.sub(r'\s*\([^)]*\)', '', x).strip())

# Handle "<100" values
df['Concentration'] = df['Concentration'].apply(lambda x: 100 if isinstance(x, str) and x == "<100" else float(x))

# Calculate log10 of concentration
df['log10_conc'] = np.log10(df['Concentration'])

# Sort data by concentration
df = df.sort_values('Concentration')

# Create figure (single axis, no broken y-axis)
fig, ax = plt.subplots(figsize=(36, 16))

# Create x positions for bars with spacing at thresholds
x_positions = []
current_pos = 0
spacing = 0.5  # Amount of space to add between threshold groups

for i in range(len(df)):
    if i > 0:  # Check if we need to add spacing
        prev_conc = df['Concentration'].iloc[i-1]
        curr_conc = df['Concentration'].iloc[i]
        # Add spacing if samples cross any threshold
        for threshold in [125, 250, 750, 1000, 1250]:
            if (prev_conc <= threshold and curr_conc > threshold) or \
               (prev_conc > threshold and curr_conc <= threshold):
                current_pos += spacing
                break
    x_positions.append(current_pos)
    current_pos += 1

# Plot bars using log10 values
ax.bar(x_positions, df['log10_conc'], color='skyblue')

# Set y-axis limits (no break)
ax.set_ylim(np.log10(100), 5.2)

# Draw horizontal lines for each threshold
thresholds = [125, 250, 750, 1000, 1250]
yields = [5, 10, 30, 40, 50]  # Corresponding yields in ng
for threshold, yield_val in zip(thresholds, yields):
    log_threshold = np.log10(threshold)
    ax.axhline(y=log_threshold, color='grey', linestyle='dashed', linewidth=1)
    # Add threshold value and yield on the right side
    label = "%d (%dng)" % (threshold, yield_val)
    ax.text(1.01, log_threshold, label, ha='left', va='center', 
            color='black', fontsize=18, transform=ax.get_yaxis_transform())

# Add solid gray line at 840
median_threshold = np.log10(840)
ax.axhline(y=median_threshold, color='grey', linestyle='solid', linewidth=2)
ax.text(1.01, median_threshold, "840 (33.6ng)", ha='left', va='center', 
        color='black', fontsize=18, weight='bold', transform=ax.get_yaxis_transform())

# ============ Annotate Group Counts Only (optional, not per sample) ============
# Define yield groups and their corresponding concentration thresholds
yield_groups = [
    (0, 5, 0, 125),      # 0-5ng: 0-125
    (5, 10, 125, 250),   # 5-10ng: 125-250
    (10, 30, 250, 750),  # 10-30ng: 250-750
    (30, 40, 750, 1000), # 30-40ng: 750-1000
    (40, 50, 1000, 1250),# 40-50ng: 1000-1250
    (50, float('inf'), 1250, float('inf')) # >50ng: >1250
]
group_counts = []
for min_yield, max_yield, min_conc, max_conc in yield_groups:
    count = len(df[(df['Concentration'] > min_conc) & (df['Concentration'] <= max_conc)])
    group_counts.append(count)
group_labels = [
    "0-5ng: {} samples".format(group_counts[0]),
    "5-10ng: {} samples".format(group_counts[1]),
    "10-30ng: {} samples".format(group_counts[2]),
    "30-40ng: {} samples".format(group_counts[3]),
    "40-50ng: {} samples".format(group_counts[4]),
    ">50ng: {} samples".format(group_counts[5])
]
y_positions = [np.log10(62.5), np.log10(187.5), np.log10(500), np.log10(875), np.log10(1125), np.log10(1375)]
for label, y_pos in zip(group_labels, y_positions):
    ax.text(1.08, y_pos, label, ha='left', va='center', 
            color='black', fontsize=20, transform=ax.get_yaxis_transform())

# ============ Final Formatting ============
ax.set_title("IDC ctDNA concentration summary", pad=20)
ax.set_xlabel("Samples", labelpad=10)
ax.set_ylabel("Log10 Concentration (ng/mL)", labelpad=10)

# Rotate x-axis labels
plt.xticks(x_positions, ['']*len(x_positions), rotation=90, fontsize=18)  # No sample names

plt.tight_layout()
plt.savefig('concentration_plot_IDC.png', dpi=300, bbox_inches='tight')
plt.close() 