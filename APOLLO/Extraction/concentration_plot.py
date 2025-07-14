import matplotlib
matplotlib.use('Agg')  # Set the backend to Agg for headless environments
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============ Data Preparation ============
file_path = "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Extraction/Concentration_summary_4.7.xlsx"

# Read Excel file
df = pd.read_excel(file_path, sheet_name="Sheet2", usecols=[0, 1], names=["Samples", "Concentration"])
df.dropna(inplace=True)  # Drop NaN values

# Exclude specific samples
df = df[~df['Samples'].isin(['TP19-M82', 'TP20-M445'])]

# Sort data by concentration
df = df.sort_values('Concentration')

# Create figure with subplots
fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, figsize=(25, 8))

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

# Plot bars on both subplots
ax_top.bar(x_positions, df['Concentration'], color='skyblue')
ax_bottom.bar(x_positions, df['Concentration'], color='skyblue')

# Set limits for the subplots
ax_top.set_ylim(1250, df['Concentration'].max() * 1.1)
ax_bottom.set_ylim(0, 1250)

# Hide spines between subplots
ax_top.spines['bottom'].set_visible(False)
ax_bottom.spines['top'].set_visible(False)

# Move x-axis ticks to bottom subplot only
ax_top.xaxis.tick_top()
ax_top.tick_params(labeltop=False)
ax_bottom.xaxis.tick_bottom()

# ============ Drawing "Break" Marks ============
d = 0.01  # Size of diagonal lines
kwargs = dict(color='k', clip_on=False)

# Diagonal lines on the top subplot (bottom edge)
ax_top.plot((-d, +d), (-d, +d), transform=ax_top.transAxes, **kwargs)
ax_top.plot((1 - d, 1 + d), (-d, +d), transform=ax_top.transAxes, **kwargs)

# Diagonal lines on the bottom subplot (top edge)
ax_bottom.plot((-d, +d), (1 - d, 1 + d), transform=ax_bottom.transAxes, **kwargs)
ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax_bottom.transAxes, **kwargs)

# ============ Draw Broken Grey Lines and Add Threshold Values ============
# Draw horizontal lines for each threshold
thresholds = [125, 250, 750, 1000, 1250]
yields = [5, 10, 30, 40, 50]  # Corresponding yields in ng

for threshold, yield_val in zip(thresholds, yields):
    if threshold <= 1250:
        ax_bottom.axhline(y=threshold, color='grey', linestyle='dashed', linewidth=1)
        # Add threshold value and yield on the right side
        label = "%d (%dng)" % (threshold, yield_val)
        ax_bottom.text(1.01, threshold, label, ha='left', va='center', 
                      color='black', fontsize=10, transform=ax_bottom.get_yaxis_transform())
    else:
        ax_top.axhline(y=threshold, color='grey', linestyle='dashed', linewidth=1)
        # Add threshold value and yield on the right side
        label = "%d (%dng)" % (threshold, yield_val)
        ax_top.text(1.01, threshold, label, ha='left', va='center', 
                   color='black', fontsize=10, transform=ax_top.get_yaxis_transform())

# Add solid gray line at 840
ax_bottom.axhline(y=840, color='grey', linestyle='solid', linewidth=2)
ax_bottom.text(1.01, 840, "840 (33.6ng)", ha='left', va='center', 
              color='black', fontsize=10, weight='bold', transform=ax_bottom.get_yaxis_transform())

# ============ Annotate Values and Count Samples ============
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

# Count samples in each group
for min_yield, max_yield, min_conc, max_conc in yield_groups:
    # Count samples
    count = len(df[(df['Concentration'] > min_conc) & (df['Concentration'] <= max_conc)])
    group_counts.append(count)

# Add group count annotations
group_labels = [
    "0-5ng: {} samples".format(group_counts[0]),
    "5-10ng: {} samples".format(group_counts[1]),
    "10-30ng: {} samples".format(group_counts[2]),
    "30-40ng: {} samples".format(group_counts[3]),
    "40-50ng: {} samples".format(group_counts[4]),
    ">50ng: {} samples".format(group_counts[5])
]

# Position for group labels (on the right side of the plot)
y_positions = [62.5, 187.5, 500, 875, 1125, 1375]  # Middle of each group
for label, y_pos in zip(group_labels, y_positions):
    if y_pos <= 1250:
        ax_bottom.text(1.08, y_pos, label, ha='left', va='center', 
                      color='black', fontsize=8, transform=ax_bottom.get_yaxis_transform())
    else:
        ax_top.text(1.08, y_pos, label, ha='left', va='center', 
                   color='black', fontsize=8, transform=ax_top.get_yaxis_transform())

# Annotate individual values
for i, (sample, conc) in enumerate(zip(df['Samples'], df['Concentration'])):
    if conc > 1250:
        ax_top.text(x_positions[i], conc, str(int(conc)), ha='center', va='bottom', 
                   color='black', fontsize=8)
    else:
        ax_bottom.text(x_positions[i], conc, str(int(conc)), ha='center', va='bottom', 
                      color='black', fontsize=8)

# ============ Final Formatting ============
ax_top.set_title("ctDNA Concentration Summary (Sorted by Concentration)")
ax_bottom.set_xlabel("Samples")
ax_top.set_ylabel("Concentration (ng/mL)")
ax_bottom.set_ylabel("Concentration (ng/mL)")

# Rotate x-axis labels
plt.xticks(x_positions, df['Samples'], rotation=90, fontsize=10)

# Tight layout to minimize overlaps
plt.tight_layout()

# Save the figure
plt.savefig('concentration_plot.png', dpi=300, bbox_inches='tight')
plt.close() 