import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os # Added for creating output directory

# Read the XGBoost results
# Path to the CSV file
csv_file_path = 'results_xgb/summaries/all_signatures_results_xgb.csv'
if not os.path.exists(csv_file_path):
    print(f"Error: CSV file not found at {csv_file_path}")
    exit()
df = pd.read_csv(csv_file_path)

# Columns to use for AUCs
auc_columns = ['auc_median', 'auc_tercile', 'auc_quartile']

# Remove any rows with NaN values in the relevant AUC columns
df = df.dropna(subset=auc_columns)

print(f"Loaded {len(df)} complete XGBoost results from {csv_file_path}")

if df.empty:
    print("No data available for plotting after NaN removal.")
    exit()

# Create a figure with a smaller width
plt.figure(figsize=(8, 6))

# Prepare data for plotting
plot_data = pd.DataFrame({
    'Median AUC (Test Set)': df['auc_median'],
    'Tercile AUC (Test Set)': df['auc_tercile'],
    'Quartile AUC (Test Set)': df['auc_quartile']
})

# Create box plot with blank filling and narrower boxes
sns.boxplot(data=plot_data, 
            palette=['white', 'white', 'white'],
            width=0.3,
            showfliers=False)  # Hide default outliers

# Add individual points in light grey
sns.stripplot(data=plot_data, 
              color='#808080',  # Medium grey color
              size=4,
              alpha=0.4,
              jitter=0.2)  # Add some jitter to better separate overlapping points

# Add horizontal line at y=0.75
plt.axhline(y=0.75, color='gray', linestyle='--', alpha=0.5)

# Highlight the three specific signatures requested (if present in the 'signature' column)
# Note: These AUCs are all from the test split of the training phase.
if 'signature' in df.columns:
highlight_signatures = {
    'UNC_RB_LOH': ('red', -0.05),
    'Basal signaling': ('blue', 0),
    'Estrogen signaling': ('green', 0.05)
}

legend_handles = []
legend_labels = []

for target_sig, (color, offset) in highlight_signatures.items():
        # Find the signature in the dataframe - need to search for partial matches in the 'signature' column
        # The 'signature' column in all_signatures_results_xgb.csv should contain the original signature names.
    matching_rows = df[df['signature'].str.contains(target_sig.replace(' ', '_'), case=False, na=False)]
    
    if not matching_rows.empty:
        row = matching_rows.iloc[0]  # Take the first match
            original_signature_name = row['signature'] # The actual name from the CSV
            values = [row['auc_median'], row['auc_tercile'], row['auc_quartile']]
        for i, value in enumerate(values):
                if pd.notna(value):
            point = plt.plot(i + offset, value, 'o', color=color, markersize=8, alpha=0.7, markeredgecolor='black')[0]
            if i == 0:  # Only add to legend once per signature
                legend_handles.append(point)
                        legend_labels.append(target_sig) # Use the user-friendly target_sig for legend
            print(f"Found and highlighted: {target_sig} (matched as '{original_signature_name}')")
        else:
            print(f"Warning: Could not find signature matching '{target_sig}' for highlighting.")
    else:
    print("Warning: 'signature' column not found in CSV. Skipping signature highlighting.")
    legend_handles = [] # Ensure it exists for the legend call
    legend_labels = []

# Customize the plot
plt.title('XGBoost Model AUC Comparison (Test Set Performance)', fontsize=14, pad=20)
plt.ylabel('AUC Score', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)

# Add legend
if legend_handles:
    plt.legend(legend_handles, legend_labels, loc='upper right', bbox_to_anchor=(1, 1))

# Set y-axis limits to ensure visibility of all points
min_val = plot_data.min().min()
max_val = plot_data.max().max()
if pd.isna(min_val) or pd.isna(max_val):
    min_val, max_val = 0.4, 1.0 # Default if data is empty
else:
    min_val = min(0.4, min_val - 0.05)
    max_val = max(1.0, max_val + 0.05)
plt.ylim(min_val, max_val)

# Adjust layout
plt.tight_layout()

# Save the plot
output_dir = 'results_xgb/summaries'
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'xgb_auc_comparison_boxplot.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved as {output_path}")

# Print summary statistics
print(f"\n=== XGBoost AUC Summary Statistics (from Test Set during Training) ===")
for auc_col, col_name in zip(auc_columns, plot_data.columns):
    if auc_col in df.columns:
        aucs = df[auc_col].dropna().values
        if len(aucs) > 0:
            print(f"{col_name:<25} - Mean: {np.mean(aucs):.3f}, Median: {np.median(aucs):.3f}, Std: {np.std(aucs):.3f}, Count: {len(aucs)}")
        else:
            print(f"{col_name:<25} - No data")
    else:
        print(f"{col_name:<25} - Column not found in CSV") 