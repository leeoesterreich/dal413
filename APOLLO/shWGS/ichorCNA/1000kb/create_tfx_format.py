import pandas as pd
import numpy as np

# Read the input files
blood_draw = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/Blood_draw.csv')
summary = pd.read_csv('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/summary.csv')

# Clean up sample names in summary to match blood draw format
summary['Sample'] = summary['Sample'].str.replace('.bam', '')

# Print some debugging information
print("First few rows of summary data:")
print(summary[['Sample', 'Tumor Fraction']].head())
print("\nUnique samples in summary:", len(summary['Sample'].unique()))
print("Total samples in summary:", len(summary))

# Create a dictionary for quick tumor fraction lookup
tf_dict = dict(zip(summary['Sample'], summary['Tumor Fraction']))
print("\nSample of tumor fraction dictionary:")
for k, v in list(tf_dict.items())[:5]:
    print(f"{k}: {v}")

# Create the output dataframe structure
max_draws = blood_draw['# of draws'].max()

# Function to get progression suffix
def get_progression_name(i):
    if i == 1:
        return '1st progression'
    elif i == 2:
        return '2nd progression'
    elif i == 3:
        return '3rd progression'
    else:
        return '%dth progression' % i

# Create column names
columns = ['Patient_Code']
for i in range(1, max_draws + 1):
    columns.extend([get_progression_name(i), 'Tumor_fraction_1000kb bin_%d' % i])

# Process each patient
all_patient_data = []

for _, row in blood_draw.iterrows():
    patient_data = {'Patient_Code': row['PT ID']}
    
    # Process each draw
    for i in range(1, max_draws + 1):
        draw_col = 'Draw %d' % i
        prog_col = get_progression_name(i)
        tf_col = 'Tumor_fraction_1000kb bin_%d' % i
        
        # Initialize with empty values
        patient_data[prog_col] = ''
        patient_data[tf_col] = np.nan
        
        if pd.notna(row[draw_col]) and row[draw_col].strip():
            sample_id = row[draw_col].strip()
            
            # Get tumor fraction from dictionary
            tf_value = tf_dict.get(sample_id, np.nan)
            
            # Add to patient data
            patient_data[prog_col] = sample_id
            patient_data[tf_col] = tf_value
    
    all_patient_data.append(patient_data)

# Create DataFrame from all patient data
output_df = pd.DataFrame(all_patient_data, columns=columns)

# Create header DataFrame with temporary column names
header_df = pd.DataFrame([['Reference genome- Hg38'] + [''] * (len(columns) - 1)], columns=columns)

# Rename the tumor fraction columns to match the desired format
rename_dict = {}
for i in range(1, max_draws + 1):
    old_col = 'Tumor_fraction_1000kb bin_%d' % i
    rename_dict[old_col] = 'Tumor_fraction_1000kb bin'

output_df = output_df.rename(columns=rename_dict)
header_df = header_df.rename(columns=rename_dict)

# Combine DataFrames
final_df = pd.concat([header_df, output_df], ignore_index=True)

# Save to file
output_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/TFx_hg38_1000kb.csv'
final_df.to_csv(output_path, index=False) 