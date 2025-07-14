import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# File paths
conc_file = "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Extraction/Concentration_summary_4.7.csv"
status_file = "/bgfs/alee/LO_LAB/Personal/Daisong/APOLLO/Extraction/ILC_APOLLO_status_4.7.csv"

# Read the data
df_conc = pd.read_csv(conc_file)
df_status = pd.read_csv(status_file)

# Print the first few rows of the concentration data
print("First few rows of concentration data:")
print(df_conc.head(10))

# Create a mapping of TP numbers to patient IDs
tp_to_patient = {}
for _, row in df_status.iterrows():
    patient_id = str(row['PT ID']).strip()
    if patient_id == 'Total samples' or pd.isna(patient_id) or patient_id == 'Reason of unreceiving samples':
        continue
    
    # Extract just the number from the patient ID
    patient_number = patient_id.split()[0].strip()
    
    # Process all columns that might contain TP numbers (Draw 1 through Draw 8)
    for i in range(1, 9):
        col_name = f'Draw {i}'
        if col_name in row and isinstance(row[col_name], str) and 'TP' in row[col_name]:
            # Extract just the TPxx-Mxxx part
            tp_parts = row[col_name].split()
            if tp_parts:
                tp_number = tp_parts[0].strip()
                if '(' in tp_number:
                    tp_number = tp_number.split('(')[0].strip()
                tp_to_patient[tp_number] = patient_number

# Print the mapping to debug
print("\nTP to Patient mapping:")
for tp, patient in tp_to_patient.items():
    print(f"{tp} -> {patient}")

# Clean up TP numbers in concentration data
df_conc['Sample'] = df_conc['Sample'].apply(lambda x: x.strip() if isinstance(x, str) else x)

# Add patient IDs to concentration dataframe
df_conc['Patient'] = df_conc['Sample'].map(tp_to_patient)

# Print samples that didn't get mapped
print("\nSamples that didn't get mapped to patients:")
unmapped_samples = df_conc[df_conc['Patient'].isna()]['Sample'].unique()
print(unmapped_samples)

# Print the concentration data to debug
print("\nConcentration data before dropping NaN:")
print(df_conc)

# Drop only rows where Patient is NaN
df_conc = df_conc.dropna(subset=['Patient'])

# Print the concentration data after dropping NaN
print("\nConcentration data after dropping NaN:")
print(df_conc)

# Count samples per patient and sort by count (descending)
patient_counts = df_conc['Patient'].value_counts()
print("\nPatient counts:")
print(patient_counts)

# Create a custom sorting order
# First, separate patients by sample count
patients_by_count = {}
for patient, count in patient_counts.items():
    if count not in patients_by_count:
        patients_by_count[count] = []
    patients_by_count[count].append(patient)

# Create a custom sorted list
sorted_patients = []
# Add patients with more than 3 samples (in descending order)
for count in sorted(patients_by_count.keys(), reverse=True):
    if count > 3:
        sorted_patients.extend(sorted(patients_by_count[count]))

# Add patients with exactly 3 samples
if 3 in patients_by_count:
    sorted_patients.extend(sorted(patients_by_count[3]))

# Add patients 4801, 4749, and 4135 if they have 2 samples
if 2 in patients_by_count:
    priority_patients = ['4801', '4749', '4135']
    other_patients = [p for p in patients_by_count[2] if p not in priority_patients]
    
    # Add priority patients first
    for patient in priority_patients:
        if patient in patients_by_count[2]:
            sorted_patients.append(patient)
    
    # Add other patients with 2 samples
    sorted_patients.extend(sorted(other_patients))

# Add patients with 1 sample
if 1 in patients_by_count:
    sorted_patients.extend(sorted(patients_by_count[1]))

# Sort data by patient sample count (descending) and then by TP number
df_conc['Patient_order'] = df_conc['Patient'].map(dict(zip(sorted_patients, range(len(sorted_patients)))))
df_conc = df_conc.sort_values(['Patient_order', 'Sample'])

# Set up the plot with broken y-axis
plt.switch_backend('Agg')  # Use non-interactive backend
threshold = 1250  # Threshold for breaking the axis
median_line = 840  # Median concentration line

# Create figure with broken y-axis - increased size
fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, figsize=(24, 14), 
                                       gridspec_kw={'height_ratios': [1, 2], 'hspace': 0.05})

# Create x-positions for bars with spacing between patients
x_positions = []
current_x = 0
current_patient = None
for patient in df_conc['Patient']:
    if current_patient != patient:
        if current_patient is not None:
            current_x += 3  # Increased spacing between patients
        current_patient = patient
    x_positions.append(current_x)
    current_x += 1

# Define highlighted patients
highlighted_patients = ['4746', '4308', '4543', '4149', '4591', '4334', '4763', '4571', 'RA125']

# Plot bars on both subplots with different colors for highlighted patients
bars_top = []
bars_bottom = []
for i, (patient, conc) in enumerate(zip(df_conc['Patient'], df_conc['Concentration(ng/mL)'])):
    color = 'darkblue' if patient in highlighted_patients else 'skyblue'
    if conc > threshold:
        # For high concentrations, plot the full bar in both subplots
        bar_top = ax_top.bar(x_positions[i], conc, color=color, width=0.8)
        bar_bottom = ax_bottom.bar(x_positions[i], threshold, color=color, width=0.8)
        bars_top.append(bar_top)
        bars_bottom.append(bar_bottom)
    else:
        # For low concentrations, only plot in the bottom subplot
        bar = ax_bottom.bar(x_positions[i], conc, color=color, width=0.8)
        bars_bottom.append(bar)

# Set y-axis limits
max_conc = df_conc['Concentration(ng/mL)'].max()
ax_top.set_ylim(threshold, max_conc * 1.1)
ax_bottom.set_ylim(0, threshold)

# Add median line to both plots
ax_top.axhline(y=median_line, color='red', linestyle='--', label=f'Median ({median_line} ng/mL)')
ax_bottom.axhline(y=median_line, color='red', linestyle='--', label=f'Median ({median_line} ng/mL)')

# Add value labels on all bars with smaller font
for i, conc in enumerate(df_conc['Concentration(ng/mL)']):
    if conc > threshold:
        ax_top.text(x_positions[i], conc, f'{int(conc)}', ha='center', va='bottom', fontsize=6)
    else:
        ax_bottom.text(x_positions[i], conc, f'{int(conc)}', ha='center', va='bottom', fontsize=6)

# Customize appearance
ax_top.spines['bottom'].set_visible(False)
ax_bottom.spines['top'].set_visible(False)
ax_top.xaxis.tick_top()
ax_top.tick_params(labeltop=False)

# Add broken axis indicators
d = 0.01
kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False)
ax_top.plot((-d, +d), (-d, +d), **kwargs)
ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)
kwargs.update(transform=ax_bottom.transAxes)
ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

# Set labels and title
ax_top.set_title("ctDNA Concentration by Patient", pad=20, fontsize=16)
ax_bottom.set_xlabel("Sample", fontsize=14)
ax_bottom.set_ylabel("Concentration (ng/mL)", fontsize=14)
ax_top.set_ylabel("Concentration (ng/mL)", fontsize=14)

# Set x-ticks and labels (vertical) - reduce font size
plt.xticks(x_positions, df_conc['Sample'], rotation=90, ha='center', fontsize=8)

# Add patient IDs below sample names (moved down a bit more)
for patient in sorted_patients:
    patient_indices = [i for i, p in enumerate(df_conc['Patient']) if p == patient]
    if patient_indices:
        mid_point = np.mean([x_positions[i] for i in patient_indices])
        ax_bottom.text(mid_point, -threshold*0.25, patient, 
                      ha='center', va='top', fontweight='bold', fontsize=10)

# Add legend
ax_bottom.legend(fontsize=10)

# Adjust layout with more space at bottom
plt.subplots_adjust(bottom=0.35)

# Save the plot with lower DPI to reduce complexity
plt.savefig('concentration_plot_4.7.png', dpi=150, bbox_inches='tight')
plt.close() 