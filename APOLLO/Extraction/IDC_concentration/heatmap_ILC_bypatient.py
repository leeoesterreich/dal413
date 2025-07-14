import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list

# File paths
conc_file = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Extraction/Concentration_summary_4.7.csv"
status_file = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Extraction/ILC_APOLLO_status_4.7.csv"

# Read the data
conc = pd.read_csv(conc_file)
status = pd.read_csv(status_file)

# Clean sample names
conc['Sample_clean'] = conc['Sample'].apply(lambda x: x.strip() if isinstance(x, str) else x)

def clean_sample_name(x):
    return re.sub(r'\s*\([^)]*\)', '', str(x)).strip()

def tp_m_key(s):
    m = re.match(r'TP(\d+)-M(\d+)', s.replace(' ', ''))
    if m:
        return (int(m.group(1)), int(m.group(2)))
    else:
        return (float('inf'), float('inf'))

# Create mapping of TP numbers to patient IDs
sample_to_patient = {}
for _, row in status.iterrows():
    patient_id = str(row['PT ID']).strip()
    if patient_id == 'Total samples' or pd.isna(patient_id) or patient_id == 'Reason of unreceiving samples':
        continue
    patient_number = patient_id.split()[0].strip()
    for i in range(1, 9):
        col_name = f'Draw {i}'
        if col_name in row and isinstance(row[col_name], str) and 'TP' in row[col_name]:
            tp_parts = row[col_name].split()
            if tp_parts:
                tp_number = tp_parts[0].strip()
                if '(' in tp_number:
                    tp_number = tp_number.split('(')[0].strip()
                sample_to_patient[clean_sample_name(tp_number)] = patient_number

# Clean up TP numbers in concentration data
conc['Sample_clean'] = conc['Sample_clean'].apply(clean_sample_name)
conc['Patient'] = conc['Sample_clean'].map(sample_to_patient)

# Drop samples that didn't get mapped
conc = conc.dropna(subset=['Patient'])

# Build patient-sample mapping
patient_samples = {}
for patient in conc['Patient'].unique():
    samples = conc[conc['Patient'] == patient]['Sample_clean'].tolist()
    patient_samples[patient] = samples

# All unique patients (sorted by sample count desc)
all_patients = sorted(conc['Patient'].unique(), key=lambda p: (-len(patient_samples[p]), p))

def get_conc_val(patient, s):
    val = conc[(conc['Patient'] == patient) & (conc['Sample_clean'] == s)]['Concentration(ng/mL)'].values[0]
    if isinstance(val, str) and str(val).strip() == '<100':
        return 100.0
    else:
        return float(val)

# Build the heatmap matrix with gaps and padding
TARGET_SAMPLES = 20  # for ILC
GAP_SIZE = 1  # number of NaN columns between patients
heatmap_rows = []
max_samples = max(TARGET_SAMPLES, max(len(patient_samples[p]) for p in all_patients))
for patient in all_patients:
    samples = patient_samples[patient]
    # Sort chronologically by TP and M number
    samples_sorted = sorted(samples, key=tp_m_key)
    concs_ordered = [get_conc_val(patient, s) for s in samples_sorted]
    # Pad to at least TARGET_SAMPLES
    row = [np.log10(c) for c in concs_ordered] + [np.nan]*(max_samples - len(concs_ordered))
    heatmap_rows.append(row)
    # Add gap after each patient except the last
    if patient != all_patients[-1]:
        heatmap_rows.append([np.nan]*max_samples)

# Convert to numpy array
heatmap_matrix = np.array(heatmap_rows)

# Patient labels (skip gap rows)
patient_labels = []
for i, patient in enumerate(all_patients):
    patient_labels.append(patient)
    if i != len(all_patients)-1:
        patient_labels.append("")  # gap row

# Calculate global min and max for consistent normalization
all_concs = []
for patient in all_patients:
    all_concs.extend([get_conc_val(patient, s) for s in patient_samples[patient]])
vmin = np.log10(min(all_concs))
vmax = np.log10(max(all_concs))

fig, ax = plt.subplots(figsize=(12, 14))
sns.heatmap(
    heatmap_matrix, 
    annot=False, 
    cmap='seismic', 
    center=2.8,
    cbar_kws={'label': 'log10(Concentration)', 'shrink': 0.3, 'pad': 0.02, 'aspect': 10},
    xticklabels=False, 
    yticklabels=patient_labels,
    square=True,
    ax=ax,
    linewidths=0,           # No gridlines
    linecolor=None
)
# Add black rectangles only for cells with data
nrows, ncols = heatmap_matrix.shape
for i in range(nrows):
    for j in range(ncols):
        if not np.isnan(heatmap_matrix[i, j]):
            rect = Rectangle((j, i), 1, 1, fill=False, edgecolor='black', linewidth=1.5)
            ax.add_patch(rect)
ax.set_xlabel('Samples')
ax.set_ylabel('Patient')
fig.suptitle('ILC ctDNA concentration by patient/sample (log10 scale)', y=0.98)
plt.tight_layout()
plt.savefig('heatmap_ILC_bypatient.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics
print("\nSummary statistics by patient:")
print(conc.groupby('Patient')['Concentration(ng/mL)'].agg(['count', 'mean', 'median', 'std'])) 