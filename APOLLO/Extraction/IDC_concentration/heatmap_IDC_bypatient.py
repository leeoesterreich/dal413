import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list

# File paths
conc_file = "Concentration_summary_IDC.csv"
patient_file = "Patient_Information_IDC.csv"

def clean_sample_name(x):
    return re.sub(r'\s*\([^)]*\)', '', str(x)).strip()

def tp_m_key(s):
    m = re.match(r'TP(\d+)-M(\d+)', s.replace(' ', ''))
    if m:
        return (int(m.group(1)), int(m.group(2)))
    else:
        return (float('inf'), float('inf'))

# Read the data
conc = pd.read_csv(conc_file)
pat = pd.read_csv(patient_file)
conc['Sample_clean'] = conc['Sample_name'].apply(clean_sample_name)

# Build patient-sample mapping
patients = []
for idx, row in pat.iterrows():
    sample_names = []
    for col in pat.columns:
        if col.startswith('Draw'):
            val = row[col]
            if pd.notnull(val) and str(val).strip():
                sample_names.append(clean_sample_name(val))
    ra_status = str(row['RA status']).strip()
    patients.append({'samples': sample_names, 'ra_status': ra_status})

# For each patient, collect available samples and concentrations
patient_data = []
for p in patients:
    samples = [s for s in p['samples'] if s in set(conc['Sample_clean'])]
    if samples:
        # Sort chronologically by TP and M number
        samples_sorted = sorted(samples, key=tp_m_key)
        concs = []
        for s in samples_sorted:
            val = conc[conc['Sample_clean'] == s]['Concentration'].values[0]
            if isinstance(val, str) and val.strip() == '<100':
                concs.append(100.0)
            else:
                concs.append(float(val))
        patient_data.append({
            'samples': samples_sorted,
            'concs': concs,
            'ra_status': p['ra_status'],
            'n_samples': len(samples_sorted)
        })

# Sort patients by number of available samples (descending)
patient_data = sorted(patient_data, key=lambda x: x['n_samples'], reverse=True)

# Build the heatmap matrix with gaps and padding
TARGET_SAMPLES = 40  # for IDC
max_samples = max(TARGET_SAMPLES, max(p['n_samples'] for p in patient_data))
heatmap_rows = []
patient_labels = []
for i, pdata in enumerate(patient_data):
    concs_ordered = pdata['concs']
    row = [np.log10(c) for c in concs_ordered] + [np.nan]*(max_samples - len(concs_ordered))
    heatmap_rows.append(row)
    patient_labels.append(f"P{i+1}")
    # Add a gap row after each patient except the last
    if i != len(patient_data) - 1:
        heatmap_rows.append([np.nan]*max_samples)
        patient_labels.append("")

heatmap_matrix = np.array(heatmap_rows)

# Calculate global min and max for consistent normalization
all_concs = []
for pdata in patient_data:
    all_concs.extend(pdata['concs'])
vmin = np.log10(min(all_concs))
vmax = np.log10(max(all_concs))

# Calculate shape for dynamic figsize
nrows, ncols = heatmap_matrix.shape
scale = 0.5  # Increase for bigger boxes
fig_width = ncols * scale
fig_height = nrows * scale
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
for i in range(nrows):
    for j in range(ncols):
        if not np.isnan(heatmap_matrix[i, j]):
            rect = Rectangle((j, i), 1, 1, fill=False, edgecolor='black', linewidth=0.7)
            ax.add_patch(rect)
ax.set_xlabel('Samples')
ax.set_ylabel('Patient')
fig.suptitle('IDC ctDNA concentration by patient/sample (log10 scale)', y=0.98)
plt.tight_layout()
plt.savefig('heatmap_IDC_bypatient.png', dpi=300, bbox_inches='tight')
plt.close() 