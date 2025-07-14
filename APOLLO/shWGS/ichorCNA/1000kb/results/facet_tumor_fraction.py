import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Read the data files
blood_draw = pd.read_csv('Blood_draw.csv')
summary = pd.read_csv('summary.csv')

# Clean up sample names in summary to match blood draw format
summary['Sample'] = summary['Sample'].str.replace('.bam', '')

# Prepare patient data
patient_trajectories = {}
for idx, row in blood_draw.iterrows():
    patient_id = row['PT ID']
    if patient_id == '/(RA)':
        continue
    draw_columns = [col for col in row.index if col.startswith('Draw')]
    samples = [row[col] for col in draw_columns]
    samples = [s.strip() for s in samples if pd.notna(s)]
    sample_data = summary[summary['Sample'].isin(samples)]
    data_points = []
    for draw_idx, sample in enumerate(samples, 1):
        tf_values = sample_data[sample_data['Sample'] == sample]['Tumor Fraction'].values
        if len(tf_values) > 0:
            data_points.append((draw_idx, tf_values[0], sample))
    if data_points:
        data_points.sort(key=lambda x: x[0])
        time_points = [d[0] for d in data_points]
        tumor_fractions = [d[1] for d in data_points]
        sample_names = [d[2] for d in data_points]
        patient_trajectories[patient_id] = (time_points, tumor_fractions, sample_names)

# Sort patients by number of draws (descending)
patient_ids = sorted(patient_trajectories.keys(), key=lambda pid: len(patient_trajectories[pid][0]), reverse=True)
n_patients = len(patient_ids)

# Assign colors: red for most draws, purple for least, using rainbow spectrum
colors = plt.cm.rainbow(np.linspace(1, 0, n_patients))
color_map = {pid: color for pid, color in zip(patient_ids, colors)}

# Determine grid size (try to make it as square as possible)
n_cols = int(np.ceil(np.sqrt(n_patients)))
n_rows = int(np.ceil(n_patients / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(3*n_cols, 4*n_rows), sharex=False, sharey=True)
axes = axes.flatten()

for ax, patient_id in zip(axes, patient_ids):
    time_points, tumor_fractions, sample_names = patient_trajectories[patient_id]
    ax.plot(time_points, tumor_fractions, 'o-', lw=2, color=color_map[patient_id])
    ax.set_title(f'Patient {patient_id}', fontsize=10)
    ax.set_xlabel('Blood Draw')
    ax.set_ylabel('Tumor Fraction')
    ax.set_ylim(0, 0.7)
    ax.axhline(y=0.03, color='grey', linestyle='--', alpha=0.5)
    # Set x-ticks to actual blood draw numbers
    ax.set_xticks(time_points)
    ax.set_xticklabels([str(tp) for tp in time_points], rotation=0, fontsize=8)
    ax.set_xlim(min(time_points)-0.5, max(time_points)+0.5)

# Hide unused subplots
for i in range(len(patient_ids), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
output_dir = 'tumor_fraction'
os.makedirs(output_dir, exist_ok=True)
plt.savefig(os.path.join(output_dir, 'facet_tumor_fraction.png'), dpi=300)
plt.close()

print('Facet plot saved as facet_tumor_fraction.png') 