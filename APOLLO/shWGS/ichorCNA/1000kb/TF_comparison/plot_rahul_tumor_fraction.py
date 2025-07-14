import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import os
import glob

# Define paths
rahul_base_dir = "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/ILC/ichorCNA/hg38_1000kb/results/ichorCNA"
blood_draw_file = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/shWGS/ichorCNA/1000kb/results/Blood_draw.csv"
output_dir = 'rahul_tumor_fraction'
os.makedirs(output_dir, exist_ok=True)

def get_tumor_fraction(sample_id):
    params_file = os.path.join(rahul_base_dir, sample_id, f"{sample_id}.params.txt")
    if not os.path.exists(params_file):
        return None
    
    with open(params_file, 'r') as f:
        for line in f:
            if line.startswith("Tumor Fraction:"):
                return float(line.strip().split("\t")[1])
    return None

# Read the blood draw data
blood_draw = pd.read_csv(blood_draw_file)

# Function to plot individual patient data
def plot_patient_data(patient_id, patient_row):
    # Get all draw columns (Draw 1 through Draw 7)
    draw_columns = [col for col in patient_row.index if col.startswith('Draw')]
    samples = [patient_row[col] for col in draw_columns]
    samples = [s.strip() for s in samples if pd.notna(s)]  # Remove any whitespace
    
    # Get tumor fractions for these samples
    data_points = []
    for idx, sample in enumerate(samples, 1):  # Start counting at 1 for first draw
        tf = get_tumor_fraction(sample)
        if tf is not None:
            data_points.append((idx, tf))
    
    if not data_points:
        return None, None
    
    # Sort by time point
    data_points.sort(key=lambda x: x[0])
    
    # Separate into time points and tumor fractions
    time_points = [d[0] for d in data_points]
    tumor_fractions = [d[1] for d in data_points]
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(time_points, tumor_fractions, 'o-')
    # Add horizontal line at TF=3%
    plt.axhline(y=0.03, color='grey', linestyle='--', alpha=0.5)
    plt.title(f'Patient {patient_id} - Tumor Fraction Over Time (Rahul\'s Analysis)')
    plt.xlabel('Blood Draw Number')
    plt.ylabel('Tumor Fraction')
    plt.grid(True)
    plt.xticks(time_points)  # Show all time points
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'patient_{patient_id}_tumor_fraction.png'))
    plt.close()
    
    return time_points, tumor_fractions

# Store data for each patient
patient_data = []

# Collect data for each patient
for idx, row in blood_draw.iterrows():
    patient_id = row['PT ID']
    if patient_id == '/(RA)':  # Skip the special case
        continue
        
    time_points, fractions = plot_patient_data(patient_id, row)
    if time_points and fractions:
        patient_data.append({
            'patient_id': patient_id,
            'time_points': time_points,
            'fractions': fractions
        })

if patient_data:
    # Sort patients by number of draws (descending)
    patient_data.sort(key=lambda x: len(x['time_points']), reverse=True)
    
    # Create integrated plot
    plt.figure(figsize=(15, 10))
    
    # Use a color palette with enough distinct colors, reversed so red is for more draws
    colors = plt.cm.rainbow(np.linspace(1, 0, len(patient_data)))
    
    # Plot each patient's line with a different color
    for patient, color in zip(patient_data, colors):
        plt.plot(patient['time_points'], patient['fractions'], 'o-', 
                color=color, label=f'Patient {patient["patient_id"]}',
                linewidth=2, markersize=8)
    
    # Add horizontal line at TF=3%
    plt.axhline(y=0.03, color='grey', linestyle='--', alpha=0.5)
    
    plt.title('Individual Patient Tumor Fraction Changes Over Time (Rahul\'s Analysis)')
    plt.xlabel('Blood Draw Number')
    plt.ylabel('Tumor Fraction')
    plt.grid(True)
    
    # Adjust legend to be outside the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # Adjust layout to accommodate the legend
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    
    plt.savefig(os.path.join(output_dir, 'integrated_tumor_fraction.png'), bbox_inches='tight')
    plt.close()
    
    # Save the data in a more readable format
    with open(os.path.join(output_dir, 'patient_trajectories.txt'), 'w') as f:
        f.write('Patient ID\tDraw Number\tTumor Fraction\n')
        for patient in patient_data:
            for t, frac in zip(patient['time_points'], patient['fractions']):
                f.write(f"{patient['patient_id']}\t{t}\t{frac:.4f}\n")

    # Calculate and print statistics
    with open(os.path.join(output_dir, 'trend_stats.txt'), 'w') as f:
        f.write('Trend Statistics:\n\n')
        f.write('Time Point\tMean\tStd Dev\tN\n')
        for t in range(max(max(p['time_points']) for p in patient_data)):
            values = [frac for p in patient_data for frac in p['fractions'] if t+1 in p['time_points']]
            if len(values) > 0:
                f.write(f'{t+1}\t{np.mean(values):.4f}\t{np.std(values):.4f}\t{len(values)}\n')

print("Analysis complete. Check the rahul_tumor_fraction directory for output files.") 