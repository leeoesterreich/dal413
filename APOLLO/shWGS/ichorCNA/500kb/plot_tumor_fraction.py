import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

# Read the data files
blood_draw = pd.read_csv('results/Blood_draw.csv')
summary = pd.read_csv('results/summary.csv')

# Clean up sample names in summary to match blood draw format
summary['Sample'] = summary['Sample'].str.replace('.bam', '')

# Function to plot individual patient data
def plot_patient_data(patient_id, patient_row, summary_data):
    # Get all draw columns (Draw 1 through Draw 7)
    draw_columns = [col for col in patient_row.index if col.startswith('Draw')]
    samples = [patient_row[col] for col in draw_columns]
    samples = [s.strip() for s in samples if pd.notna(s)]  # Remove any whitespace
    
    # Get tumor fractions for these samples
    sample_data = summary_data[summary_data['Sample'].isin(samples)]
    if len(sample_data) == 0:
        return None, None
    
    # Create pairs of (time_point, tumor_fraction) and filter out any missing data
    data_points = []
    for idx, sample in enumerate(samples, 1):  # Start counting at 1 for first draw
        tf_values = sample_data[sample_data['Sample'] == sample]['Tumor Fraction'].values
        if len(tf_values) > 0:
            data_points.append((idx, tf_values[0]))
    
    if not data_points:
        return None, None
    
    # Sort by time point (should already be sorted, but just in case)
    data_points.sort(key=lambda x: x[0])
    
    # Separate into time points and tumor fractions
    time_points = [d[0] for d in data_points]
    tumor_fractions = [d[1] for d in data_points]
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(time_points, tumor_fractions, 'o-')
    plt.title(f'Patient {patient_id} - Tumor Fraction Over Time')
    plt.xlabel('Blood Draw Number')
    plt.ylabel('Tumor Fraction')
    plt.grid(True)
    plt.xticks(time_points)  # Show all time points
    plt.tight_layout()
    plt.savefig(f'results/patient_{patient_id}_tumor_fraction.png')
    plt.close()
    
    return time_points, tumor_fractions

# Store data for each patient
patient_data = []

# Collect data for each patient
for idx, row in blood_draw.iterrows():
    patient_id = row['PT ID']
    if patient_id == '/(RA)':  # Skip the special case
        continue
        
    time_points, fractions = plot_patient_data(patient_id, row, summary)
    if time_points and fractions:
        patient_data.append({
            'patient_id': patient_id,
            'time_points': time_points,
            'fractions': fractions
        })

if patient_data:
    # Create integrated plot
    plt.figure(figsize=(15, 10))
    
    # Use a color palette with enough distinct colors
    colors = plt.cm.rainbow(np.linspace(0, 1, len(patient_data)))
    
    # Plot each patient's line with a different color
    for patient, color in zip(patient_data, colors):
        plt.plot(patient['time_points'], patient['fractions'], 'o-', 
                color=color, label=f'Patient {patient["patient_id"]}',
                linewidth=2, markersize=8)
    
    plt.title('Individual Patient Tumor Fraction Changes Over Time')
    plt.xlabel('Blood Draw Number')
    plt.ylabel('Tumor Fraction')
    plt.grid(True)
    
    # Adjust legend to be outside the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # Adjust layout to accommodate the legend
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    
    plt.savefig('results/integrated_tumor_fraction.png', bbox_inches='tight')
    plt.close()
    
    # Save the data in a more readable format
    with open('results/patient_trajectories.txt', 'w') as f:
        f.write('Patient ID\tDraw Number\tTumor Fraction\n')
        for patient in patient_data:
            for t, frac in zip(patient['time_points'], patient['fractions']):
                f.write(f"{patient['patient_id']}\t{t}\t{frac:.4f}\n")

    # Calculate and print statistics
    with open('results/trend_stats.txt', 'w') as f:
        f.write('Trend Statistics:\n\n')
        f.write('Time Point\tMean\tStd Dev\tN\n')
        for t in range(max(max(p['time_points']) for p in patient_data)):
            values = [frac for p in patient_data for frac in p['fractions'] if t+1 in p['time_points']]
            if len(values) > 0:
                f.write(f'{t+1}\t{np.mean(values):.4f}\t{np.std(values):.4f}\t{len(values)}\n') 