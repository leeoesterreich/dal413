import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

def load_clinical_data():
    """Load clinical data"""
    clinical_data = pd.read_csv('data_clinical_patient.txt', sep='\t', skiprows=4)
    clinical_data.set_index('PATIENT_ID', inplace=True)
    return clinical_data

def load_scores():
    """Load both direct and predicted signature scores"""
    direct_scores = pd.read_csv('direct_signature_scores.csv', index_col=0)
    predicted_scores = pd.read_csv('predicted_signature_scores.csv', index_col=0)
    return direct_scores, predicted_scores

def plot_survival_curves(clinical_data, scores, score_type, output_file):
    """Plot survival curves for high and low risk groups"""
    # Merge clinical data with scores
    data = clinical_data.join(scores)
    
    # Calculate cutoff for top 1/3
    cutoff_score = data['score'].quantile(2/3)
    print(f"\nTop 1/3 cutoff {score_type} score: {cutoff_score}")
    
    # Assign risk groups
    data['risk_group'] = np.where(data['score'] > cutoff_score, 'High Risk', 'Low Risk')
    
    # Create Kaplan-Meier fitter
    kmf = KaplanMeierFitter()
    
    # Plot survival curves
    plt.figure(figsize=(10, 6))
    
    # High risk group
    high_risk = data[data['risk_group'] == 'High Risk']
    kmf.fit(high_risk['OS_MONTHS'], high_risk['OS_STATUS'], label='High Risk')
    kmf.plot(ci_show=True)
    
    # Low risk group
    low_risk = data[data['risk_group'] == 'Low Risk']
    kmf.fit(low_risk['OS_MONTHS'], low_risk['OS_STATUS'], label='Low Risk')
    kmf.plot(ci_show=True)
    
    # Add labels and title
    plt.xlabel('Time (months)')
    plt.ylabel('Survival Probability')
    plt.title(f'Survival Curves by {score_type} Risk Group (Top 1/3)')
    plt.grid(True)
    
    # Save plot
    plt.savefig(output_file)
    plt.close()
    
    # Print group sizes
    print(f"\nNumber of samples in each risk group ({score_type}):")
    print(data['risk_group'].value_counts())

def main():
    # Load data
    clinical_data = load_clinical_data()
    direct_scores, predicted_scores = load_scores()
    
    # Filter clinical data to include only samples with valid OS_MONTHS
    clinical_data = clinical_data[clinical_data['OS_MONTHS'].notna() & (clinical_data['OS_MONTHS'] > 0)]
    clinical_data['OS_MONTHS'] = pd.to_numeric(clinical_data['OS_MONTHS'], errors='coerce')
    clinical_data = clinical_data.dropna(subset=['OS_MONTHS'])
    
    # Convert OS_STATUS to binary (1 for events, 0 for censored)
    clinical_data['OS_STATUS'] = (clinical_data['OS_STATUS'] == '1:DECEASED').astype(int)
    
    # Plot survival curves for direct scores
    direct_scores = direct_scores.rename(columns={direct_scores.columns[0]: 'score'})
    plot_survival_curves(clinical_data, direct_scores, 'Direct Signature', 'direct_survival_curves.png')
    
    # Plot survival curves for predicted scores
    predicted_scores = predicted_scores.rename(columns={predicted_scores.columns[0]: 'score'})
    plot_survival_curves(clinical_data, predicted_scores, 'Predicted Signature', 'predicted_survival_curves.png')
    
    print("\nSurvival analysis complete. Plots have been saved as:")
    print("- direct_survival_curves.png")
    print("- predicted_survival_curves.png")

if __name__ == "__main__":
    main() 