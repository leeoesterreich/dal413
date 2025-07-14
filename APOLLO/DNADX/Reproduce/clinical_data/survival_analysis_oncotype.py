import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def load_clinical_data():
    """Load clinical data"""
    clinical_data = pd.read_csv('data_clinical_patient.txt', sep='\t', skiprows=4)
    clinical_data.set_index('PATIENT_ID', inplace=True)
    return clinical_data

def load_oncotype_scores():
    """Load Oncotype DX scores"""
    print("Loading Oncotype DX scores...")
    scores = pd.read_csv('oncotype_dx_score.csv', index_col=0)
    return scores

def plot_survival_curves(clinical_data, scores, output_file='oncotype_survival_curves.png'):
    """Plot survival curves for high and low risk groups based on Oncotype DX scores"""
    # Merge clinical data with scores
    data = clinical_data.join(scores)

    # Calculate cutoff for top 1/3
    cutoff_score = data['GHI_RS_Model_NJEM.2004_PMID.15591335'].quantile(2/3)
    print(f"\nCutoff score (top 1/3): {cutoff_score}")

    # Assign risk groups
    data['risk_group'] = 'Low'
    data.loc[data['GHI_RS_Model_NJEM.2004_PMID.15591335'] > cutoff_score, 'risk_group'] = 'High'

    # Convert OS_STATUS to binary (1 for event, 0 for censored)
    data['event'] = (data['OS_STATUS'] == '1:DECEASED').astype(int)

    # Convert OS_MONTHS to OS_YEARS
    data['OS_YEARS'] = data['OS_MONTHS'] / 12.0

    # Censor at 10 years
    over_10 = data['OS_YEARS'] > 10
    data.loc[over_10, 'OS_YEARS'] = 10
    data.loc[over_10, 'event'] = 0

    # Create Kaplan-Meier fitter
    kmf = KaplanMeierFitter()

    # Plot survival curves (square plot)
    plt.figure(figsize=(8, 8))

    # Fit and plot for each risk group
    for risk_group in ['Low', 'High']:
        group_data = data[data['risk_group'] == risk_group]
        kmf.fit(group_data['OS_YEARS'], 
                group_data['event'],
                label=f'{risk_group} Risk (n={len(group_data)})')
        kmf.plot(ci_show=True)

    # Perform log-rank test
    high_risk_data = data[data['risk_group'] == 'High']
    low_risk_data = data[data['risk_group'] == 'Low']

    results = logrank_test(high_risk_data['OS_YEARS'],
                          low_risk_data['OS_YEARS'],
                          high_risk_data['event'],
                          low_risk_data['event'])

    # Add labels and title
    plt.xlabel('Time (years)')
    plt.ylabel('Survival Probability')
    plt.title(f'Survival Curves by Oncotype DX Risk Group (10-Year)\nLog-rank p-value: {results.p_value:.3e}')
    plt.grid(True)
    plt.ylim(0, 1.0)

    # Print statistics
    print(f"\nNumber of high-risk patients: {len(high_risk_data)}")
    print(f"Number of low-risk patients: {len(low_risk_data)}")
    print(f"Log-rank test p-value: {results.p_value:.3e}")

    # Calculate median survival times (in years)
    high_risk_median = high_risk_data['OS_YEARS'].median()
    low_risk_median = low_risk_data['OS_YEARS'].median()
    print(f"Median survival time - High risk: {high_risk_median:.2f} years")
    print(f"Median survival time - Low risk: {low_risk_median:.2f} years")

    # Save plot
    plt.savefig(output_file)
    print(f"\nSurvival curves saved to {output_file}")

def main():
    # Load data
    print("Loading clinical data...")
    clinical_data = load_clinical_data()
    
    print("Loading Oncotype DX scores...")
    scores = load_oncotype_scores()
    
    # Filter clinical data for valid OS_MONTHS
    clinical_data = clinical_data[clinical_data['OS_MONTHS'].notna() & (clinical_data['OS_MONTHS'] > 0)]
    clinical_data['OS_MONTHS'] = pd.to_numeric(clinical_data['OS_MONTHS'], errors='coerce')
    clinical_data = clinical_data.dropna(subset=['OS_MONTHS'])
    
    print(f"\nNumber of samples with valid survival data: {len(clinical_data)}")
    
    # Plot survival curves
    plot_survival_curves(clinical_data, scores)

if __name__ == "__main__":
    main() 