import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from joblib import load
import os

def load_clinical_data():
    """Load clinical data"""
    clinical_data = pd.read_csv('data_clinical_patient.txt', sep='\t', skiprows=4)
    clinical_data.set_index('PATIENT_ID', inplace=True)
    return clinical_data

def load_metabric_data():
    """Load METABRIC CNA and mRNA data"""
    print("Loading METABRIC data...")
    cna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_cna_METABRIC.txt", sep='\t')
    mrna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_mrna_agilent_microarray.txt", sep='\t')
    return cna_data, mrna_data

def calculate_direct_signature_score(mrna_data):
    """Calculate signature score directly using the provided code"""
    import sys
    sys.path.append('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/clinical_data')
    from signature_score_and_segment_score_calculation import GHI_RS
    
    # Calculate Oncotype DX score
    print("Calculating Oncotype DX score directly...")
    oncotype = GHI_RS(mrna_data)
    return oncotype

def predict_with_model(cna_data):
    """Predict scores using the trained model"""
    # Load the model
    model_path = "GHI_RS_Model_NJEM.2004_PMID.15591335_model.joblib"
    print("\nLoading GHI_RS model for prediction")
    model_info = load(model_path)
    
    # Clean feature names
    feature_names = [name.replace('\xa0', ' ') for name in model_info['feature_names']]
    
    # Prepare data
    X = cna_data.reindex(feature_names, fill_value=np.nan).values.T
    
    # Convert to numeric types and impute missing values
    X = X.astype(np.float64)
    X = impute_preserve_all_columns(X)
    
    # Make predictions
    predictions = model_info['model'].predict(X)
    return predictions

def impute_preserve_all_columns(X):
    """Impute each column: if all-NaN, fill with 0; else fill NaN with mean"""
    X_imputed = np.empty_like(X)
    for i in range(X.shape[1]):
        col = X[:, i]
        if np.all(np.isnan(col)):
            X_imputed[:, i] = 0
        else:
            mean_val = np.nanmean(col)
            col_imputed = np.where(np.isnan(col), mean_val, col)
            X_imputed[:, i] = col_imputed
    return X_imputed

def load_predicted_scores():
    """Load predicted signature scores"""
    print("Loading predicted scores...")
    scores = pd.read_csv('predicted_signature_scores.csv', index_col=0)
    scores.columns = ['score']  # Rename the column to 'score'
    return scores

def plot_survival_curves(clinical_data, scores, output_file='survival_curves.png'):
    """Plot survival curves for high and low risk groups (10-year analysis, square plot)"""
    # Merge clinical data with scores
    data = clinical_data.join(scores)

    # Calculate cutoff for top 1/3
    cutoff_score = data['score'].quantile(2/3)
    print(f"\nCutoff score (top 1/3): {cutoff_score}")

    # Assign risk groups
    data['risk_group'] = 'Low'
    data.loc[data['score'] > cutoff_score, 'risk_group'] = 'High'

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
    plt.title(f'Survival Curves by Risk Group (Predicted Scores, 10-Year)\nLog-rank p-value: {results.p_value:.3e}')
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
    
    print("Loading predicted scores...")
    scores = load_predicted_scores()
    
    # Filter clinical data for valid OS_MONTHS
    clinical_data = clinical_data[clinical_data['OS_MONTHS'].notna() & (clinical_data['OS_MONTHS'] > 0)]
    clinical_data['OS_MONTHS'] = pd.to_numeric(clinical_data['OS_MONTHS'], errors='coerce')
    clinical_data = clinical_data.dropna(subset=['OS_MONTHS'])
    
    print(f"\nNumber of samples with valid survival data: {len(clinical_data)}")
    
    # Plot survival curves
    plot_survival_curves(clinical_data, scores)

if __name__ == "__main__":
    main() 