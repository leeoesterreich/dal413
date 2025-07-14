import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from joblib import load
import os
from helper import GHI_RS, calc_segments
from datetime import datetime

def load_data():
    """Load all required data"""
    print("Loading data...")
    # Load clinical data
    clinical_data = pd.read_csv('data_clinical_patient.txt', sep='\t', skiprows=4)
    clinical_data.set_index('PATIENT_ID', inplace=True)
    
    # Load RNA data
    mrna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_mrna_agilent_microarray.txt", 
                           sep='\t', index_col=0)
    
    # Load CNA data
    cna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_cna_METABRIC.txt", 
                          sep='\t')
    
    return clinical_data, mrna_data, cna_data

def calculate_rna_based_score(mrna_data, clinical_data):
    """Calculate Oncotype DX score directly from RNA data and clean index"""
    print("\nCalculating Oncotype DX score from RNA data...")
    score = GHI_RS(mrna_data)
    score.name = 'RNA_Based_Score'
    # Remove any non-sample rows (e.g., 'Entrez_Gene_Id')
    valid_ids = set(clinical_data.index)
    score = score[score.index.isin(valid_ids)]
    return score

def predict_from_dna(cna_data, clinical_data):
    """Predict Oncotype DX score from DNA segment scores and clean index"""
    print("\nPredicting Oncotype DX score from DNA data...")
    # Calculate segment scores
    segment_score = calc_segments(cna_data, 
                                "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/DNA-based-predictors-of-non-genetic-cancer-phenotypes/data/CNA_segments.utf8.gmt", 
                                method='mean')
    
    # Load the model
    model_path = "GHI_RS_Model_NJEM.2004_PMID.15591335_model.joblib"
    print("Loading GHI_RS model for prediction")
    model_info = load(model_path)
    
    # Clean feature names
    feature_names = [name.replace('\xa0', ' ') for name in model_info['feature_names']]
    
    # Prepare data
    X = segment_score.reindex(feature_names, fill_value=np.nan).values.T
    
    # Convert to numeric types and impute missing values
    X = X.astype(np.float64)
    X = impute_preserve_all_columns(X)
    
    # Make predictions
    predictions = model_info['model'].predict(X)
    score = pd.Series(predictions, index=segment_score.columns)
    score.name = 'DNA_Based_Score'
    # Remove any non-sample rows
    valid_ids = set(clinical_data.index)
    score = score[score.index.isin(valid_ids)]
    return score

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

def plot_survival_curves(clinical_data, scores_dict, output_prefix='survival_curves'):
    """Plot survival curves for both RNA-based and DNA-based predictions"""
    # Add timestamp for unique filenames
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    # Filter clinical data for valid OS_MONTHS
    clinical_data = clinical_data[clinical_data['OS_MONTHS'].notna() & (clinical_data['OS_MONTHS'] > 0)]
    clinical_data['OS_MONTHS'] = pd.to_numeric(clinical_data['OS_MONTHS'], errors='coerce')
    clinical_data = clinical_data.dropna(subset=['OS_MONTHS'])
    
    print(f"\nNumber of samples with valid survival data: {len(clinical_data)}")
    
    # Create plots for each score type
    for score_type, scores in scores_dict.items():
        # Merge clinical data with scores
        data = clinical_data.join(scores)
        
        # Calculate cutoff for top 1/3
        cutoff_score = data[score_type].quantile(2/3)
        print(f"\n{score_type} cutoff score (top 1/3): {cutoff_score}")
        
        # Assign risk groups
        data['risk_group'] = 'Low'
        data.loc[data[score_type] > cutoff_score, 'risk_group'] = 'High'
        
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
        
        # Format p-value
        p_value = results.p_value
        if p_value < 0.0001:
            p_value_str = "<0.0001"
        else:
            p_value_str = f"{p_value:.3e}"
        
        # Add labels and title
        plt.xlabel('Time (years)')
        plt.ylabel('Survival Probability')
        plt.title(f'Survival Curves by {score_type} Risk Group (10-Year)\nLog-rank p-value: {p_value_str}')
        plt.grid(True)
        plt.ylim(0, 1.0)
        
        # Print statistics
        print(f"\n{score_type} Analysis:")
        print(f"Number of high-risk patients: {len(high_risk_data)}")
        print(f"Number of low-risk patients: {len(low_risk_data)}")
        print(f"Log-rank test p-value: {p_value_str}")
        
        # Calculate median survival times (in years)
        high_risk_median = high_risk_data['OS_YEARS'].median()
        low_risk_median = low_risk_data['OS_YEARS'].median()
        print(f"Median survival time - High risk: {high_risk_median:.2f} years")
        print(f"Median survival time - Low risk: {low_risk_median:.2f} years")
        
        # Save plot
        output_file = f"{output_prefix}_{score_type}_{timestamp}.png"
        plt.savefig(output_file)
        print(f"\nSurvival curves saved to {output_file}")

def main():
    # Load all data
    clinical_data, mrna_data, cna_data = load_data()
    
    # Calculate scores
    rna_based_score = calculate_rna_based_score(mrna_data, clinical_data)
    dna_based_score = predict_from_dna(cna_data, clinical_data)
    
    # Create scores dictionary
    scores_dict = {
        'RNA_Based_Score': rna_based_score,
        'DNA_Based_Score': dna_based_score
    }
    
    # Plot survival curves
    plot_survival_curves(clinical_data, scores_dict)

if __name__ == "__main__":
    main() 