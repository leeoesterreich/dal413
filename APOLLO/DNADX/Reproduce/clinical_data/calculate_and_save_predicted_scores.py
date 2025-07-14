import pandas as pd
import numpy as np
from joblib import load
import os
import pickle
from helper import GHI_RS
import gseapy as gp

def load_clinical_data():
    """Load clinical data and use OS_MONTHS directly"""
    clinical_data = pd.read_csv('data_clinical_patient.txt', sep='\t', skiprows=4)
    clinical_data.set_index('PATIENT_ID', inplace=True)
    return clinical_data

def load_metabric_data():
    """Load METABRIC CNA and mRNA data"""
    print("Loading METABRIC data...")
    cna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_cna_METABRIC.txt", sep='\t')
    mrna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_mrna_agilent_microarray.txt", sep='\t')
    return cna_data, mrna_data

def convert_symbol_to_entrez(symbols, gmt_file):
    """Convert gene symbols to Entrez IDs using the GMT file"""
    # Read GMT file
    with open(gmt_file, 'r') as f:
        gmt_lines = f.readlines()
    
    # Create symbol to EntrezID mapping
    symbol_to_entrez = {}
    for line in gmt_lines:
        parts = line.strip().split('\t')
        if len(parts) > 2:
            # The genes are in the third column onwards
            for gene in parts[2:]:
                if gene.isdigit():  # If it's an Entrez ID
                    symbol_to_entrez[gene] = gene  # Map the ID to itself
    
    print("\nGene mapping statistics:")
    print(f"Number of genes in mapping: {len(symbol_to_entrez)}")
    print("Sample of mapped genes:", list(symbol_to_entrez.items())[:5])
    
    # Convert symbols to EntrezIDs
    mapped_genes = [symbol_to_entrez.get(symbol, symbol) for symbol in symbols]
    print(f"\nNumber of input genes: {len(symbols)}")
    print(f"Number of mapped genes: {sum(1 for x in mapped_genes if x in symbol_to_entrez)}")
    print("Sample of input genes:", symbols[:5])
    print("Sample of mapped genes:", mapped_genes[:5])
    
    return mapped_genes

def predict_with_model(segment_scores):
    """Predict scores using the trained model"""
    # Load the model
    model_path = "GHI_RS_Model_NJEM.2004_PMID.15591335_model.joblib"
    print("\nLoading GHI_RS model for prediction")
    model_info = load(model_path)
    
    # Print model coefficients
    print("Model coefficients:", model_info['model'].coef_)
    print("Model intercept:", model_info['model'].intercept_)
    
    # Clean feature names
    feature_names = [name.replace('\xa0', ' ') for name in model_info['feature_names']]
    print("\nNumber of features in model:", len(feature_names))
    
    # Prepare data - use segment scores directly
    X = segment_scores.reindex(feature_names, fill_value=np.nan).values.T
    # Convert to float type
    X = X.astype(np.float64)
    print(f"X shape (after reindex and transpose): {X.shape}")
    print("\nFirst few rows of X:")
    print(X[:5, :5])
    
    # Print some statistics about X
    print("\nX statistics:")
    print("Number of NaN values:", np.isnan(X).sum())
    print("Number of non-NaN values:", (~np.isnan(X)).sum())
    
    # Make predictions manually, skipping NaN values
    predictions = []
    coef = model_info['model'].coef_
    intercept = model_info['model'].intercept_
    
    for i in range(X.shape[0]):
        sample = X[i]
        # Only use non-NaN values for prediction
        valid_mask = ~np.isnan(sample)
        valid_coef = coef[valid_mask]
        valid_sample = sample[valid_mask]
        pred = np.dot(valid_sample, valid_coef) + intercept
        predictions.append(pred)
    
    predictions = np.array(predictions)
    
    # Print prediction calculation for first sample
    first_sample = X[0]
    valid_mask = ~np.isnan(first_sample)
    valid_coef = coef[valid_mask]
    valid_sample = first_sample[valid_mask]
    manual_pred = np.dot(valid_sample, valid_coef) + intercept
    print("\nPrediction calculation for first sample:")
    print("Number of valid values used:", len(valid_sample))
    print("Valid values used:", valid_sample)
    print("Corresponding coefficients:", valid_coef)
    print("Dot product:", np.dot(valid_sample, valid_coef))
    print("Intercept:", intercept)
    print("Manual prediction:", manual_pred)
    print("Model prediction:", predictions[0])
    
    return predictions

def main():
    # Load clinical data
    clinical_data = load_clinical_data()
    
    # Load METABRIC data
    cna_data, mrna_data = load_metabric_data()
    
    # The sample IDs are columns in cna_data and mrna_data, starting from the third column
    cna_sample_ids = set(cna_data.columns[2:])
    mrna_sample_ids = set(mrna_data.columns[2:])
    clinical_sample_ids = set(clinical_data.index)
    common_samples = clinical_sample_ids & cna_sample_ids & mrna_sample_ids
    print(f"\nFound {len(common_samples)} common samples across all datasets")
    
    # Filter data to common samples
    clinical_data = clinical_data.loc[common_samples]
    cna_data = cna_data[['Hugo_Symbol', 'Entrez_Gene_Id'] + list(common_samples)]
    mrna_data = mrna_data[['Hugo_Symbol', 'Entrez_Gene_Id'] + list(common_samples)]
    
    # Filter out samples with non-zero, non-blank OS_MONTHS and ensure numeric values
    clinical_data = clinical_data[clinical_data['OS_MONTHS'].notna() & (clinical_data['OS_MONTHS'] > 0)]
    clinical_data['OS_MONTHS'] = pd.to_numeric(clinical_data['OS_MONTHS'], errors='coerce')
    clinical_data = clinical_data.dropna(subset=['OS_MONTHS'])
    print(f"\nFiltered to {len(clinical_data)} samples with usable OS_MONTHS.")
    
    # Set sample order to match clinical_data
    sample_order = list(clinical_data.index)
    
    # Load pre-calculated segment scores
    with open('/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/clinical_data/segment_score.pkl', 'rb') as f:
        segment_scores = pickle.load(f)
    
    # Filter segment scores to match our sample order
    segment_scores = segment_scores[sample_order]
    
    # Predict scores using the model
    model_scores = predict_with_model(segment_scores)
    print(f"Length of model_scores: {len(model_scores)}")
    print(f"Unique values in model_scores: {np.unique(model_scores)}")
    
    # Convert model_scores to a pandas Series with the appropriate index
    model_scores_series = pd.Series(model_scores, index=sample_order)
    
    # Save scores to CSV files
    model_scores_series.to_csv('predicted_signature_scores.csv')
    print("\nScores have been saved to predicted_signature_scores.csv.")

    # Calculate Oncotype DX score
    print("Calculating Oncotype DX score...")
    # Load mRNA data directly from the specified file
    mrna_data = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/data_mrna_agilent_microarray.txt", sep='\t')
    print("\nFirst few rows of raw mRNA data:")
    print(mrna_data.head())
    
    # Convert gene symbols to Entrez IDs using the GMT file
    gmt_file = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/DNA-based-predictors-of-non-genetic-cancer-phenotypes/data/gene_signatures_20170111.gmt"
    mrna_data['Entrez_Gene_Id'] = convert_symbol_to_entrez(mrna_data['Hugo_Symbol'], gmt_file)
    
    # Set Entrez IDs as index
    mrna_data.set_index('Entrez_Gene_Id', inplace=True)
    # Drop the Hugo_Symbol column
    mrna_data = mrna_data.drop('Hugo_Symbol', axis=1)
    
    # Print required gene IDs and check if they exist in the data
    required_genes = ['2597', '2990', '60', '7037', '6175', '2886', '2064', '596', '5241', '57758', '2099', 
                     '6790', '4605', '891', '332', '4288', '4320', '1515', '968', '2944', '573']
    print("\nRequired gene IDs:", required_genes)
    print("Available gene IDs in data:", mrna_data.index.tolist()[:10])
    missing_genes = [gene for gene in required_genes if gene not in mrna_data.index]
    print("Missing genes:", missing_genes)
    
    if missing_genes:
        print("\nWARNING: Some required genes are missing from the data!")
        print("This will affect the score calculation.")
    
    print(f"\nmRNA data shape: {mrna_data.shape}")
    print("First few rows of processed mRNA data:")
    print(mrna_data.head())
    
    # Check if we have any of the required genes
    found_genes = [gene for gene in required_genes if gene in mrna_data.index]
    print("\nFound genes in data:", found_genes)
    if found_genes:
        print("Expression values for found genes:")
        print(mrna_data.loc[found_genes].head())
    
    oncotype = GHI_RS(mrna_data)
    print(f"\nOncotype DX scores shape: {oncotype.shape}")
    print("First few Oncotype DX scores:")
    print(oncotype.head())
    
    # Create DataFrame with scores (samples as columns)
    signature_score = pd.DataFrame(oncotype, index=['GHI_RS_Model_NJEM.2004_PMID.15591335'])
    print("\nFinal signature score shape:", signature_score.shape)
    print("First few columns of signature score:")
    print(signature_score.iloc[:, :5])

    # Save signature scores
    print("\nSaving signature scores...")
    signature_score.to_csv("direct_signature_scores.csv")

if __name__ == "__main__":
    main() 