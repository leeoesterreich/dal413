import pandas as pd
import numpy as np
import os

def load_aligned_data(rna_path="scaling/training_data/rna_signature_score_median_no_norm.pkl", 
                      cna_path="scaling/training_data/cna_segment_score_mean_no_norm.pkl"):
    """Loads RNA signature scores and CNA segment scores, aligns them by common sample names."""
    print(f"Loading RNA data from: {rna_path}")
    print(f"Loading CNA data from: {cna_path}")
    
    try:
        signature_score_df = pd.read_pickle(rna_path)
        segment_score_df = pd.read_pickle(cna_path)
    except FileNotFoundError as e:
        print(f"Error loading data: {e}. Please ensure paths are correct.")
        return None, None

    print("\nOriginal sample name formats (first 5):")
    print(f"RNA signature score samples: {list(signature_score_df.columns[:5])}")
    print(f"CNA segment score samples: {list(segment_score_df.columns[:5])}")

    # Standardize CNA sample names (e.g., TCGA.A1.0001 -> TCGA-A1-0001)
    # Assuming RNA data uses '-' and CNA data might use '.'
    segment_score_df.columns = segment_score_df.columns.str.replace('.', '-', regex=False)
    
    common_samples = sorted(list(set(signature_score_df.columns) & set(segment_score_df.columns)))
    
    if not common_samples:
        print("No common samples found between RNA and CNA data after attempting to align names.")
        return None, None
        
    print(f"\nFound {len(common_samples)} common samples after name standardization.")
    print(f"First 5 common samples: {common_samples[:5]}")

    aligned_signature_score_df = signature_score_df[common_samples]
    aligned_segment_score_df = segment_score_df[common_samples]
    
    print(f"Aligned RNA signature scores shape: {aligned_signature_score_df.shape}")
    print(f"Aligned CNA segment scores shape: {aligned_segment_score_df.shape}")
    
    return aligned_signature_score_df, aligned_segment_score_df

def main():
    print("--- Calculating Standard Deviations for Input Data ---")
    
    rna_signatures_df, cna_features_df = load_aligned_data()

    if rna_signatures_df is None or cna_features_df is None:
        print("Exiting due to data loading or alignment issues.")
        return

    # Specific RNA signatures to analyze
    target_rna_signatures = [
        "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP",
        "GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN",
        "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450"
    ]

    print("\n--- Standard Deviations for Target RNA Signatures (sd(Y)) ---")
    for sig_name in target_rna_signatures:
        if sig_name in rna_signatures_df.index:
            y_values = rna_signatures_df.loc[sig_name].values
            y_std = np.std(y_values)
            print(f"Std dev for '{sig_name}': {y_std:.4f}")
        else:
            print(f"Warning: RNA signature '{sig_name}' not found in the data.")

    # CNA features (X data)
    # In cna_features_df, features are rows, samples are columns. Transpose for X.
    X_features_cna = cna_features_df.T 
    
    print("\n--- Standard Deviations for CNA Features (sd(X)) ---")
    if not X_features_cna.empty:
        # Calculate std dev for each feature (column in X_features_cna)
        # Ensure data is numeric and handle potential NaNs before std calculation
        feature_std_devs = X_features_cna.apply(lambda x: np.std(x.astype(float).dropna())).values
        
        if feature_std_devs.size > 0:
            print(f"Number of CNA features: {len(feature_std_devs)}")
            print(f"Mean of feature std devs: {np.mean(feature_std_devs):.4f}")
            print(f"Median of feature std devs: {np.median(feature_std_devs):.4f}")
            print(f"Min of feature std devs: {np.min(feature_std_devs):.4f}")
            print(f"Max of feature std devs: {np.max(feature_std_devs):.4f}")
            
            print("\nExample feature std devs (first 5):")
            for i, feature_name in enumerate(X_features_cna.columns[:5]):
                # Recalculate for the specific feature to ensure correct handling if needed for print
                feature_data = X_features_cna[feature_name].astype(float).dropna()
                if not feature_data.empty:
                     print(f"  Std dev for feature '{feature_name}': {np.std(feature_data):.4f}")
                else:
                     print(f"  Std dev for feature '{feature_name}': (No valid data after dropna)")
        else:
            print("No feature standard deviations could be calculated (feature_std_devs is empty).")
    else:
        print("CNA feature data (X_features_cna) is empty.")

if __name__ == "__main__":
    main() 