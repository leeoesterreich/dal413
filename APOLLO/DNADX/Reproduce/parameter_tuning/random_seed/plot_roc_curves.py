import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from joblib import load
from sklearn.metrics import roc_curve, auc
import os

def load_data(data_type='metabric'):
    """Load signature and segment scores"""
    print(f"Loading {data_type} data...")
    if data_type == 'metabric':
        signature_score = pd.read_pickle("score/signature_score.pkl")
        segment_score = pd.read_pickle("score/segment_score.pkl")
    else:  # testing data
        signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/signature_score.pkl")
        segment_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/segment_score.pkl")
    
    # Find common samples
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples in {data_type}")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

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

def plot_roc_curves():
    # Define the signatures to plot
    signatures = [
        'UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450',
        'GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP',
        'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN'
    ]
    
    # Load both testing and validation data
    test_signature_score, test_segment_score = load_data('testing')
    metabric_signature_score, metabric_segment_score = load_data('metabric')
    
    # Plot ROC curves for each signature
    for signature in signatures:
        # Create figure for this signature with square aspect ratio
        plt.figure(figsize=(8, 8))
        
        # Load the model
        model_path = f"results/models/{signature}_model.joblib"
        print(f"\nLoading model for {signature}")
        model_info = load(model_path)
        
        # Clean feature names
        feature_names = [name.replace('\xa0', ' ') for name in model_info['feature_names']]
        
        # Prepare test data
        test_segment_score_aligned = test_segment_score.reindex(feature_names, fill_value=np.nan)
        X_test = test_segment_score_aligned.values.T
        y_test = test_signature_score.loc[signature].values
        
        # Prepare METABRIC data
        metabric_segment_score_aligned = metabric_segment_score.reindex(feature_names, fill_value=np.nan)
        X_metabric = metabric_segment_score_aligned.values.T
        y_metabric = metabric_signature_score.loc[signature].values
        
        # Convert to numeric types
        X_test = X_test.astype(np.float64)
        y_test = y_test.astype(np.float64)
        X_metabric = X_metabric.astype(np.float64)
        y_metabric = y_metabric.astype(np.float64)
        
        # Impute missing values
        X_test = impute_preserve_all_columns(X_test)
        X_metabric = impute_preserve_all_columns(X_metabric)
        y_test = np.nan_to_num(y_test, nan=np.nanmean(y_test))
        y_metabric = np.nan_to_num(y_metabric, nan=np.nanmean(y_metabric))
        
        # Make predictions
        y_pred_test = model_info['model'].predict(X_test)
        y_pred_metabric = model_info['model'].predict(X_metabric)
        
        # Calculate thresholds for binary classification
        threshold_test = np.percentile(y_test, 66.67)
        threshold_metabric = np.percentile(y_metabric, 66.67)
        
        # Convert to binary labels
        y_test_bin = (y_test > threshold_test).astype(int)
        y_metabric_bin = (y_metabric > threshold_metabric).astype(int)
        
        # Calculate ROC curves
        fpr_test, tpr_test, _ = roc_curve(y_test_bin, y_pred_test)
        fpr_metabric, tpr_metabric, _ = roc_curve(y_metabric_bin, y_pred_metabric)
        
        # Calculate AUC scores
        auc_test = auc(fpr_test, tpr_test)
        auc_metabric = auc(fpr_metabric, tpr_metabric)
        
        # Plot ROC curves
        plt.plot(fpr_test, tpr_test, color='red', lw=2,
                label=f'TCGA Testing (AUC = {auc_test:.3f})')
        plt.plot(fpr_metabric, tpr_metabric, color='blue', lw=2,
                label=f'METABRIC Validation (AUC = {auc_metabric:.3f})')
        
        # Customize plot
        plt.plot([0, 1], [0, 1], color='gray', linestyle='--', alpha=0.5)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        
        # Create a shorter title from the signature name
        if 'UNC_RB_LOH' in signature:
            title = 'UNC_RB_LOH'
        elif 'Basal_signaling' in signature:
            title = 'Basal signaling'
        else:
            title = 'Estrogen signaling'
        plt.title(f'ROC Curves for {title}')
        
        plt.legend(loc="lower right")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout(pad=2.0)  # Add more padding around the plot
        
        # Save the plot
        plt.savefig(f'roc_curves_{title.lower().replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    plot_roc_curves() 