import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle # Ensures pickle is used
from sklearn.metrics import roc_curve, auc
import os

def load_data(data_type='metabric'):
    """Load signature and segment scores"""
    print(f"Loading {data_type} data...")
    if data_type == 'metabric': # This is validation
        signature_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/validation_data/metabric_scores/rna_signature_score_cleaned.pkl"
        segment_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/validation_data/metabric_scores/cna_signature_score_cleaned.pkl"
    else:  # testing data
        signature_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/training_data/rna_signature_score_median_no_norm.pkl"
        segment_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/training_data/cna_segment_score_mean_no_norm.pkl"
    
    print(f"Loading signature scores from: {signature_path}")
    print(f"Loading segment scores from: {segment_path}")
    
    signature_score = pd.read_pickle(signature_path)
    segment_score = pd.read_pickle(segment_path)
    
    print(f"Signature score shape: {signature_score.shape}")
    print(f"Segment score shape: {segment_score.shape}")
    
    # Normalize sample IDs for testing data
    if data_type == 'testing':
        # Assuming signature_score uses hyphens, e.g., TCGA-XX-YYYY-ZZ
        # And segment_score might use dots, e.g., TCGA.XX.YYYY.ZZ
        # We will convert segment_score columns to use hyphens for consistency
        original_segment_cols = segment_score.columns.tolist()
        segment_score.columns = [col.replace('.', '-') for col in segment_score.columns]
        # Check if replacement actually happened for a few columns to debug
        if len(original_segment_cols) > 0 and len(segment_score.columns) > 0:
            print(f"Sample ID normalization for testing segment scores:")
            for i in range(min(3, len(original_segment_cols))):
                print(f"  Original: {original_segment_cols[i]} -> New: {segment_score.columns[i]}")
    
    # Find common samples
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples in {data_type}")
    
    if len(common_samples) == 0:
        print("WARNING: No common samples found!")
        print("Signature score columns:", signature_score.columns.tolist()[:5], "...")
        print("Segment score columns:", segment_score.columns.tolist()[:5], "...")
    
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
    
    # Define model paths
    model_paths_map = {
        'UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450': '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450_model.pkl',
        'GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP': '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP_model.pkl',
        'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN': '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN_model.pkl'
    }
    
    # Load both testing and validation data
    test_signature_score, test_segment_score = load_data('testing')
    metabric_signature_score, metabric_segment_score = load_data('metabric')
    
    # Plot ROC curves for each signature
    for signature in signatures:
        # Create figure for this signature with square aspect ratio
        plt.figure(figsize=(8, 8))
        
        # Load the model
        model_path = model_paths_map[signature]
        print(f"\nLoading model for {signature} from {model_path} using pickle.load with encoding='latin1'")
        try:
            with open(model_path, 'rb') as f:
                model_info = pickle.load(f, encoding='latin1')
        except Exception as e:
            print(f"Error loading model with pickle (encoding='latin1'): {e}")
            print("This model might have been saved with joblib or a different library, or might be corrupted.")
            continue # Skip to the next signature if a model fails to load
        
        print("Model info type:", type(model_info))
        
        # Get feature names and model from model_info
        if isinstance(model_info, dict):
            print("Model loaded as a dictionary. Keys:", model_info.keys())
            if 'feature_names' in model_info:
                feature_names = model_info['feature_names']
            else:
                print("WARNING: 'feature_names' not found in model_info dictionary. Using segment score index as fallback.")
                # Ensure test_segment_score is available and not empty
                if 'test_segment_score' in locals() and not test_segment_score.empty:
                    feature_names = test_segment_score.index.tolist()
                elif 'metabric_segment_score' in locals() and not metabric_segment_score.empty:
                     print("Using metabric_segment_score index for feature names as test_segment_score is not available.")
                     feature_names = metabric_segment_score.index.tolist()
                else:
                    print("ERROR: Cannot determine feature names. Both test and metabric segment scores are unavailable or empty.")
                    continue # Skip this signature
            
            if 'model' in model_info:
                model = model_info['model']
            else:
                print("WARNING: 'model' not found in model_info dictionary. Assuming model_info itself is the model.")
                model = model_info # This case might be problematic if feature_names were also missing
        else:
            print("Model loaded directly (not a dictionary). Assuming model_info is the model object.")
            model = model_info
            # Attempt to get feature names from common attributes if model is a scikit-learn pipeline or estimator
            if hasattr(model, 'feature_names_in_'):
                feature_names = model.feature_names_in_
            elif hasattr(model, 'feature_names_') : # Some older versions or custom estimators
                 feature_names = model.feature_names_
            elif 'test_segment_score' in locals() and not test_segment_score.empty:
                print("WARNING: Could not infer feature names from the model object. Using segment score index as fallback.")
                feature_names = test_segment_score.index.tolist()
            elif 'metabric_segment_score' in locals() and not metabric_segment_score.empty:
                print("WARNING: Could not infer feature names from the model object. Using metabric segment score index as fallback.")
                feature_names = metabric_segment_score.index.tolist()
            else:
                print("ERROR: Cannot determine feature names. Model is not a dict and segment scores are unavailable.")
                continue # Skip this signature

        # Clean feature names
        feature_names = [name.replace('\xa0', ' ') for name in feature_names]
        
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
        
        # Make predictions - adapt based on model structure
        y_pred_test = model.predict(X_test)
        y_pred_metabric = model.predict(X_metabric)
        
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
            title = 'RB LOH Signature'
            output_name = 'rb_loh'
        elif 'Basal_signaling' in signature:
            title = 'Basal Signaling Signature'
            output_name = 'basal_signaling'
        else:
            title = 'Estrogen Signaling Signature'
            output_name = 'estrogen_signaling'
        plt.title(f'ROC Curves for {title}')
        
        plt.legend(loc="lower right")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout(pad=2.0)  # Add more padding around the plot
        
        # Save the plot with new naming
        output_path = f'roc_curve_{output_name}.png'
        print(f"Saving plot to {output_path}")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    plot_roc_curves() 