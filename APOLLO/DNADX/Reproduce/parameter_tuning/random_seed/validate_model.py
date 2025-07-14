import pandas as pd
import numpy as np
from joblib import load
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, roc_auc_score
import os

def clean_feature_names(feature_names):
    """Clean feature names by replacing non-breaking spaces with regular spaces"""
    return [name.replace('\xa0', ' ') for name in feature_names]

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

def load_metabric_data():
    """Load METABRIC signature and segment scores"""
    print("Loading METABRIC data...")
    signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/score/signature_score.pkl")
    segment_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/score/segment_score.pkl")
    
    # Find common samples
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples in METABRIC")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

def validate_model():
    """Validate the trained model on METABRIC dataset"""
    # Load the trained model
    model_path = "results/models/UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450_model.joblib"
    print(f"Loading model from {model_path}")
    model_info = load(model_path)
    
    # Load METABRIC data
    signature_score, segment_score = load_metabric_data()
    
    # Get the target signature
    target_signature = "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450"
    
    # Clean feature names in model_info
    model_info['feature_names'] = clean_feature_names(model_info['feature_names'])
    
    # Align METABRIC segment features to match training
    feature_names = model_info['feature_names']
    segment_score = segment_score.reindex(feature_names, fill_value=np.nan)
    print(f"Number of features in model_info['feature_names']: {len(feature_names)}")
    print(f"Number of features in METABRIC segment_score after reindexing: {segment_score.shape[0]}")
    missing_features = [f for f in feature_names if f not in segment_score.index]
    if missing_features:
        print(f"Missing features in METABRIC segment_score: {missing_features}")
    else:
        print("No missing features after reindexing.")
    
    # Only keep columns that are in feature_names (should be all, but robust)
    segment_score = segment_score.loc[feature_names]
    print(f"Number of features after selecting training-used features in validation: {len(segment_score)}")
    
    # Prepare data
    X = segment_score.values.T
    y = signature_score.loc[target_signature].values
    
    # Convert to numeric types
    X = X.astype(np.float64)
    y = y.astype(np.float64)
    
    print(f"METABRIC sample count: {len(X)}")
    # Debug: print number of columns that are all NaN before imputation
    nan_cols = np.isnan(X).all(axis=0)
    print(f"Number of columns that are all NaN: {np.sum(nan_cols)}")
    if np.sum(nan_cols) > 0:
        print(f"Indices of all-NaN columns: {np.where(nan_cols)[0]}")
    
    # Print NaN statistics before imputation
    total_nan = np.isnan(X).sum()
    print(f"Total number of NaN values before imputation: {total_nan}")
    print(f"Percentage of NaN values before imputation: {(total_nan / X.size) * 100:.2f}%")
    
    # Apply the same imputation strategy as training
    X = impute_preserve_all_columns(X)
    y = np.nan_to_num(y, nan=np.nanmean(y))  # Simple imputation for y
    
    # Print NaN statistics after imputation
    total_nan = np.isnan(X).sum()
    print(f"Total number of NaN values after imputation: {total_nan}")
    print(f"Percentage of NaN values after imputation: {(total_nan / X.size) * 100:.2f}%")
    
    # Make predictions
    y_pred = model_info['model'].predict(X)

    # Remove extreme outliers in y_pred (e.g., abs(y_pred) > 1000)
    mask = np.abs(y_pred) <= 1000
    y_pred = y_pred[mask]
    y = y[mask]
    X = X[mask]

    # Calculate threshold for upper 1/3 vs lower 2/3
    threshold = np.percentile(y, 66.67)
    y_bin = (y > threshold).astype(int)

    # Calculate ROC curve and AUC
    fpr, tpr, _ = roc_curve(y_bin, y_pred)
    roc_auc = auc(fpr, tpr)

    # Calculate correlation
    correlation = np.corrcoef(y, y_pred)[0, 1]

    # Create plots directory if it doesn't exist
    os.makedirs('results/validation_plots', exist_ok=True)

    # Plot ROC curve
    plt.figure(figsize=(10, 6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, 
             label=f'ROC curve (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve for METABRIC Validation\n(Upper 1/3 vs Lower 2/3)')
    plt.legend(loc="lower right")
    plt.savefig('results/validation_plots/metabric_roc_curve.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Plot actual vs predicted
    plt.figure(figsize=(10, 6))
    plt.scatter(y, y_pred, alpha=0.5)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'METABRIC Validation: Actual vs Predicted\nCorrelation: {correlation:.3f}')
    plt.savefig('results/validation_plots/metabric_actual_vs_predicted.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Print performance metrics
    print("\nMETABRIC Validation Performance:")
    print(f"Correlation: {correlation:.3f}")
    print(f"AUC: {roc_auc:.3f}")

    # Save metrics
    metrics = {
        'correlation': correlation,
        'auc': roc_auc,
        'threshold': threshold,
        'total_nan_before': total_nan,
        'nan_percentage_before': (total_nan / X.size) * 100
    }
    pd.DataFrame([metrics]).to_csv('results/validation_plots/metabric_metrics.csv', index=False)
    print("\nMetrics saved to results/validation_plots/metabric_metrics.csv")

if __name__ == "__main__":
    validate_model() 