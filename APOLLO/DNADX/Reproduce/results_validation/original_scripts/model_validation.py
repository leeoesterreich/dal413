import pandas as pd
import numpy as np
from joblib import load
from sklearn.metrics import roc_auc_score, mean_squared_error
import matplotlib.pyplot as plt
import os

def load_metabric_data():
    """Load the METABRIC validation dataset"""
    print("Loading METABRIC validation data...")
    # Load the signature scores and segment scores for METABRIC
    signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/results_validation/signature_score.pkl")
    segment_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/results_validation/segment_score.pkl")
    
    # Find common samples
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples in METABRIC dataset")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

def impute_preserve_all_columns(X):
    # Impute each column: if all-NaN, fill with 0; else fill NaN with mean
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

def validate_model(signature_score, segment_score, target_signature):
    """Validate the model on METABRIC dataset"""
    print(f"\nValidating model for {target_signature}")
    
    # Load the trained model and preprocessing information
    safe_name = target_signature.replace('/', '_').replace('\\', '_')
    model_info = load(f'results/models/{safe_name}_model.joblib')
    
    # Prepare validation data
    # Filter and order features to match training
    feature_names = model_info['feature_names']
    
    # Normalize whitespace in feature names and segment_score index
    def normalize_name(name):
        return ' '.join(str(name).replace('\xa0', ' ').split())

    normalized_feature_names = [normalize_name(f) for f in feature_names]
    normalized_index = [normalize_name(idx) for idx in segment_score.index]

    # Build a mapping from normalized name to original index
    index_map = {normalize_name(idx): idx for idx in segment_score.index}

    # Find which normalized feature names are missing from the normalized index
    missing_features = [f for f in normalized_feature_names if f not in index_map]
    if missing_features:
        print(f"ERROR: The following features are missing from segment_score after normalization: {missing_features}")
        raise ValueError("Missing features in validation data after normalization. Cannot proceed.")

    # Select features using the mapped original index
    selected_indices = [index_map[f] for f in normalized_feature_names]

    # Print type and value for each feature in selected_indices before building X_val_df
    print("\nChecking each feature in selected_indices:")
    for idx in selected_indices:
        val = segment_score.loc[idx]
        print(f"Feature: {repr(idx)} | type: {type(val)} | shape: {getattr(val, 'shape', None)}")
        if not hasattr(val, 'shape') or (hasattr(val, 'shape') and (1 in val.shape or 0 in val.shape)):
            print(f"  Value: {val}")

    # Select features using the mapped original index
    X_val_df = segment_score.loc[selected_indices]
    X_val = X_val_df.values.T
    y_val = signature_score.loc[target_signature].values
    
    # Print diagnostics for segment_score and X_val (now using normalized selection)
    print("\nDiagnostics for segment_score and X_val (after normalization):")
    print("selected_indices sample:", selected_indices[:5])
    print("X_val_df shape:", X_val_df.shape)
    print("X_val type:", type(X_val), "X_val shape:", X_val.shape)

    # Check for missing features in validation data
    missing_features = [f for f in feature_names if f not in segment_score.index]
    if missing_features:
        print(f"ERROR: The following features are missing from segment_score: {missing_features}")
        raise ValueError("Missing features in validation data. Cannot proceed.")

    # Debug: Check if target is accidentally in features
    print("\nChecking for target leakage in features...")
    for i, name in enumerate(feature_names):
        try:
            feature_values = X_val[:, i]
            print(f"Feature {name}: type={type(feature_values)}, shape={getattr(feature_values, 'shape', None)}")
            if not hasattr(feature_values, 'shape') or feature_values.shape[0] != len(y_val):
                print(f"Skipping feature {name} due to shape mismatch.")
                continue
            r = np.corrcoef(feature_values, y_val)[0, 1]
            if abs(r) > 0.99:
                print(f"WARNING: Feature {name} is nearly identical to target (corr={r:.3f})")
        except Exception as e:
            print(f"Error at feature index {i}, name {name}: {str(e)}")
            raise
    
    # Convert to numeric types
    X_val = X_val.astype(np.float64)
    y_val = y_val.astype(np.float64)
    
    # Apply the same preprocessing as training
    X_val = impute_preserve_all_columns(X_val)
    # Do not impute y_val, just use as is
    
    # Make predictions
    y_val_pred = model_info['model'].predict(X_val)
    
    # Debug: Print sample values and check for identity
    print("\nDebugging predictions:")
    print("Sample y_val:", y_val[:5])
    print("Sample y_val_pred:", y_val_pred[:5])
    print("Are y_val and y_val_pred identical?", np.allclose(y_val, y_val_pred))
    
    # Calculate metrics
    correlation = np.corrcoef(y_val, y_val_pred)[0, 1]
    mse = mean_squared_error(y_val, y_val_pred)
    
    # Calculate AUC using upper 1/3 vs lower 2/3
    threshold = np.percentile(y_val, model_info['threshold_percentile'])
    y_val_bin = (y_val > threshold).astype(int)
    auc = roc_auc_score(y_val_bin, y_val_pred)
    
    # Plot validation results
    plt.figure(figsize=(10, 5))
    
    # Plot 1: Actual vs Predicted
    plt.subplot(1, 2, 1)
    plt.scatter(y_val, y_val_pred, alpha=0.5)
    plt.plot([y_val.min(), y_val.max()], [y_val.min(), y_val.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'METABRIC Validation\nCorrelation: {correlation:.3f}, MSE: {mse:.3f}')
    
    # Plot 2: ROC Curve
    plt.subplot(1, 2, 2)
    from sklearn.metrics import roc_curve
    fpr, tpr, _ = roc_curve(y_val_bin, y_val_pred)
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {auc:.3f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve (Upper 1/3 vs Lower 2/3)')
    plt.legend(loc="lower right")
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs('results/validation_plots', exist_ok=True)
    plt.savefig(f'results/validation_plots/{safe_name}_metabric_validation.png')
    plt.close()
    
    # Save validation metrics
    validation_metrics = {
        'signature': target_signature,
        'correlation': correlation,
        'mse': mse,
        'auc': auc
    }
    
    metrics_file = 'results/validation_metrics.csv'
    metrics_df = pd.DataFrame([validation_metrics])
    
    if os.path.exists(metrics_file):
        existing_metrics = pd.read_csv(metrics_file)
        metrics_df = pd.concat([existing_metrics, metrics_df], ignore_index=True)
    
    metrics_df.to_csv(metrics_file, index=False)
    
    print(f"\nMETABRIC Validation Results:")
    print(f"Correlation: {correlation:.3f}")
    print(f"MSE: {mse:.3f}")
    print(f"AUC: {auc:.3f}")
    
    return correlation, mse, auc

def main():
    # Load METABRIC data
    signature_score, segment_score = load_metabric_data()
    
    # Process RB_LOH signature
    target_signature = "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450"
    print(f"\nProcessing RB_LOH signature")
    
    try:
        correlation, mse, auc = validate_model(signature_score, segment_score, target_signature)
        validation_results = [{
            'signature': target_signature,
            'correlation': correlation,
            'mse': mse,
            'auc': auc
        }]
        
        # Save validation results
        results_df = pd.DataFrame(validation_results)
        results_file = 'results/metabric_validation_results.csv'
        
        # Save results
        results_df.to_csv(results_file, index=False)
        print("\nValidation results saved to results/metabric_validation_results.csv")
        
        # Print summary
        print("\nMETABRIC Validation Performance:")
        print(results_df)
    except Exception as e:
        print(f"Error validating model: {str(e)}")

if __name__ == "__main__":
    main() 