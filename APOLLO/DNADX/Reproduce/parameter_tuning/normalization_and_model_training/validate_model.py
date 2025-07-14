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
    signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/validation_data/metabric_scores/rna_signature_score.pkl")
    segment_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/validation_data/metabric_scores/cna_signature_score.pkl")
    
    print("\nInitial data shapes:")
    print(f"signature_score shape: {signature_score.shape}")
    print(f"segment_score shape: {segment_score.shape}")
    
    print("\nNaN values in signature_score:")
    print(f"Total NaN values: {signature_score.isna().sum().sum()}")
    print(f"Percentage of NaN values: {(signature_score.isna().sum().sum() / signature_score.size) * 100:.2f}%")
    
    print("\nNaN values in segment_score:")
    print(f"Total NaN values: {segment_score.isna().sum().sum()}")
    print(f"Percentage of NaN values: {(segment_score.isna().sum().sum() / segment_score.size) * 100:.2f}%")
    
    # Find common samples
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples in METABRIC")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    print("\nNaN values after aligning common samples:")
    print(f"signature_score NaN values: {signature_score.isna().sum().sum()}")
    print(f"segment_score NaN values: {segment_score.isna().sum().sum()}")
    
    return signature_score, segment_score

def validate_model(model_path, target_signature):
    """Validate the trained model on METABRIC dataset"""
    # Load the trained model
    print(f"Loading model from {model_path}")
    model_info = load(model_path)
    
    # Load METABRIC data
    signature_score, segment_score = load_metabric_data()
    
    print("\nData shapes after loading:")
    print(f"signature_score shape: {signature_score.shape}")
    print(f"segment_score shape: {segment_score.shape}")
    print(f"Target signature '{target_signature}' exists: {target_signature in signature_score.index}")
    
    # Clean feature names in model_info
    model_info['feature_names'] = clean_feature_names(model_info['feature_names'])
    print(f"\nNumber of features in model: {len(model_info['feature_names'])}")
    print("First few feature names in model:", model_info['feature_names'][:5])
    
    # Align METABRIC segment features to match training
    feature_names = model_info['feature_names']
    segment_score = segment_score.reindex(feature_names, fill_value=np.nan)
    print(f"\nNumber of features in model_info['feature_names']: {len(feature_names)}")
    print(f"Number of features in METABRIC segment_score after reindexing: {segment_score.shape[0]}")
    # Only keep columns that are in feature_names (should be all, but robust)
    segment_score = segment_score.loc[feature_names]
    print(f"Number of features after selecting training-used features in validation: {len(segment_score)}")
    # Prepare data
    X = segment_score.values.T
    y = signature_score.loc[target_signature].values
    print("\nData shapes before processing:")
    print(f"X shape: {X.shape}")
    print(f"y shape: {y.shape}")
    # Convert to numeric types
    X = X.astype(np.float64)
    y = y.astype(np.float64)
    print(f"\nMETABRIC sample count: {len(X)}")
    # Debug: print number of columns that are all NaN before imputation
    nan_cols = np.isnan(X).all(axis=0)
    print(f"Number of columns that are all NaN: {np.sum(nan_cols)}")
    if np.sum(nan_cols) > 0:
        print(f"Indices of all-NaN columns: {np.where(nan_cols)[0]}")
        print("First few feature names corresponding to NaN columns:", [feature_names[i] for i in np.where(nan_cols)[0][:5]])
    # Print NaN statistics before imputation
    total_nan = np.isnan(X).sum()
    print(f"Total number of NaN values before imputation: {total_nan}")
    print(f"Percentage of NaN values before imputation: {(total_nan / X.size) * 100:.2f}%")
    # Apply the same imputation strategy as training
    X = impute_preserve_all_columns(X)
    y = np.nan_to_num(y, nan=np.nanmean(y))  # Simple imputation for y
    print("\nData shapes after imputation:")
    print(f"X shape: {X.shape}")
    print(f"y shape: {y.shape}")
    # Print NaN statistics after imputation
    total_nan = np.isnan(X).sum()
    print(f"Total number of NaN values after imputation: {total_nan}")
    print(f"Percentage of NaN values after imputation: {(total_nan / X.size) * 100:.2f}%")
    # Make predictions
    y_pred = model_info['model'].predict(X)
    print(f"\nPredictions shape: {y_pred.shape}")
    print(f"Predictions range: [{y_pred.min():.2f}, {y_pred.max():.2f}]")
    # Remove extreme outliers in y_pred (e.g., abs(y_pred) > 1000)
    mask = np.abs(y_pred) <= 1000
    y_pred = y_pred[mask]
    y = y[mask]
    X = X[mask]
    print(f"\nData shapes after outlier removal:")
    print(f"X shape: {X.shape}")
    print(f"y shape: {y.shape}")
    print(f"y_pred shape: {y_pred.shape}")
    
    # Filter out samples where actual value > 40
    filter_mask = (y <= 40)
    y = y[filter_mask]
    y_pred = y_pred[filter_mask]
    X = X[filter_mask]
    print(f"\nData shapes after filtering actual values > 40:")
    print(f"X shape: {X.shape}")
    print(f"y shape: {y.shape}")
    print(f"y_pred shape: {y_pred.shape}")
    
    # Save actual vs predicted values for plotting
    safe_signature = target_signature.replace('/', '_').replace('\\', '_')
    plot_df = pd.DataFrame({'actual_value': y, 'predicted_value': y_pred})
    plot_csv_path = f'results/validation_plots/{safe_signature}_actual_vs_predicted.csv'
    plot_df.to_csv(plot_csv_path, index=False)
    print(f"Saved actual vs predicted CSV for {target_signature} to {plot_csv_path}")
    # Plot actual vs predicted
    plt.figure(figsize=(10, 6))
    plt.scatter(y, y_pred, alpha=0.5)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'METABRIC Validation: Actual vs Predicted\nCorrelation: {np.corrcoef(y, y_pred)[0, 1]:.3f}')
    plot_png_path = f'results/validation_plots/{safe_signature}_actual_vs_predicted.png'
    plt.savefig(plot_png_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate threshold for upper 1/3 vs lower 2/3
    threshold = np.percentile(y, 66.67)
    y_bin = (y > threshold).astype(int)
    
    print(f"\nBinary classification stats:")
    print(f"Number of positive samples: {np.sum(y_bin)}")
    print(f"Number of negative samples: {len(y_bin) - np.sum(y_bin)}")
    
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
    return metrics

if __name__ == "__main__":
    validate_model()

# Directory containing model files
model_dir = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/results/models"

# List to store validation results
validation_results = []

# Iterate over each model file
for model_file in os.listdir(model_dir):
    if model_file.endswith('.pkl'):
        model_path = os.path.join(model_dir, model_file)
        model = load(model_path)
        
        # Extract signature name from the model file name
        signature_name = model_file.replace('_model.pkl', '')
        
        # Validate the model
        metrics = validate_model(model_path, signature_name)
        
        # Store the results
        validation_results.append({
            'signature': signature_name,
            **metrics
        })

# Convert results to a DataFrame for easier analysis
results_df = pd.DataFrame(validation_results)

# Save results to a CSV file
results_df.to_csv('results/validation_summary/validation_metrics.csv', index=False) 