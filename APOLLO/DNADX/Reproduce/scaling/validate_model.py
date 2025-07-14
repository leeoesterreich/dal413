import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc, confusion_matrix
import matplotlib.pyplot as plt
import os
from joblib import load
from model_training import clean_feature_names

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
    total_nan = np.isnan(X).sum()
    print(f"Number of columns that are all NaN: {sum(nan_cols)}")
    print(f"Total NaN values in X: {total_nan}")
    print(f"Percentage of NaN values in X: {(total_nan / X.size) * 100:.2f}%")
    
    # Impute NaN values with column means
    col_means = np.nanmean(X, axis=0)
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(col_means, np.where(nan_mask)[1])
    
    # Make predictions
    print("\nMaking predictions...")
    y_pred = model_info['model'].predict(X)
    
    # Calculate correlation
    correlation = np.corrcoef(y, y_pred)[0, 1]
    print(f"Correlation between actual and predicted values: {correlation:.3f}")
    
    # Calculate MSE
    mse = np.mean((y - y_pred)**2)
    print(f"Mean squared error: {mse:.3f}")
    
    # Create safe signature name for file paths
    safe_signature = target_signature.replace('/', '_').replace('\\', '_')
    
    # Save actual vs predicted values
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
    plt.title(f'METABRIC Validation: Actual vs Predicted\nCorrelation: {correlation:.3f}')
    plot_png_path = f'results/validation_plots/{safe_signature}_actual_vs_predicted.png'
    plt.savefig(plot_png_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate threshold for upper 1/3 vs lower 2/3
    threshold = np.percentile(y, 66.67)
    y_bin = (y > threshold).astype(int)
    y_pred_bin = (y_pred > threshold).astype(int)
    
    # Calculate confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_bin, y_pred_bin).ravel()
    
    # Calculate metrics
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\nBinary classification stats:")
    print(f"Number of positive samples: {np.sum(y_bin)}")
    print(f"Number of negative samples: {len(y_bin) - np.sum(y_bin)}")
    print(f"True Positives: {tp}")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"Accuracy: {accuracy:.3f}")
    print(f"Precision: {precision:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1 Score: {f1:.3f}")
    
    # Calculate ROC curve and AUC
    fpr, tpr, _ = roc_curve(y_bin, y_pred)
    roc_auc = auc(fpr, tpr)
    
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
        'nan_percentage_before': (total_nan / X.size) * 100,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'true_positives': tp,
        'true_negatives': tn,
        'false_positives': fp,
        'false_negatives': fn
    }
    pd.DataFrame([metrics]).to_csv('results/validation_plots/metabric_metrics.csv', index=False)
    print("\nMetrics saved to results/validation_plots/metabric_metrics.csv")
    
    return metrics 