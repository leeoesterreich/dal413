import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from joblib import load
from sklearn.metrics import roc_curve, auc
from model_training_rs100 import load_data, build_prediction_model
from validate_model import clean_feature_names, impute_preserve_all_columns

def load_all_signatures():
    """Load all signatures from the signature_score.pkl file"""
    print("Loading signature scores...")
    signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/signature_score.pkl")
    all_signatures = signature_score.index.tolist()
    print("Found {} signatures in total".format(len(all_signatures)))
    return pd.DataFrame({'signature': all_signatures})

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

def validate_model_with_outlier_removal(model_info, signature_score, segment_score, target_signature):
    """Validate model using the approach from validate_model.py with outlier removal"""
    print("\nValidating model for {}".format(target_signature))
    
    # Clean feature names in model_info
    model_info['feature_names'] = clean_feature_names(model_info['feature_names'])
    
    # Align METABRIC segment features to match training
    feature_names = model_info['feature_names']
    segment_score = segment_score.reindex(feature_names, fill_value=np.nan)
    
    # Only keep columns that are in feature_names
    segment_score = segment_score.loc[feature_names]
    
    # Prepare data
    X = segment_score.values.T
    y = signature_score.loc[target_signature].values
    
    # Convert to numeric types
    X = X.astype(np.float64)
    y = y.astype(np.float64)
    
    # Apply imputation
    X = impute_preserve_all_columns(X)
    y = np.nan_to_num(y, nan=np.nanmean(y))
    
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
    
    # Create plots directory
    os.makedirs('results_rs100/validation_plots', exist_ok=True)
    
    # Create safe filename
    safe_name = target_signature.replace('/', '_').replace('\\', '_')
    
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
    plt.savefig(f'results_rs100/validation_plots/{safe_name}_roc_curve.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot actual vs predicted
    plt.figure(figsize=(10, 6))
    plt.scatter(y, y_pred, alpha=0.5)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title(f'METABRIC Validation: Actual vs Predicted\nCorrelation: {correlation:.3f}')
    plt.savefig(f'results_rs100/validation_plots/{safe_name}_actual_vs_predicted.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return correlation, roc_auc

def train_and_validate_signatures():
    """Train and validate models for all signatures"""
    # Load all signatures
    all_signatures = load_all_signatures()
    
    # Load training data
    signature_score, segment_score = load_data()
    
    # Load METABRIC validation data
    metabric_signature_score, metabric_segment_score = load_metabric_data()
    
    # Create results directory
    os.makedirs('results_rs100/validation_summary', exist_ok=True)
    os.makedirs('results_rs100/models', exist_ok=True)
    os.makedirs('results_rs100/plots', exist_ok=True)
    os.makedirs('results_rs100/validation_plots', exist_ok=True)
    
    # Store results for visualization
    results = []
    
    # Process each signature
    for _, row in all_signatures.iterrows():
        signature = row['signature']
        print("\nProcessing signature: {}".format(signature))
        
        try:
            # Train model
            print("\nTraining model...")
            model, _, _, _ = build_prediction_model(signature_score, segment_score, signature)
            
            # Load the trained model
            safe_name = signature.replace('/', '_').replace('\\', '_')
            model_info = load('results_rs100/models/{}_model.joblib'.format(safe_name))
            
            # Validate model with outlier removal
            print("\nValidating model...")
            correlation, auc = validate_model_with_outlier_removal(
                model_info, metabric_signature_score, metabric_segment_score, signature
            )
            
            # Store results
            results.append({
                'signature': signature,
                'train_auc': row['train_auc'],
                'test_auc': row['test_auc'],
                'validation_auc': auc,
                'validation_correlation': correlation
            })
            
        except Exception as e:
            print("Error processing signature {}: {}".format(signature, str(e)))
            continue
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results
    results_df.to_csv('results_rs100/validation_summary/auc_comparison.csv', index=False)
    
    # Create visualization
    plot_auc_comparison(results_df)

def plot_auc_comparison(results_df):
    """Create visualization comparing train, test, and validation AUCs"""
    plt.figure(figsize=(15, 8))
    
    # Sort results by validation AUC
    results_df = results_df.sort_values('validation_auc', ascending=False)
    
    # Create bar plot
    x = np.arange(len(results_df))
    width = 0.25
    
    plt.bar(x - width, results_df['train_auc'], width, label='Train AUC')
    plt.bar(x, results_df['test_auc'], width, label='Test AUC')
    plt.bar(x + width, results_df['validation_auc'], width, label='Validation AUC')
    
    # Customize plot
    plt.xlabel('Signatures')
    plt.ylabel('AUC')
    plt.title('Comparison of Train, Test, and Validation AUCs\n(with outlier removal)')
    plt.xticks(x, results_df['signature'], rotation=45, ha='right')
    plt.legend()
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig('results_rs100/validation_summary/auc_comparison.png')
    plt.close()
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print("\nTrain AUC:")
    print(results_df['train_auc'].describe())
    print("\nTest AUC:")
    print(results_df['test_auc'].describe())
    print("\nValidation AUC:")
    print(results_df['validation_auc'].describe())
    print("\nValidation Correlation:")
    print(results_df['validation_correlation'].describe())

if __name__ == "__main__":
    train_and_validate_signatures() 