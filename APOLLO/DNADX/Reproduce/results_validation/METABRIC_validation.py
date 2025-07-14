import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from joblib import load
from sklearn.metrics import roc_curve, auc
from model_training import load_data, build_prediction_model
from validate_model import load_metabric_data, clean_feature_names, impute_preserve_all_columns

def load_all_signatures():
    """Load all signatures from the metrics file"""
    print("Loading model metrics...")
    metrics_df = pd.read_csv("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/plots/model_metrics.csv")
    print(f"\nFound {len(metrics_df)} signatures in total")
    return metrics_df

def validate_model_with_outlier_removal(model_info, signature_score, segment_score, target_signature):
    """Validate model using the approach from validate_model.py with outlier removal (actual or predicted > 40)"""
    print(f"\nValidating model for {target_signature}")
    # Clean feature names in model_info
    model_info['feature_names'] = clean_feature_names(model_info['feature_names'])
    # Align METABRIC segment features to match training
    feature_names = model_info['feature_names']
    segment_score = segment_score.reindex(feature_names, fill_value=np.nan)
    segment_score = segment_score.loc[feature_names]
    # Prepare data
    X = segment_score.values.T
    y = signature_score.loc[target_signature].values
    X = X.astype(np.float64)
    y = y.astype(np.float64)
    X = impute_preserve_all_columns(X)
    y = np.nan_to_num(y, nan=np.nanmean(y))
    y_pred = model_info['model'].predict(X)
    # Remove samples where actual or predicted > 40
    mask = (y <= 40) & (y_pred <= 40)
    y = y[mask]
    y_pred = y_pred[mask]
    X = X[mask]
    # Save actual vs predicted CSV
    safe_signature = target_signature.replace('/', '_').replace('\\', '_')
    plot_df = pd.DataFrame({'actual_value': y, 'predicted_value': y_pred})
    plot_csv_path = f'results/validation_plots/{safe_signature}_actual_vs_predicted.csv'
    plot_df.to_csv(plot_csv_path, index=False)
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
    # Calculate ROC curve and AUC
    fpr, tpr, _ = roc_curve(y_bin, y_pred)
    roc_auc = auc(fpr, tpr)
    # Calculate correlation
    correlation = np.corrcoef(y, y_pred)[0, 1]
    return correlation, roc_auc, len(y)

def train_and_validate_signatures():
    """Train and validate models for all signatures"""
    all_signatures = load_all_signatures()
    signature_score, segment_score = load_data()
    metabric_signature_score, metabric_segment_score = load_metabric_data()
    os.makedirs('results/validation_summary', exist_ok=True)
    os.makedirs('results/validation_plots', exist_ok=True)
    results = []
    for _, row in all_signatures.iterrows():
        signature = row['signature']
        print(f"\nProcessing signature: {signature}")
        try:
            model, _, _, _ = build_prediction_model(signature_score, segment_score, signature)
            safe_name = signature.replace('/', '_').replace('\\', '_')
            model_info = load(f'results/models/{safe_name}_model.joblib')
            correlation, auc, n_samples = validate_model_with_outlier_removal(
                model_info, metabric_signature_score, metabric_segment_score, signature
            )
            results.append({
                'signature': signature,
                'train_auc': row['train_auc'],
                'test_auc': row['test_auc'],
                'validation_auc': auc,
                'validation_correlation': correlation,
                'n_samples_after_filtering': n_samples
            })
        except Exception as e:
            print(f"Error processing signature {signature}: {str(e)}")
            continue
    results_df = pd.DataFrame(results)
    results_df.to_csv('results/validation_summary/validation_metrics.csv', index=False)
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
    plt.savefig('results/validation_summary/auc_comparison.png')
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