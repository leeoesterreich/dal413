# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from joblib import load
from sklearn.metrics import roc_curve, auc, mean_squared_error
from model_training import load_data, clean_data, plot_model_performance
from validate_model import load_metabric_data, clean_feature_names



def load_all_signatures():
    """Load all signatures from the metrics file"""
    print("Loading model metrics...")
    metrics_df = pd.read_csv("results/all_model_metrics.csv")
    print("\nFound {} signatures in total".format(len(metrics_df)))
    return metrics_df

def load_metabric_scores():
    """Load METABRIC signature and segment scores"""
    validation_dir = "validation_data/metabric_scores"
    rna_path = "validation_data/metabric_scores/rna_signature_score.pkl"
    cna_path = "validation_data/metabric_scores/cna_signature_score.pkl"
    
    print("Loading scores from: {}".format(validation_dir))
    
    # Load RNA signature scores
    print("Loading RNA scores from: {}".format(rna_path))
    rna_scores = pd.read_pickle(rna_path)
    
    # Load CNA segment scores
    print("Loading CNA scores from: {}".format(cna_path))
    cna_scores = pd.read_pickle(cna_path)
    
    print("\nRNA Scores Info:")
    print("Shape: {}".format(rna_scores.shape))
    print("NaN values: {} ({:.2f}%)".format(rna_scores.isna().sum().sum(), (rna_scores.isna().sum().sum() / rna_scores.size) * 100))
    print("Value range: [{:.3f}, {:.3f}]".format(rna_scores.min().min(), rna_scores.max().max()))
    
    print("\nCNA Scores Info:")
    print("Shape: {}".format(cna_scores.shape))
    print("NaN values: {} ({:.2f}%)".format(cna_scores.isna().sum().sum(), (cna_scores.isna().sum().sum() / cna_scores.size) * 100))
    print("Value range: [{:.3f}, {:.3f}]".format(cna_scores.min().min(), cna_scores.max().max()))
    
    return rna_scores, cna_scores

def validate_model(model, X, y, signature_name):
    """Validate model on METABRIC data"""
    # Make predictions
    y_pred = model.predict(X)
    
    # Handle NaN and infinity values
    mask = np.isfinite(y.flatten()) & np.isfinite(y_pred.flatten())
    if not np.any(mask):
        raise ValueError("No valid (finite) values found in predictions or actual values")
    
    y_valid = y.flatten()[mask]
    y_pred_valid = y_pred.flatten()[mask]
    
    # Remove points where actual value > 40 or predicted value > 40
    filter_mask = (y_valid <= 40)
    y_valid = y_valid[filter_mask]
    y_pred_valid = y_pred_valid[filter_mask]
    
    if len(y_valid) < 10:  # Require at least 10 samples after filtering
        raise ValueError("Not enough valid samples after filtering (need at least 10, got {})".format(len(y_valid)))
    
    # Save actual and predicted values used for plotting
    safe_signature_name = signature_name.replace('/', '_').replace('\\', '_')
    plot_df = pd.DataFrame({
        'actual_value': y_valid,
        'predicted_value': y_pred_valid
    })
    plot_csv_path = "results/validation_plots/{}_actual_vs_predicted.csv".format(safe_signature_name)
    plot_df.to_csv(plot_csv_path, index=False)
    print("Saved actual vs predicted CSV for {} to {}".format(signature_name, plot_csv_path))
    
    # Find the sample causing high correlation (on filtered data only)
    correlations = []
    for i in range(len(y_valid)):
        y_temp = np.delete(y_valid, i)
        y_pred_temp = np.delete(y_pred_valid, i)
        corr = np.corrcoef(y_temp, y_pred_temp)[0, 1]
        correlations.append(corr)
    base_corr = np.corrcoef(y_valid, y_pred_valid)[0, 1]
    corr_changes = [abs(base_corr - corr) for corr in correlations]
    outlier_idx = np.argmax(corr_changes)
    print("\nOutlier Analysis for {} (after filtering):".format(signature_name))
    print("Base correlation (with all filtered samples): {:.3f}".format(base_corr))
    print("Correlation without sample {}: {:.3f}".format(outlier_idx, correlations[outlier_idx]))
    print("Change in correlation: {:.3f}".format(corr_changes[outlier_idx]))
    print("Outlier sample values - Actual: {:.3f}, Predicted: {:.3f}".format(y_valid[outlier_idx], y_pred_valid[outlier_idx]))
    
    # Remove the outlier sample (optional, as before)
    y_valid_wo = np.delete(y_valid, outlier_idx)
    y_pred_valid_wo = np.delete(y_pred_valid, outlier_idx)
    
    # Calculate metrics (on filtered data only)
    correlation = np.corrcoef(y_valid, y_pred_valid)[0, 1]
    mse = mean_squared_error(y_valid, y_pred_valid)
    fpr, tpr, _ = roc_curve(y_valid, y_pred_valid)
    roc_auc = auc(fpr, tpr)
    n_samples = len(y_valid)
    top_threshold = np.percentile(y_valid, 66.67)
    pred_threshold = np.percentile(y_pred_valid, 66.67)
    y_binary = (y_valid >= top_threshold).astype(int)
    y_pred_binary = (y_pred_valid >= pred_threshold).astype(int)
    tp = np.sum((y_binary == 1) & (y_pred_binary == 1))
    tn = np.sum((y_binary == 0) & (y_pred_binary == 0))
    fp = np.sum((y_binary == 0) & (y_pred_binary == 1))
    fn = np.sum((y_binary == 1) & (y_pred_binary == 0))
    accuracy = (tp + tn) / n_samples
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    os.makedirs('results/validation_plots', exist_ok=True)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    # Plot actual vs predicted (filtered data only)
    ax1.scatter(y_valid, y_pred_valid, alpha=0.5)
    ax1.axhline(y=pred_threshold, color='r', linestyle='--', label='Prediction threshold')
    ax1.axvline(x=top_threshold, color='g', linestyle='--', label='Actual threshold')
    ax1.plot([y_valid.min(), y_valid.max()], [y_valid.min(), y_valid.max()], 'k--', lw=2)
    ax1.set_xlabel('Actual')
    ax1.set_ylabel('Predicted')
    ax1.set_title('Validation: Actual vs Predicted\nCorrelation: {:.3f}, MSE: {:.3f}'.format(correlation, mse))
    ax1.legend()
    # Plot ROC curve
    ax2.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (AUC = {:.3f})'.format(roc_auc))
    ax2.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('False Positive Rate')
    ax2.set_ylabel('True Positive Rate')
    ax2.set_title('Validation: Receiver Operating Characteristic')
    ax2.legend(loc="lower right")
    # Plot distribution of actual values
    ax3.hist(y_valid, bins=30, alpha=0.5, label='Actual')
    ax3.axvline(x=top_threshold, color='r', linestyle='--', label='Top 1/3 threshold')
    ax3.set_xlabel('Value')
    ax3.set_ylabel('Count')
    ax3.set_title('Distribution of Actual Values')
    ax3.legend()
    # Plot distribution of predicted values
    ax4.hist(y_pred_valid, bins=30, alpha=0.5, label='Predicted')
    ax4.axvline(x=pred_threshold, color='r', linestyle='--', label='Top 1/3 threshold')
    ax4.set_xlabel('Value')
    ax4.set_ylabel('Count')
    ax4.set_title('Distribution of Predicted Values')
    ax4.legend()
    plt.tight_layout()
    plt.savefig('results/validation_plots/{}_validation_performance.png'.format(safe_signature_name))
    plt.close()
    return {
        'correlation': correlation,
        'mse': mse,
        'roc_auc': roc_auc,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'true_positives': tp,
        'true_negatives': tn,
        'false_positives': fp,
        'false_negatives': fn,
        'outlier_sample_idx': outlier_idx,
        'outlier_actual_value': y_valid_wo[outlier_idx] if len(y_valid_wo) > outlier_idx else None,
        'outlier_predicted_value': y_pred_valid_wo[outlier_idx] if len(y_pred_valid_wo) > outlier_idx else None,
        'correlation_with_outlier': base_corr,
        'correlation_without_outlier': correlations[outlier_idx],
        'n_samples_after_filtering': len(y_valid)
    }

def validate_all_models():
    """Validate all trained models on METABRIC dataset"""
    # Load METABRIC scores
    print("Loading METABRIC scores...")
    metabric_rna_scores, metabric_cna_scores = load_metabric_scores()
    
    # Print data information
    print("\nMETABRIC Data Information:")
    print("="*50)
    print("RNA scores shape: {}".format(metabric_rna_scores.shape))
    print("CNA scores shape: {}".format(metabric_cna_scores.shape))
    print("Number of samples: {}".format(metabric_rna_scores.shape[1]))
    
    # Get all model files
    import glob
    model_files = glob.glob('results/models/*.pkl')
    print("\nFound {} model files to validate".format(len(model_files)))
    
    results = []
    successful_validations = 0
    failed_validations = 0
    
    for i, model_file in enumerate(model_files):
        # Extract signature name from filename
        filename = os.path.basename(model_file)
        signature = filename.replace('_model.pkl', '')
        
        print("\n" + "="*80)
        print("Validating model {}/{}: {}".format(i+1, len(model_files), signature))
        print("="*80)
        
        try:
            # Load the trained model
            if not os.path.exists(model_file):
                print("Warning: Model file not found: {}".format(model_file))
                failed_validations += 1
                continue
                
            model_info = load(model_file)
            model = model_info['model']
            
            # Check if signature exists in RNA scores
            if signature not in metabric_rna_scores.index:
                print("Warning: Signature {} not found in RNA scores".format(signature))
                failed_validations += 1
                continue
                
            # Prepare METABRIC data
            # Reindex CNA data to match model's expected features
            feature_names = model_info['feature_names']
            cna_df = pd.DataFrame(metabric_cna_scores, columns=metabric_cna_scores.index)
            cna_df = cna_df.reindex(columns=feature_names, fill_value=np.nan)
            X = cna_df.values
            
            # Calculate mean of X for imputation
            X_means = np.nanmean(X, axis=0)
            # Handle NaN values in X by imputing with column means
            X = np.where(np.isnan(X), X_means, X)
            y = metabric_rna_scores.loc[signature].values.reshape(-1, 1)  # RNA score as target
            
            # Print data shapes and check for NaN values
            print("Data shapes - X: {}, y: {}".format(X.shape, y.shape))
            print("NaN values in X: {}".format(np.isnan(X).sum()))
            print("NaN values in y: {}".format(np.isnan(y).sum()))
            
            # Remove samples with NaN values in y
            mask = ~np.isnan(y.flatten())
            X = X[mask]
            y = y[mask]
            
            if len(X) < 10:  # Require at least 10 samples
                print("Warning: Not enough valid samples for {} (need at least 10, got {})".format(signature, len(X)))
                failed_validations += 1
                continue
                
            # Validate model
            metrics = validate_model(model, X, y, signature)
            
            # Store results
            result = {
                'signature': signature,
                'correlation': metrics['correlation'],
                'roc_auc': metrics['roc_auc'],
                'mse': metrics['mse'],
                'accuracy': metrics['accuracy'],
                'precision': metrics['precision'],
                'recall': metrics['recall'],
                'f1_score': metrics['f1_score'],
                'true_positives': metrics['true_positives'],
                'true_negatives': metrics['true_negatives'],
                'false_positives': metrics['false_positives'],
                'false_negatives': metrics['false_negatives'],
                'n_samples_after_filtering': metrics['n_samples_after_filtering']
            }
            results.append(result)
            successful_validations += 1
            
            print("Validation successful! Correlation: {:.3f}, AUC: {:.3f}".format(
                metrics['correlation'], metrics['roc_auc']))
                
        except Exception as e:
            print("Error validating signature {}: {}".format(signature, str(e)))
            failed_validations += 1
            continue
    
    # Create results directory if it doesn't exist
    os.makedirs('results/validation_summary', exist_ok=True)
    
    # Save results
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv('results/validation_summary/validation_metrics.csv', index=False)
        print("\n" + "="*80)
        print("VALIDATION SUMMARY")
        print("="*80)
        print("Total models: {}".format(len(model_files)))
        print("Successful validations: {}".format(successful_validations))
        print("Failed validations: {}".format(failed_validations))
        print("Success rate: {:.1f}%".format(100 * successful_validations / len(model_files)))
        print("\nTop 10 models by AUC:")
        top_models = results_df.nlargest(10, 'roc_auc')[['signature', 'correlation', 'roc_auc']]
        for idx, row in top_models.iterrows():
            print("  {}: Corr={:.3f}, AUC={:.3f}".format(row['signature'][:60], row['correlation'], row['roc_auc']))
        print("\nValidation complete! Results saved to results/validation_summary/validation_metrics.csv")
        print("New results have replaced the old validation data.")
    else:
        print("\nNo models were successfully validated")

if __name__ == "__main__":
    validate_all_models() 