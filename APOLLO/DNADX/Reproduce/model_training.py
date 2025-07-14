import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import ElasticNet
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, mean_squared_error
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from joblib import Parallel, delayed
import warnings
import os
from sklearn.impute import SimpleImputer
warnings.filterwarnings('ignore')

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("results/signature_score.pkl")
    segment_score = pd.read_pickle("results/segment_score.pkl")
    
    # Find common samples using exact names
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

def get_lambda_range(X, y, alpha):
    """Determine min and max lambda values for a given alpha using glmnet-like (ElasticNet) approach"""
    # Fit ElasticNet to get lambda range
    model = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=10000)
    model.fit(X, y)
    lambda_max = np.max(np.abs(model.coef_))
    lambda_min = lambda_max * 0.01
    return np.logspace(np.log10(lambda_min), np.log10(lambda_max), 100)  # 100 lambda values

def monte_carlo_cv(X, y, alpha, lambda_val, n_splits=200, train_size=0.7):  # Changed to 0.7 train size (30% test)
    """Perform Monte Carlo cross validation with 200 splits"""
    accuracies = []
    for _ in range(n_splits):
        # Split data with consistent random state
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, train_size=train_size, random_state=42  # Fixed random state for consistency
        )
        
        # Train model
        model = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=10000)
        model.fit(X_train, y_train)
        
        # Make predictions and calculate accuracy
        y_pred = model.predict(X_test)
        threshold = np.percentile(y_test, 66.67)  # Upper 1/3 threshold
        y_test_bin = (y_test > threshold).astype(int)
        y_pred_bin = (y_pred > threshold).astype(int)
        accuracy = accuracy_score(y_test_bin, y_pred_bin)
        accuracies.append(accuracy)
    
    return np.mean(accuracies)

def evaluate_params(args):
    """Helper function for parallel parameter evaluation"""
    X, y, alpha, lambda_val = args
    return (alpha, lambda_val, monte_carlo_cv(X, y, alpha, lambda_val))

def plot_model_performance(y_true, y_pred, signature_name, correlation, mse):
    """Plot model performance metrics"""
    # Create plots directory if it doesn't exist
    os.makedirs('results/plots', exist_ok=True)
    
    # Calculate ROC curve using upper 1/3 vs lower 2/3
    threshold = np.percentile(y_true, 66.67)  # Upper 1/3 threshold
    y_true_bin = (y_true > threshold).astype(int)
    
    fpr, tpr, _ = roc_curve(y_true_bin, y_pred)
    roc_auc = auc(fpr, tpr)
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Plot 1: Actual vs Predicted
    ax1.scatter(y_true, y_pred, alpha=0.5)
    ax1.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', lw=2)
    ax1.set_xlabel('Actual Values')
    ax1.set_ylabel('Predicted Values')
    ax1.set_title(f'Actual vs Predicted\nCorrelation: {correlation:.3f}, MSE: {mse:.3f}')
    
    # Plot 2: ROC Curve
    ax2.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
    ax2.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('False Positive Rate')
    ax2.set_ylabel('True Positive Rate')
    ax2.set_title('ROC Curve (Upper 1/3 vs Lower 2/3)')
    ax2.legend(loc="lower right")
    
    plt.tight_layout()
    
    # Save plot
    safe_name = signature_name.replace('/', '_').replace('\\', '_')
    plt.savefig(f'results/plots/{safe_name}_performance.png')
    plt.close()
    
    return roc_auc

def build_prediction_model(signature_score, segment_score, target_signature):
    """Build prediction model for a signature using Elastic Net"""
    print(f"\nBuilding prediction model for {target_signature}")
    
    # Skip problematic signature
    if target_signature == "UNC_Differentiation.Score_Model_BCR.2010_PMID.20813035":
        print("Skipping problematic signature")
        return None, None, None, None
    
    # Prepare data
    X = segment_score.values.T  # segments as features
    y = signature_score.loc[target_signature].values
    
    # Convert to numeric types
    X = X.astype(np.float64)
    y = y.astype(np.float64)
    
    print(f"Initial sample count: {len(X)}")
    
    # Handle missing values by imputation
    imputer = SimpleImputer(strategy='mean')
    X = imputer.fit_transform(X)
    y = imputer.fit_transform(y.reshape(-1, 1)).ravel()
    
    print(f"Sample count after imputation: {len(X)}")
    
    # Split data with 30% test size
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    # Parameter tuning
    print("Performing parameter tuning...")
    alphas = np.arange(0.1, 1.1, 0.1)  # 0.1 to 1.0 by 0.1
    
    # Prepare parameter combinations for parallel processing
    param_combinations = []
    for alpha in alphas:
        lambda_range = get_lambda_range(X_train, y_train, alpha)
        for lambda_val in lambda_range:
            param_combinations.append((X_train, y_train, alpha, lambda_val))
    
    # Run parameter tuning in parallel
    results = Parallel(n_jobs=-1)(delayed(evaluate_params)(params) for params in param_combinations)
    
    # Find best parameters
    best_accuracy = -np.inf
    best_params = None
    for alpha, lambda_val, accuracy in results:
        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_params = (alpha, lambda_val)
    
    print(f"Best parameters - alpha: {best_params[0]:.3f}, lambda: {best_params[1]:.3f}")
    
    # Train final model with best parameters
    model = ElasticNet(alpha=best_params[0], l1_ratio=0.5, max_iter=10000)
    model.fit(X_train, y_train)
    
    # Get predictions for both train and test sets
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)
    
    # Calculate metrics for both sets
    train_correlation = np.corrcoef(y_train, y_train_pred)[0, 1]
    test_correlation = np.corrcoef(y_test, y_test_pred)[0, 1]
    train_mse = mean_squared_error(y_train, y_train_pred)
    test_mse = mean_squared_error(y_test, y_test_pred)
    
    # Calculate AUC for both sets using upper 1/3 vs lower 2/3
    train_threshold = np.percentile(y_train, 66.67)
    test_threshold = np.percentile(y_test, 66.67)
    y_train_bin = (y_train > train_threshold).astype(int)
    y_test_bin = (y_test > test_threshold).astype(int)
    train_auc = roc_auc_score(y_train_bin, y_train_pred)
    test_auc = roc_auc_score(y_test_bin, y_test_pred)
    
    # Save metrics
    metrics = {
        'signature': target_signature,
        'train_correlation': train_correlation,
        'test_correlation': test_correlation,
        'train_mse': train_mse,
        'test_mse': test_mse,
        'train_auc': train_auc,
        'test_auc': test_auc
    }
    
    # Create plots directory if it doesn't exist
    os.makedirs('results/plots', exist_ok=True)
    
    # Save metrics
    metrics_file = 'results/plots/model_metrics.csv'
    metrics_df = pd.DataFrame([metrics])
    
    # Append to existing results file
    if os.path.exists(metrics_file):
        existing_metrics = pd.read_csv(metrics_file)
        metrics_df = pd.concat([existing_metrics, metrics_df], ignore_index=True)
    
    metrics_df.to_csv(metrics_file, index=False)
    print(f"Model metrics saved to {metrics_file}")
    
    # Plot performance
    plot_model_performance(y_test, y_test_pred, target_signature, test_correlation, test_mse)
    
    print(f"Train correlation: {train_correlation:.3f}")
    print(f"Test correlation: {test_correlation:.3f}")
    print(f"Train MSE: {train_mse:.3f}")
    print(f"Test MSE: {test_mse:.3f}")
    print(f"Train AUC: {train_auc:.3f}")
    print(f"Test AUC: {test_auc:.3f}")
    
    return model, test_correlation, test_mse, test_auc

def main():
    # Load data
    signature_score, segment_score = load_data()
    
    # Process each signature
    print("\nBuilding prediction models for all signatures...")
    all_results = []
    total_signatures = len(signature_score.index)
    
    for idx, signature in enumerate(signature_score.index, 1):
        print(f"\nProcessing signature {idx}/{total_signatures}: {signature}")
        try:
            model, correlation, mse, roc_auc = build_prediction_model(signature_score, segment_score, signature)
            
            if model is not None:
                # Store results
                all_results.append({
                    'signature': signature,
                    'correlation': correlation,
                    'mse': mse,
                    'roc_auc': roc_auc
                })
        except Exception as e:
            print(f"Error processing signature {signature}: {str(e)}")
            print("Continuing with next signature...")
            continue
    
    if not all_results:
        print("\nNo successful models were built!")
        return
        
    # Save all results
    results_summary = pd.DataFrame(all_results)
    results_summary.to_csv('results/modeling_results.csv', index=False)
    print(f"\nSaved modeling results to results/modeling_results.csv")
    print(f"Successfully built {len(all_results)} out of {total_signatures} models")
    
    # Print summary of results
    print("\nModeling Results Summary:")
    print(results_summary.sort_values('correlation', ascending=False).head(10))

if __name__ == "__main__":
    main() 


    #def perform_association_test(signature_score, segment_score):
#    """Perform association tests between signatures and segments"""
#    print("\nPerforming association tests...")
#    results = []
#    
#    # Convert to numpy arrays for faster computation
#    sig_array = signature_score.values
#    seg_array = segment_score.values
#    
#    print(f"Signature array shape: {sig_array.shape}")
#    print(f"Segment array shape: {seg_array.shape}")
#    
#    for i, signature in enumerate(signature_score.index):
#        if i % 10 == 0:  # Progress update
#            print(f"Processing signature {i+1}/{len(signature_score.index)}")
#        
#        score = sig_array[i].astype(float)
#        
#        for j, segment in enumerate(segment_score.index):
#            CN = seg_array[j].astype(float)
#            
#            # Remove missing values
#            mask = (~np.isnan(score)) & (~np.isnan(CN))
#            if np.sum(mask) < 3:
#                continue
#            score_valid = score[mask]
#            CN_valid = CN[mask]
#            
#            try:
#                # Spearman correlation
#                spearman_cor, spearman_p = stats.spearmanr(score_valid, CN_valid)
#                
#                # Linear model
#                X = sm.add_constant(CN_valid)
#                y = score_valid
#                model = sm.OLS(y, X).fit()
#                beta = model.params[1]  # index 1 for CN
#                p_value = model.pvalues[1]
#                
#                results.append({
#                    'signature': signature,
#                    'segment': segment,
#                    'spearman_cor': spearman_cor,
#                    'spearman_p': spearman_p,
#                    'beta': beta,
#                    'p_value': p_value,
#                    'n_samples': len(score_valid)
#                })
#            except Exception as e:
#                print(f"Error processing signature {signature} and segment {segment}: {str(e)}")
#                continue
#    
#    if not results:
#        print("No valid associations found!")
#        return pd.DataFrame()
#    
#    # Convert to DataFrame and save
#    results_df = pd.DataFrame(results)
#    print(f"\nFound {len(results_df)} valid associations")
#    
#    # Check if results contain required columns
#    required_columns = ['signature', 'segment', 'spearman_cor', 'spearman_p', 'beta', 'p_value', 'n_samples']
#    missing_columns = [col for col in required_columns if col not in results_df.columns]
#    if missing_columns:
#        print(f"Warning: Missing columns in results: {missing_columns}")
#        return results_df
#    
#    results_df = results_df.sort_values('spearman_p')
#    
#    # Save complete results
#    results_df.to_csv('results/association_results.csv', index=False)
#    print(f"\nSaved {len(results_df)} association results")
#    
#    # Print most significant correlations
#    print("\nTop 10 most significant correlations:")
#    print(results_df.head(10))
#    
#    return results_df