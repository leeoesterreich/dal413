import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, mean_squared_error
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed, dump, load
from sklearn.impute import SimpleImputer
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """Load the signature scores and segment scores"""
    print("Loading data...")
    signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/signature_score.pkl")
    segment_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/segment_score.pkl")
    
    # Find common samples
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

def get_lambda_range(X, y, alpha):
    """Determine min and max lambda values for a given alpha"""
    model = ElasticNet(alpha=alpha, l1_ratio=0.5, max_iter=10000)
    model.fit(X, y)
    lambda_max = np.max(np.abs(model.coef_))
    lambda_min = lambda_max * 0.01
    return np.logspace(np.log10(lambda_min), np.log10(lambda_max), 100)  # 100 lambda values

def monte_carlo_cv(X, y, alpha, lambda_val, n_splits=200, train_size=0.7):
    """Perform Monte Carlo cross validation with 200 splits"""
    accuracies = []
    for _ in range(n_splits):
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, train_size=train_size, random_state=42
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
    os.makedirs('results/plots', exist_ok=True)
    
    # Calculate ROC curve using upper 1/3 vs lower 2/3
    threshold = np.percentile(y_true, 66.67)
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

def build_prediction_model(signature_score, segment_score, target_signature):
    """Build prediction model for a signature using Elastic Net"""
    print(f"\nBuilding prediction model for {target_signature}")
    
    # Prepare data
    X = segment_score.values.T
    y = signature_score.loc[target_signature].values
    
    # Convert to numeric types
    X = X.astype(np.float64)
    y = y.astype(np.float64)
    
    print(f"Initial sample count: {len(X)}")
    
    # Handle missing values by imputation
    X = impute_preserve_all_columns(X)
    y = np.nan_to_num(y, nan=np.nanmean(y))  # Simple imputation for y
    
    print(f"Sample count after imputation: {len(X)}")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    # Parameter tuning
    print("Performing parameter tuning...")
    alphas = np.arange(0.1, 1.1, 0.1)  # Range from 0.1 to 1.0 by 0.1
    
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
    
    # Save the model and preprocessing information
    model_info = {
        'model': model,
        'feature_names': segment_score.index.tolist(),  # Save all feature names
        'best_params': best_params,
        'threshold_percentile': 66.67  # Upper 1/3 threshold
    }
    
    # Create directory if it doesn't exist
    os.makedirs('results/models', exist_ok=True)
    
    # Save the model and preprocessing information
    safe_name = target_signature.replace('/', '_').replace('\\', '_')
    dump(model_info, f'results/models/{safe_name}_model.joblib')
    print(f"Saved model and preprocessing information to results/models/{safe_name}_model.joblib")
    
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
        'test_auc': test_auc,
        'best_alpha': best_params[0],
        'best_lambda': best_params[1]
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
    
    # Define the target signature
    target_signature = "GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN"
    print(f"\nBuilding model for {target_signature}")
    
    try:
        # Build prediction model for the target signature
        model, correlation, mse, auc_score = build_prediction_model(signature_score, segment_score, target_signature)
        
        if auc_score is not None:
            print(f"\nModel built successfully for {target_signature}")
            print(f"Test Correlation: {correlation:.3f}")
            print(f"Test MSE: {mse:.3f}")
            print(f"Test AUC: {auc_score:.3f}")
            
            # Save model results
            model_results = [{
                'signature': target_signature,
                'correlation': correlation,
                'mse': mse,
                'auc': auc_score
            }]
            
            # Save results
            results_df = pd.DataFrame(model_results)
            results_file = f'results/{target_signature.replace("/", "_")}_model_results.csv'
            results_df.to_csv(results_file, index=False)
            print(f"\nModel results saved to {results_file}")
            
        else:
            print(f"\nFailed to build model for {target_signature}")
            
    except Exception as e:
        print(f"Error building model for {target_signature}: {str(e)}")

if __name__ == "__main__":
    main() 