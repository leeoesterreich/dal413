# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet, ElasticNetCV, MultiTaskElasticNetCV
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed, dump, load
from sklearn.impute import SimpleImputer
import warnings
import pickle
from sklearn.linear_model import LogisticRegressionCV
from scipy.stats import pearsonr

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("training_data/rna_signature_score_median_no_norm.pkl")
    segment_score = pd.read_pickle("training_data/cna_segment_score_mean_no_norm.pkl")
    
    print("\nSample name formats:")
    print("RNA signature score sample names (first 5):", list(signature_score.columns[:5]))
    print("CNA segment score sample names (first 5):", list(segment_score.columns[:5]))
    
    # Convert CNA sample names to match RNA format (replace . with -)
    segment_score.columns = segment_score.columns.str.replace('.', '-')
    
    # Find common samples using exact names
    common_samples = sorted(set(signature_score.columns) & set(segment_score.columns))
    print("\nFound {} common samples".format(len(common_samples)))
    if len(common_samples) > 0:
        print("First 5 common samples:", common_samples[:5])
    
    # Align by common samples
    signature_score = signature_score[common_samples]
    segment_score = segment_score[common_samples]
    
    return signature_score, segment_score

def get_lambda_range(X, y, alpha):
    """Determine min and max lambda values for a given alpha using glmnet-like (ElasticNet) approach"""
    # Fit ElasticNetCV to get lambda range
    model = ElasticNetCV(alphas=[alpha], l1_ratio=0.5, cv=5, max_iter=10000)
    model.fit(X, y)
    lambda_max = np.max(np.abs(model.coef_))
    lambda_min = lambda_max * 0.01
    return np.logspace(np.log10(lambda_min), np.log10(lambda_max), 100)  # 100 lambda values

def normalize_data(X_train, X_test=None):
    """Normalize data to match R's glmnet standardize=TRUE (ddof=1)"""
    X_centered = X_train - X_train.mean(axis=0)
    std = X_train.std(axis=0, ddof=1)
    # Handle zero standard deviations
    std[std == 0] = 1.0
    X_scaled = X_centered / std
    
    if X_test is not None:
        X_test_centered = X_test - X_train.mean(axis=0)
        X_test_scaled = X_test_centered / std
        return X_scaled, X_test_scaled
    return X_scaled

def clean_data(X):
    """Clean data by handling NaN and infinite values"""
    # Convert to numpy array if not already and ensure float type
    X = np.array(X, dtype=float)
    
    # Replace inf with NaN
    X = np.where(np.isinf(X), np.nan, X)
    
    # Check for all-NaN columns
    nan_cols = np.all(np.isnan(X), axis=0)
    if np.any(nan_cols):
        print("Warning: Found {} columns with all NaN values".format(np.sum(nan_cols)))
        # Replace all-NaN columns with zeros
        X[:, nan_cols] = 0
    
    # Replace remaining NaN with column means
    imputer = SimpleImputer(strategy='mean')
    X = imputer.fit_transform(X)
    
    # Final check for any remaining NaN or inf
    if np.any(np.isnan(X)) or np.any(np.isinf(X)):
        print("Warning: NaN or inf values remain after cleaning")
        # Replace any remaining NaN or inf with 0
        X = np.where(np.isnan(X) | np.isinf(X), 0, X)
    
    return X

def caret_wrap(trainX, trainY, testX, testY, alphas=np.logspace(-3, 1, 100), l1_ratio=0.5, cv=5):
    """Train an ElasticNet model with cross-validation"""
    print("  Cleaning training data...")
    trainX = clean_data(trainX)
    testX = clean_data(testX)
    
    print("  Normalizing data...")
    scaler = StandardScaler()
    trainX_scaled = scaler.fit_transform(trainX)
    testX_scaled = scaler.transform(testX)
    
    print("  Training ElasticNetCV model...")
    model = ElasticNetCV(
        l1_ratio=l1_ratio,
        alphas=alphas,
        cv=cv,
        random_state=42,
        max_iter=10000,  # Increased from 1000
        tol=1e-7,  # More stringent tolerance
        selection='random'  # Can help with convergence
    )
    model.fit(trainX_scaled, trainY)
    print("  Model training completed. Best alpha:", model.alpha_)
    
    return model, trainX_scaled, testX_scaled

def monte_carlo_cv(X, y, signature_name, n_splits=200, train_size=0.7, alphas=None):
    """Perform Monte Carlo cross validation with 200 splits for a single signature"""
    if alphas is None:
        alphas = np.logspace(-1, 0, 10)  # 0.1 to 1.0
    
    print("\nProcessing signature: {}".format(signature_name))
    print("Starting Monte Carlo cross-validation with {} splits".format(n_splits))
    
    # Ensure y is in the correct shape (n_samples,)
    y = y.ravel()
    
    cv_scores = []
    
    for i in range(n_splits):
        # Split data with consistent random state
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, train_size=train_size, random_state=i
        )
        
        # Clean and normalize data
        X_train = clean_data(X_train)
        X_test = clean_data(X_test)
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train model
        model = ElasticNetCV(
            alphas=alphas,
            l1_ratio=0.5,
            cv=5,
            max_iter=10000,
            tol=1e-4,
            selection='random',
            random_state=i
        )
        model.fit(X_train_scaled, y_train)
        
        # Calculate score
        score = model.score(X_test_scaled, y_test)
        cv_scores.append({
            'split': i,
            'alpha': model.alpha_,
            'score': score
        })
        
        if (i + 1) % 20 == 0:
            print("Completed {}/{} splits".format(i + 1, n_splits))
    
    # Calculate final statistics
    mean_score = np.mean([s['score'] for s in cv_scores])
    best_alpha = cv_scores[np.argmax([s['score'] for s in cv_scores])]['alpha']
    
    print("Monte Carlo cross-validation completed for signature: {}".format(signature_name))
    print("Mean score: {:.3f}, Best alpha: {}".format(mean_score, best_alpha))
    
    return pd.DataFrame(cv_scores)

def evaluate_params(args):
    """Helper function for parallel parameter evaluation"""
    X, y, alpha, lambda_val = args
    return (alpha, lambda_val, monte_carlo_cv(X, y, alpha, lambda_val))

def plot_model_performance(y_true, y_pred, signature_name, correlation, mse):
    """Plot model performance metrics"""
    # Calculate ROC curve using upper 1/3 and lower 2/3 thresholds
    thresholds = np.percentile(y_true, [33.33, 66.67])
    y_true_binary = np.where(y_true >= thresholds[1], 1, 0)  # Upper 1/3 is positive class
    
    fpr, tpr, _ = roc_curve(y_true_binary, y_pred)
    roc_auc = auc(fpr, tpr)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot actual vs predicted
    ax1.scatter(y_true, y_pred, alpha=0.5)
    ax1.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', lw=2)
    ax1.set_xlabel('Actual')
    ax1.set_ylabel('Predicted')
    ax1.set_title('Actual vs Predicted\nCorrelation: {:.3f}, MSE: {:.3f}'.format(correlation, mse))
    
    # Plot ROC curve
    ax2.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (AUC = {:.3f})'.format(roc_auc))
    ax2.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('False Positive Rate')
    ax2.set_ylabel('True Positive Rate')
    ax2.set_title('Receiver Operating Characteristic')
    ax2.legend(loc="lower right")
    
    plt.tight_layout()
    plt.savefig('results/plots/{}_performance.png'.format(signature_name.replace('/', '_').replace('\\', '_')))
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

def build_prediction_model(X, y, signature_name, best_alpha, feature_names):
    """Build and evaluate the final model"""
    print("\nBuilding final model for signature: {}".format(signature_name))
    
    # Ensure y is in the correct shape (n_samples, 1)
    y = y.reshape(-1, 1)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Clean and normalize data
    X_train = clean_data(X_train)
    X_test = clean_data(X_test)
    y_train = clean_data(y_train)
    y_test = clean_data(y_test)
    
    # Train model
    model = ElasticNetCV(alphas=[best_alpha], cv=5, random_state=42)
    model.fit(X_train, y_train.ravel())  # Use ravel() when fitting
    
    # Make predictions
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)
    
    # Calculate metrics
    train_correlation, _ = pearsonr(y_train.ravel(), y_train_pred.ravel())
    test_correlation, _ = pearsonr(y_test.ravel(), y_test_pred.ravel())
    train_mse = mean_squared_error(y_train.ravel(), y_train_pred.ravel())
    test_mse = mean_squared_error(y_test.ravel(), y_test_pred.ravel())
    
    # Plot performance and get ROC AUC
    train_roc_auc = plot_model_performance(y_train.ravel(), y_train_pred, signature_name + '_train', train_correlation, train_mse)
    test_roc_auc = plot_model_performance(y_test.ravel(), y_test_pred, signature_name + '_test', test_correlation, test_mse)
    
    # Create directories if they don't exist
    os.makedirs('results/models', exist_ok=True)
    os.makedirs('results/plots', exist_ok=True)
    
    # Save model
    model_info = {
        'model': model,
        'feature_names': feature_names,
        'best_alpha': best_alpha,
        'train_correlation': train_correlation,
        'test_correlation': test_correlation,
        'train_mse': train_mse,
        'test_mse': test_mse,
        'train_roc_auc': train_roc_auc,
        'test_roc_auc': test_roc_auc
    }
    
    # Save model
    safe_name = signature_name.replace('/', '_').replace('\\', '_')
    dump(model_info, 'results/models/{}_model.pkl'.format(safe_name))
    
    # Create metrics dictionary
    metrics = {
        'signature': signature_name,
        'train_correlation': train_correlation,
        'test_correlation': test_correlation,
        'train_mse': train_mse,
        'test_mse': test_mse,
        'train_roc_auc': train_roc_auc,
        'test_roc_auc': test_roc_auc,
        'best_alpha': best_alpha
    }
    
    # Save metrics to CSV
    metrics_df = pd.DataFrame([metrics])
    metrics_file = 'results/all_model_metrics.csv'
    
    # Append to existing file or create new one
    if os.path.exists(metrics_file):
        metrics_df.to_csv(metrics_file, mode='a', header=False, index=False)
    else:
        metrics_df.to_csv(metrics_file, index=False)
    
    print("Saved metrics for signature: {}".format(signature_name))
    
    return metrics

def main():
    """Main function to run the model training pipeline"""
    # Create results directory if it doesn't exist
    results_dir = "results_v2_original_data"
    os.makedirs(results_dir, exist_ok=True)
    
    # Load and align data
    signature_score, segment_score = load_data()
    print("\nData loaded and aligned:")
    print("RNA signature scores shape:", signature_score.shape)
    print("CNA segment scores shape:", segment_score.shape)
    
    # Save aligned data
    signature_score.to_pickle(os.path.join(results_dir, "aligned_rna_signatures.pkl"))
    segment_score.to_pickle(os.path.join(results_dir, "aligned_cna_segments.pkl"))
    
    # Process each RNA signature
    for signature_name in signature_score.index:
        print("\nProcessing signature:", signature_name)
        
        # Get target values (RNA signature scores)
        y = signature_score.loc[signature_name].values
        
        # Get features (CNA segment scores)
        X = segment_score.values.T  # Transpose to get samples x features
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
        
        # Train model
        print("Training model...")
        model, X_train_scaled, X_test_scaled = caret_wrap(
            X_train, y_train, X_test, y_test,
            alphas=np.logspace(-3, 1, 100),  # More granular alpha range
            l1_ratio=0.5,
            cv=5
        )
        
        # Make predictions
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        
        # Calculate performance metrics
        train_corr, _ = pearsonr(y_train, y_pred_train)
        test_corr, _ = pearsonr(y_test, y_pred_test)
        train_mse = mean_squared_error(y_train, y_pred_train)
        test_mse = mean_squared_error(y_test, y_pred_test)
        
        print("Train correlation: {:.3f}, MSE: {:.3f}".format(train_corr, train_mse))
        print("Test correlation: {:.3f}, MSE: {:.3f}".format(test_corr, test_mse))
        
        # Save model and results
        model_file = os.path.join(results_dir, "{}_model.joblib".format(signature_name))
        results_file = os.path.join(results_dir, "{}_results.pkl".format(signature_name))
        
        dump(model, model_file)
        
        # Calculate ROC curves
        # For training set
        train_thresh = np.percentile(y_train, [33.33, 66.67])
        y_train_binary = np.where(y_train >= train_thresh[1], 1, 0)
        fpr_train, tpr_train, _ = roc_curve(y_train_binary, y_pred_train)
        roc_auc_train = auc(fpr_train, tpr_train)
        
        # For test set
        test_thresh = np.percentile(y_test, [33.33, 66.67])
        y_test_binary = np.where(y_test >= test_thresh[1], 1, 0)
        fpr_test, tpr_test, _ = roc_curve(y_test_binary, y_pred_test)
        roc_auc_test = auc(fpr_test, tpr_test)
        
        results = {
            'train_corr': train_corr,
            'test_corr': test_corr,
            'train_mse': train_mse,
            'test_mse': test_mse,
            'train_roc_auc': roc_auc_train,
            'test_roc_auc': roc_auc_test,
            'feature_importance': pd.Series(model.coef_, index=segment_score.index),
            'y_train': y_train,
            'y_pred_train': y_pred_train,
            'y_test': y_test,
            'y_pred_test': y_pred_test
        }
        with open(results_file, 'wb') as f:
            pickle.dump(results, f)
        
        # Plot results with ROC curves
        plt.figure(figsize=(15, 5))
        
        # Training set scatter plot
        plt.subplot(1, 3, 1)
        plt.scatter(y_train, y_pred_train, alpha=0.5)
        plt.plot([min(y_train), max(y_train)], [min(y_train), max(y_train)], 'r--')
        plt.xlabel('Actual')
        plt.ylabel('Predicted')
        plt.title('Training Set\nCorr: {:.3f}, MSE: {:.3f}'.format(train_corr, train_mse))
        
        # Test set scatter plot
        plt.subplot(1, 3, 2)
        plt.scatter(y_test, y_pred_test, alpha=0.5)
        plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--')
        plt.xlabel('Actual')
        plt.ylabel('Predicted')
        plt.title('Test Set\nCorr: {:.3f}, MSE: {:.3f}'.format(test_corr, test_mse))
        
        # ROC curves
        plt.subplot(1, 3, 3)
        plt.plot(fpr_train, tpr_train, 'b-', label='Train (AUC = {:.3f})'.format(roc_auc_train))
        plt.plot(fpr_test, tpr_test, 'r-', label='Test (AUC = {:.3f})'.format(roc_auc_test))
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curves')
        plt.legend(loc="lower right")
        
        plt.tight_layout()
        plt.savefig(os.path.join(results_dir, "{}_performance.png".format(signature_name)))
        plt.close()
        
    print("\nModel training completed. Results saved in:", results_dir)

if __name__ == '__main__':
    main() 