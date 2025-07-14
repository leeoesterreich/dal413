# -*- coding: utf-8 -*-

"""
Best-performing model settings (2024-03-xx):
- No scaling for both X and Y data
- L1 ratio = 0.02
- Alpha = 0.154413 (found via CV)
- These settings produced feature landscapes most similar to the original paper
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split, StratifiedKFold, ShuffleSplit
from sklearn.linear_model import ElasticNet, ElasticNetCV, MultiTaskElasticNetCV, LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed, dump, load
from sklearn.impute import SimpleImputer
import warnings
import pickle
from scipy.stats import pearsonr
from sklearn.model_selection import GridSearchCV
import time

warnings.filterwarnings('ignore')

# Configuration
USE_TRIAL_SIGNATURES = True
TRIAL_SIGNATURES = [
    "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP",
    "GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN",
    "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450"
] # User-provided exact names
CV_FOLDS_OUTER = 1 # For train/test split, not Monte Carlo for this script
CV_FOLDS_INNER_ELASTICNET = 4 # 4-fold CV for ElasticNetCV hyperparameter tuning
L1_RATIOS_TO_TUNE = np.arange(0.1, 1.01, 0.1).tolist() # L1 ratios from 0.1 to 1.0 by 0.1
N_ALPHAS_CV = 100 # Ensure 100 alphas (lambdas) are used in CV for each L1 ratio
RANDOM_STATE = 42
TEST_SIZE_SPLIT = 0.25
MANUALLY_SET_ALPHA_FOR_FINAL_MODEL = False # Let CV find best alpha for each L1 ratio
MANUAL_ALPHA_VALUE = 0.4 # This will not be used as MANUALLY_SET_ALPHA_FOR_FINAL_MODEL is False

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    # Ensure these paths are correct for your environment
    try:
        signature_score = pd.read_pickle("scaling/training_data/rna_signature_score_median_no_norm.pkl")
        segment_score = pd.read_pickle("scaling/training_data/cna_segment_score_mean_no_norm.pkl")
    except FileNotFoundError as e:
        print(f"Error loading data: {e}")
    
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

def binarize_score(score_array, percentile_threshold=66.67):
    """Binarizes a score array based on a percentile threshold."""
    if score_array.ndim > 1 and score_array.shape[1] == 1:
        score_array = score_array.ravel()
    threshold_val = np.percentile(score_array[~np.isnan(score_array)], percentile_threshold)
    return (score_array >= threshold_val).astype(int)

def monte_carlo_cv_classification(X_input, y_input_continuous, signature_name, 
                                  n_splits=200, train_proportion=0.7, 
                                  l1_ratios=[0.5], C_values=np.logspace(-2, 2, 10)):
    """
    Perform Monte Carlo cross-validation for Logistic Regression with ElasticNet penalty.
    Finds the best C (inverse of regularization strength) for a given l1_ratio.
    X_input: features (samples x features)
    y_input_continuous: continuous target variable (samples,)
    """
    print(f"\nStarting Monte Carlo CV for {signature_name} (Classification)")
    
    y_input_binary = binarize_score(y_input_continuous) # Binarize y for classification

    # Stratified ShuffleSplit is good for maintaining class proportions
    # The R script used balancedstratification, StratifiedShuffleSplit is a good sklearn equivalent
    # For direct n_splits without grouping, ShuffleSplit might be closer if R's LGOCV didn't use groups from balancing_variables
    # If R's balancedstratification was complex, an exact match might be hard.
    # Using StratifiedKFold if we want distinct folds for each of the n_splits.
    # However, R's LGOCV with number=200 means 200 random 70/30 splits. ShuffleSplit is best for this.
    
    cv_splitter = ShuffleSplit(n_splits=n_splits, train_size=train_proportion, random_state=42)
    
    best_avg_auc = -1
    best_C_overall = None
    
    # Loop over l1_ratios if you want to tune it as well, current R script seems to fix it at 0.5
    # For simplicity, assuming fixed l1_ratio as per R script's glmnet call (often default 0.5 for caret if not specified broadly)
    # The R script uses alpha = seq(0.1,0.9,by=0.1) in caret_wrap's expand.grid for glmnet,
    # which is the l1_ratio in scikit-learn.
    # So we should iterate over l1_ratios as well.
    
    # Let's refine this: the R `Elastic_Net_modeling.R` calls `caret_wrap` which has a loop for alpha (l1_ratio)
    # and lambda (C). `glmnet_obj$bestTune` would have the best alpha and lambda.
    # For Python, LogisticRegressionCV can tune C for fixed l1_ratios, or GridSearchCV can tune both.
    
    # To simplify and align with finding one "best_alpha" (lambda in glmnet context) as in previous Python structure:
    # We find the best C for each l1_ratio, then pick the l1_ratio/C pair that gives best mean AUC.

    # Simplified: Tune C for a fixed l1_ratio (e.g., 0.5) or iterate l1_ratios too.
    # The R script Elastic_Net_modeling.R seems to let caret find best alpha (l1_ratio) and lambda.
    # So, we should use GridSearchCV.

    param_grid = {'C': C_values, 'l1_ratio': l1_ratios}
    
    # Using StratifiedShuffleSplit for the inner CV of GridSearchCV
    inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42) 
    
    model_for_grid = LogisticRegression(penalty='elasticnet', solver='saga', max_iter=5000, tol=1e-4, random_state=42) # Increased max_iter

    # Accumulate scores from n_splits outer loops
    # The R script does 200 rounds of LGOCV for parameter tuning.
    # This means it gets 200 estimates of performance for each param combo.
    # GridSearchCV does one set of CV folds.
    # To truly match R's 200 rounds for tuning, we'd need to wrap GridSearchCV in a loop
    # or use a more complex CV setup.
    
    # Simpler approach: GridSearchCV finds best params based on its inner CV.
    # If R's caret_wrap uses number=200 for LGOCV *within* its tuning for each (alpha,lambda) pair,
    # that's very computationally intensive.
    # More likely, caret::train uses number=200 to mean 200 resamples to evaluate each hyperparameter point in the grid.

    # Let's assume one robust GridSearchCV is sufficient to find good C and l1_ratio
    
    X_cleaned = clean_data(X_input) # Clean once
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_cleaned)

    grid_search = GridSearchCV(model_for_grid, param_grid, cv=inner_cv, scoring='roc_auc', n_jobs=-1)
    grid_search.fit(X_scaled, y_input_binary)
    
    best_params_from_grid = grid_search.best_params_
    best_avg_auc_from_grid = grid_search.best_score_

    print(f"GridSearchCV for {signature_name} completed. Best params: {best_params_from_grid}, Best CV AUC: {best_avg_auc_from_grid:.4f}")
    
    return best_params_from_grid['C'], best_params_from_grid['l1_ratio']

def build_final_model_and_evaluate_classification(X_train_scaled, X_test_scaled, 
                                                  y_train_binary, y_test_binary, # These are binary
                                                  signature_name, best_C, best_l1_ratio, 
                                                  feature_names):
    """
    Builds and evaluates the final Logistic Regression with ElasticNet penalty model.
    """
    print(f"\nBuilding final classification model for {signature_name} with C={best_C}, l1_ratio={best_l1_ratio}")

    final_model = LogisticRegression(C=best_C, l1_ratio=best_l1_ratio, penalty='elasticnet', 
                                     solver='saga', max_iter=10000, tol=1e-7, random_state=42) # Increased max_iter
    final_model.fit(X_train_scaled, y_train_binary)

    # Predictions (probabilities for class 1)
    y_train_pred_proba = final_model.predict_proba(X_train_scaled)[:, 1]
    y_test_pred_proba = final_model.predict_proba(X_test_scaled)[:, 1]

    # Metrics
    train_auc = roc_auc_score(y_train_binary, y_train_pred_proba)
    test_auc = roc_auc_score(y_test_binary, y_test_pred_proba)
    
    # For "correlation" and "MSE", it's less direct for classification.
    # R might have reported pseudo-R2 or other metrics.
    # We can report accuracy or log-loss. For now, focusing on AUC and coeffs.
    train_accuracy = accuracy_score(y_train_binary, final_model.predict(X_train_scaled))
    test_accuracy = accuracy_score(y_test_binary, final_model.predict(X_test_scaled))

    print(f"  Train AUC: {train_auc:.3f}, Test AUC: {test_auc:.3f}")
    print(f"  Train Accuracy: {train_accuracy:.3f}, Test Accuracy: {test_accuracy:.3f}")

    coefficients = final_model.coef_[0] # Coeffs for class 1

    metrics = {
        'signature': signature_name,
        'train_auc': train_auc,
        'test_auc': test_auc,
        'train_accuracy': train_accuracy,
        'test_accuracy': test_accuracy,
        # Reporting dummy values for correlation and MSE for now, or remove them
        'train_correlation': np.nan, 
        'test_correlation': np.nan,
        'train_mse': np.nan,
        'test_mse': np.nan,
        'coefficients': coefficients,
        'intercept': final_model.intercept_[0]
    }
    
    # Your plot_model_performance might need adjustment for classification
    # It was plotting y_true vs y_pred (continuous) and ROC.
    # For classification, scatter plot of y_true vs y_pred_proba might not be as intuitive as ROC.
    # plot_model_performance_classification(y_test_binary, y_test_pred_proba, signature_name, test_auc)
    
    return metrics

def clean_data_array(X_arr):
    """Clean numpy array data by handling NaN and infinite values for features (X)."""
    X_arr = np.array(X_arr, dtype=float)
    X_arr[np.isinf(X_arr)] = np.nan # Replace infinities with NaN
    
    # Impute NaNs column by column using column mean
    # If a whole column is NaN, it will be filled with 0 after this loop fails to find a mean.
    imputer = SimpleImputer(strategy='mean', missing_values=np.nan)
    X_arr_imputed = imputer.fit_transform(X_arr)
    
    # If any NaNs remain (e.g., a column was all NaN, so mean was NaN, imputer might leave it or error)
    # SimpleImputer should handle all-NaN columns by imputing with the mean of other columns if fit on full data,
    # or 0 if fit on a single all-NaN column. A final check is good.
    if np.any(np.isnan(X_arr_imputed)):
        # print("Warning: NaNs remain after SimpleImputer. Filling remaining NaNs with 0.")
        X_arr_imputed[np.isnan(X_arr_imputed)] = 0 # Fill any residual NaNs with 0
    
    return X_arr_imputed

def clean_target_array(y_arr):
    """Clean numpy array for target variable (y). Simple NaN to 0 for 1D array."""
    y_arr = np.array(y_arr, dtype=float)
    y_arr[np.isinf(y_arr) | np.isnan(y_arr)] = 0 # Replace inf or NaN with 0
    return y_arr

def train_evaluate_elasticnet_regression(X_train, y_train, X_test, y_test, 
                                         l1_ratios_to_tune, n_alphas_cv, cv_folds_inner, random_state, 
                                         rna_signature_name, scale_x_flag, scale_y_flag):
    print(f"    Preprocessing and Scaling (scale_x={scale_x_flag}, scale_y={scale_y_flag})...")
    # Clean data: drop columns with all NaNs or all zeros from training set, then apply to test
    na_cols_train = X_train.columns[X_train.isna().all()].tolist()
    zero_cols_train = X_train.columns[(X_train == 0).all()].tolist()
    cols_to_drop = list(set(na_cols_train + zero_cols_train))
    
    X_train_clean = X_train.drop(columns=cols_to_drop)
    X_test_clean = X_test.drop(columns=cols_to_drop)
    print(f"      Dropped {len(cols_to_drop)} columns with all NaNs/zeros based on training data.")
    print(f"      X_train_clean shape: {X_train_clean.shape}, X_test_clean shape: {X_test_clean.shape}")

    X_train_to_fit = X_train_clean.copy()
    X_test_to_transform = X_test_clean.copy()
    y_train_to_fit = y_train.copy()

    scaler_x_fitted = None
    if scale_x_flag:
        scaler_x_fitted = StandardScaler()
        X_train_to_fit = scaler_x_fitted.fit_transform(X_train_clean)
        X_test_to_transform = scaler_x_fitted.transform(X_test_clean)
        print("      X (CNA features) scaled.")
    else:
        X_train_to_fit = X_train_clean.values # ElasticNetCV expects numpy array
        X_test_to_transform = X_test_clean.values
        print("      X (CNA features) NOT scaled.")

    scaler_y_fitted = None
    if scale_y_flag:
        scaler_y_fitted = StandardScaler()
        y_train_to_fit = scaler_y_fitted.fit_transform(y_train.values.reshape(-1, 1)).ravel()
        print("      Y (RNA scores) scaled.")
    else:
        y_train_to_fit = y_train.values.ravel()
        print("      Y (RNA scores) NOT scaled.")

    print(f"    Tuning ElasticNetCV for {rna_signature_name}...")
    best_l1_ratio_cv = None
    best_alpha_cv = None
    best_mse_cv = float('inf')
    
    for i, l1_ratio in enumerate(l1_ratios_to_tune):
        print(f"      Tuning for L1 ratio: {l1_ratio:.2f} ({i+1}/{len(l1_ratios_to_tune)})...", flush=True)
        start_time_l1 = time.time()
        # Inner CV to find the best alpha for the current L1 ratio
        model_cv = ElasticNetCV(
            l1_ratio=l1_ratio,
            n_alphas=n_alphas_cv, # Number of alphas along the regularization path
            cv=cv_folds_inner,    # Use KFold object or integer for KFold
            random_state=random_state,
            max_iter=100000, # R glmnet default: 1e5
            tol=1e-7,      # R glmnet default: 1e-7
            eps=1e-4,    # R glmnet lambda.min.ratio default for nobs >= nvars (808 > 534)
            n_jobs=-1,     # Use all available CPUs
            selection='cyclic' # Default is 'cyclic'. 'random' can be faster for large datasets
        )
        model_cv.fit(X_train_to_fit, y_train_to_fit)
        current_best_mse_for_l1 = np.min(np.mean(model_cv.mse_path_, axis=1))
        current_best_alpha_for_l1 = model_cv.alpha_
        l1_time = (time.time() - start_time_l1)
        print(f"        L1 ratio {l1_ratio:.2f} completed in {l1_time:.2f}s. Best alpha for this L1: {current_best_alpha_for_l1:.6f}, MSE: {current_best_mse_for_l1:.6f}", flush=True)

        if current_best_mse_for_l1 < best_mse_cv:
            best_mse_cv = current_best_mse_for_l1
            best_l1_ratio_cv = l1_ratio
            best_alpha_cv = current_best_alpha_for_l1
    
    print(f"    CV Best parameters: l1_ratio = {best_l1_ratio_cv:.2f}, alpha = {best_alpha_cv:.6f}")

    # Determine alpha for final model
    alpha_for_final_model = MANUAL_ALPHA_VALUE if MANUALLY_SET_ALPHA_FOR_FINAL_MODEL else best_alpha_cv
    l1_ratio_for_final_model = best_l1_ratio_cv # Use CV best l1_ratio, or could also be made manual
    
    if MANUALLY_SET_ALPHA_FOR_FINAL_MODEL:
        print(f"    MANUALLY SETTING alpha for final model to: {alpha_for_final_model:.6f} (l1_ratio from CV: {l1_ratio_for_final_model:.2f})")
    else:
        print(f"    Using CV-selected alpha for final model: {alpha_for_final_model:.6f}")

    print(f"    Refitting final ElasticNet model with: l1_ratio={l1_ratio_for_final_model:.2f}, alpha={alpha_for_final_model:.6f}")
    final_model = ElasticNet(
        alpha=alpha_for_final_model, 
        l1_ratio=l1_ratio_for_final_model, 
        random_state=random_state,
        tol=1e-7, 
        max_iter=100000 
    )
    final_model.fit(X_train_to_fit, y_train_to_fit)
    print("      Final model refitted.")

    coeffs_from_model = final_model.coef_
    intercept_from_model = final_model.intercept_

    y_pred_test_on_model_scale = final_model.predict(X_test_to_transform)
    y_pred_train_on_model_scale = final_model.predict(X_train_to_fit)

    if scaler_y_fitted:
        y_pred_test_orig_scale = scaler_y_fitted.inverse_transform(y_pred_test_on_model_scale.reshape(-1,1)).ravel()
        y_pred_train_orig_scale = scaler_y_fitted.inverse_transform(y_pred_train_on_model_scale.reshape(-1,1)).ravel()
    else:
        y_pred_test_orig_scale = y_pred_test_on_model_scale
        y_pred_train_orig_scale = y_pred_train_on_model_scale

    coeffs_final_orig_xy = coeffs_from_model.copy()
    intercept_final_orig_xy = intercept_from_model 

    if scaler_x_fitted: 
        safe_scale_x = np.where(scaler_x_fitted.scale_ == 0, 1, scaler_x_fitted.scale_)
        coeffs_final_orig_xy = coeffs_final_orig_xy / safe_scale_x
        intercept_final_orig_xy = intercept_final_orig_xy - np.sum((coeffs_from_model * scaler_x_fitted.mean_) / safe_scale_x)

    if scaler_y_fitted: 
        coeffs_final_orig_xy = coeffs_final_orig_xy * scaler_y_fitted.scale_[0]
        intercept_final_orig_xy = intercept_final_orig_xy * scaler_y_fitted.scale_[0] + scaler_y_fitted.mean_[0]

    train_r2 = r2_score(y_train, y_pred_train_orig_scale)
    test_r2 = r2_score(y_test, y_pred_test_orig_scale)
    train_mse = mean_squared_error(y_train, y_pred_train_orig_scale)
    test_mse = mean_squared_error(y_test, y_pred_test_orig_scale)
    
    train_corr = np.corrcoef(y_train, y_pred_train_orig_scale)[0, 1]
    test_corr = np.corrcoef(y_test, y_pred_test_orig_scale)[0, 1]

    y_train_threshold = np.percentile(y_train, 100 * (2/3))
    y_test_threshold = np.percentile(y_test, 100 * (2/3))
    
    y_train_binary = (y_train >= y_train_threshold).astype(int)
    y_test_binary = (y_test >= y_test_threshold).astype(int)
    
    auc_train = roc_auc_score(y_train_binary, y_pred_train_orig_scale) if len(np.unique(y_train_binary)) > 1 else np.nan
    auc_test = roc_auc_score(y_test_binary, y_pred_test_orig_scale) if len(np.unique(y_test_binary)) > 1 else np.nan

    num_selected_features = np.sum(final_model.coef_ != 0)
    print(f"      Number of selected features by final model: {num_selected_features}")
    
    results = {
        'model': final_model,
        'coefficients_original_scale': coeffs_final_orig_xy, 
        'intercept_original_scale': intercept_final_orig_xy, 
        'scaler_x': scaler_x_fitted,
        'scaler_y': scaler_y_fitted,
        'feature_names': X_train_clean.columns.tolist(),
        'best_l1_ratio': l1_ratio_for_final_model, # Report the l1_ratio used for the final model
        'best_alpha': alpha_for_final_model,    # Report the alpha used for the final model
        'train_r2': train_r2, 'test_r2': test_r2,
        'train_mse': train_mse, 'test_mse': test_mse,
        'train_corr': train_corr, 'test_corr': test_corr,
        'auc_train': auc_train, 'auc_test': auc_test,
        'num_selected_features': num_selected_features,
        'y_pred_train': y_pred_train_orig_scale,
        'y_pred_test': y_pred_test_orig_scale,
        'y_train_true': y_train.values, 
        'y_test_true': y_test.values
    }
    return results

def plot_predictions_vs_actual(y_true, y_pred, title_prefix, output_dir, filename_suffix=""):
    plt.figure(figsize=(8, 8))
    plt.scatter(y_true, y_pred, alpha=0.5)
    plt.plot([min(y_true.min(),y_pred.min()), max(y_true.max(),y_pred.max())], [min(y_true.min(),y_pred.min()), max(y_true.max(),y_pred.max())], 'r--')
    plt.xlabel("Actual RNA Scores")
    plt.ylabel("Predicted RNA Scores")
    corr = np.corrcoef(y_true, y_pred)[0, 1]
    plt.title(f"{title_prefix} (Corr: {corr:.3f})")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f'{title_prefix}_{filename_suffix}_pred_vs_actual.png'))
    plt.close()

def hyperparameter_tuning_monte_carlo(X_input_cleaned_np, y_input_cleaned_np, signature_name,
                                      l1_ratio_grid, alpha_grid,
                                      n_monte_carlo_splits=200, train_size=0.75,
                                      random_state_base=42, results_dir="results/tuning"):
    """
    Performs hyperparameter tuning for ElasticNet using Monte Carlo cross-validation.
    THIS FUNCTION IS NO LONGER THE PRIMARY WORKFLOW IN main() but kept for potential future use.

    Args:
        X_input_cleaned_np (np.ndarray): Cleaned input features (samples, features).
        y_input_cleaned_np (np.ndarray): Cleaned input target variable (samples,).
        signature_name (str): Name of the signature being processed.
        l1_ratio_grid (list or np.ndarray): Grid of l1_ratios to test.
        alpha_grid (list or np.ndarray): Grid of alphas (penalty strength) to test.
        n_monte_carlo_splits (int): Number of Monte Carlo iterations.
        train_size (float): Proportion of the dataset to include in the train split.
        random_state_base (int): Base for random state generation for reproducibility.
        results_dir (str): Directory to save detailed tuning results per signature.

    Returns:
        tuple: (pd.DataFrame containing all tuning results, 
                dict containing the best parameters found)
    """
    results_data = []
    overall_best_auc_for_sig = -1.0
    overall_best_params_for_sig = {}

    print(f"\nStarting hyperparameter tuning for signature: {signature_name}")
    print(f"  L1 ratios: {l1_ratio_grid}")
    print(f"  Number of alphas: {len(alpha_grid)} (from ~{alpha_grid[0]:.2e} to ~{alpha_grid[-1]:.2e})")
    print(f"  Monte Carlo splits: {n_monte_carlo_splits}, Train size: {train_size}")

    total_combinations = len(l1_ratio_grid) * len(alpha_grid)
    current_combination_num = 0

    for l1_r in l1_ratio_grid:
        for alpha_v in alpha_grid:
            current_combination_num += 1
            
            mc_aucs_for_combo = []
            for i in range(n_monte_carlo_splits):
                current_random_state = random_state_base + i
                
                X_train_fold, X_test_fold, y_train_cont_fold, y_test_cont_fold = train_test_split(
                    X_input_cleaned_np, y_input_cleaned_np, train_size=train_size, random_state=current_random_state
                )

                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X_train_fold)
                X_test_scaled = scaler.transform(X_test_fold)
                
                try:
                    model = ElasticNet(alpha=alpha_v, l1_ratio=l1_r, 
                                       random_state=current_random_state,
                                       max_iter=10000, tol=1e-4)
                    model.fit(X_train_scaled, y_train_cont_fold)
                    y_pred_test_cont = model.predict(X_test_scaled)

                    if len(y_test_cont_fold) < 3 or len(np.unique(y_test_cont_fold)) < 2:
                        mc_aucs_for_combo.append(np.nan)
                        continue
                    
                    try:
                        threshold_low = np.percentile(y_test_cont_fold, 33.33)
                        threshold_high = np.percentile(y_test_cont_fold, 66.67)

                        if threshold_low == threshold_high:
                            median_val = np.median(y_test_cont_fold)
                            y_test_binary = np.where(y_test_cont_fold > median_val, 1, 0)
                        else:
                            y_test_binary = np.where(y_test_cont_fold >= threshold_high, 1, 0)

                        if len(np.unique(y_test_binary)) < 2:
                            mc_aucs_for_combo.append(np.nan)
                            continue
                        
                        auc_score = roc_auc_score(y_test_binary, y_pred_test_cont)
                        mc_aucs_for_combo.append(auc_score)
                    except ValueError: 
                        mc_aucs_for_combo.append(np.nan)
                        
                except Exception: 
                    mc_aucs_for_combo.append(np.nan)

            valid_aucs = [auc for auc in mc_aucs_for_combo if not np.isnan(auc)]
            avg_auc_for_combo = np.mean(valid_aucs) if valid_aucs else np.nan
            num_valid_splits = len(valid_aucs)

            if current_combination_num % 10 == 0 or current_combination_num == total_combinations:
                 print(f"    Processed combination {current_combination_num}/{total_combinations}: l1={l1_r:.2f}, alpha={alpha_v:.2e} -> Avg AUC: {avg_auc_for_combo if not np.isnan(avg_auc_for_combo) else 'NaN'} ({num_valid_splits} valid splits)")

            results_data.append({
                'l1_ratio': l1_r,
                'alpha': alpha_v,
                'average_auc': avg_auc_for_combo,
                'num_valid_mc_splits': num_valid_splits
            })

            if not np.isnan(avg_auc_for_combo) and avg_auc_for_combo > overall_best_auc_for_sig:
                overall_best_auc_for_sig = avg_auc_for_combo
                overall_best_params_for_sig = {'l1_ratio': l1_r, 'alpha': alpha_v, 'auc': avg_auc_for_combo, 'num_valid_splits': num_valid_splits}
    
    print() 
    results_df = pd.DataFrame(results_data)
    
    os.makedirs(results_dir, exist_ok=True)
    signature_file_name = signature_name.replace('/', '_').replace('\\\\', '_').replace(':', '_') # Corrected replacement for backslash
    detailed_results_path = os.path.join(results_dir, f"{signature_file_name}_tuning_details.csv")
    results_df.sort_values(by='average_auc', ascending=False).to_csv(detailed_results_path, index=False)
    print(f"  Detailed tuning results saved to: {detailed_results_path}")

    if overall_best_params_for_sig:
        print(f"  Best parameters for {signature_name}: L1={overall_best_params_for_sig['l1_ratio']:.2f}, Alpha={overall_best_params_for_sig['alpha']:.4e}, Avg AUC={overall_best_params_for_sig['auc']:.4f} ({overall_best_params_for_sig['num_valid_splits']} splits)")
    else:
        print(f"  No optimal parameters found for {signature_name} (all AUCs were NaN or no valid splits).")
        
    return results_df, overall_best_params_for_sig

def main():
    overall_start_time = time.time()
    
    # Define base output directory for new regression models
    base_output_dir = "scaling/adjusted_output_scripted" 
    os.makedirs(base_output_dir, exist_ok=True)
    
    # General results directories (can be consolidated or kept if other parts of script use them)
    os.makedirs("results/models", exist_ok=True) # Kept for now, though new main saves elsewhere
    os.makedirs("results/plots", exist_ok=True)  # Kept for now
    os.makedirs("results/metrics", exist_ok=True) # Kept for now

    signature_score, segment_score = load_data()
    X_source_df = segment_score.T  # CNA features: samples as rows, segments as columns

    if USE_TRIAL_SIGNATURES:
        # Filter signature_score to only include trial signatures if they exist in its index
        target_rna_signatures = [sig for sig in TRIAL_SIGNATURES if sig in signature_score.index]
        if not target_rna_signatures:
            print(f"ERROR: None of the TRIAL_SIGNATURES ({TRIAL_SIGNATURES}) found in the loaded signature_score data. Exiting.")
            print("First 20 available signatures in data:", signature_score.index.tolist()[:20])
            return
        elif len(target_rna_signatures) < len(TRIAL_SIGNATURES):
            print(f"WARNING: Not all TRIAL_SIGNATURES were found. Only processing: {target_rna_signatures}")
            found_set = set(target_rna_signatures)
            missing_set = [s for s in TRIAL_SIGNATURES if s not in found_set]
            print(f"Missing signatures: {missing_set}")
            print("First 20 available signatures in data for reference:", signature_score.index.tolist()[:20])
            
        print(f"\n--- Running REGRESSION MODEL TRAINING for TRIAL signatures: {target_rna_signatures} ---")
    else:
        target_rna_signatures = signature_score.index.tolist()
        print(f"\n--- Running REGRESSION MODEL TRAINING for ALL {len(target_rna_signatures)} signatures ---")

    # Scaling flags for the main regression task
    SCALE_X_FLAG = True # Keep scaling X enabled
    SCALE_Y_FLAG = True # Keep scaling Y enabled
    scale_info_str = f"{'scaleX' if SCALE_X_FLAG else 'noScaleX'}_{'scaleY' if SCALE_Y_FLAG else 'noScaleY'}"

    overall_metrics_list = []

    for rna_signature_name in target_rna_signatures:
        signature_process_start_time = time.time()
        print(f"\n===== Processing RNA Signature for Regression: {rna_signature_name} =====")

        if rna_signature_name not in signature_score.index:
            print(f"Signature {rna_signature_name} not found in signature_score.index. Skipping.")
            continue

        y_target_series = signature_score.loc[rna_signature_name]
        
        # Align X and Y by samples (index) before splitting
        X_aligned, y_aligned = X_source_df.align(y_target_series, join='inner', axis=0)

        if X_aligned.empty or y_aligned.empty:
            print(f"  No common samples found for signature {rna_signature_name} after alignment. Skipping.")
            continue
        
        print(f"  Data shapes for {rna_signature_name} after alignment: X_df={X_aligned.shape}, y_series=({len(y_aligned)},)")

        # Train/test split
        X_train, X_test, y_train, y_test = train_test_split(
            X_aligned, y_aligned, 
            test_size=TEST_SIZE_SPLIT, 
            random_state=RANDOM_STATE
        )
        print(f"  Train/Test split: X_train={X_train.shape}, X_test={X_test.shape}, y_train=({len(y_train)}), y_test=({len(y_test)})")

        # Train and evaluate the regression model
        regression_results = train_evaluate_elasticnet_regression(
            X_train, y_train, X_test, y_test,
            L1_RATIOS_TO_TUNE, # Global, user-specified
            N_ALPHAS_CV,
            CV_FOLDS_INNER_ELASTICNET,
            RANDOM_STATE,
            rna_signature_name,
            scale_x_flag=SCALE_X_FLAG,
            scale_y_flag=SCALE_Y_FLAG
        )

        test_corr = regression_results['test_corr']
        test_corr_str = f"{test_corr:.3f}".replace('.', 'p') # Format for directory name, e.g., 0p958

        # Sanitize rna_signature_name for use in file/directory paths
        safe_rna_name = rna_signature_name.replace('/', '_').replace('\\\\', '_').replace(':', '_').replace('=', '_eq_')
        
        # Construct specific output directory for this model run
        model_output_dirname = f"{safe_rna_name}.r={test_corr_str}_{scale_info_str}"
        model_output_path = os.path.join(base_output_dir, model_output_dirname)
        os.makedirs(model_output_path, exist_ok=True)
        print(f"  Outputs will be saved to: {model_output_path}")

        # Save model components
        model_components_to_save = {
            'model_object': regression_results['model'], # The trained sklearn.linear_model.ElasticNet object
            'coefficients_on_original_scale': regression_results['coefficients_original_scale'],
            'intercept_on_original_scale': regression_results['intercept_original_scale'],
            'selected_feature_names': regression_results['feature_names'],
            'l1_ratio_used': regression_results['best_l1_ratio'],
            'alpha_used': regression_results['best_alpha'],
            'scale_x_flag': SCALE_X_FLAG,
            'scale_y_flag': SCALE_Y_FLAG,
            'scaler_x_info': str(regression_results['scaler_x']) if regression_results['scaler_x'] else "None",
            'scaler_y_info': str(regression_results['scaler_y']) if regression_results['scaler_y'] else "None"
        }
        model_pkl_filename = f"{safe_rna_name}.r_{test_corr_str}_{scale_info_str}_elasticnet_model_components.pkl"
        dump(model_components_to_save, os.path.join(model_output_path, model_pkl_filename))
        print(f"    Saved model components to: {model_pkl_filename}")

        # Save features and coefficients to CSV
        features_df = pd.DataFrame({
            'feature': regression_results['feature_names'],
            'coefficient_original_scale': regression_results['coefficients_original_scale']
        })
        # Filter out features with zero coefficient if desired, or save all
        # features_df = features_df[features_df['coefficient_original_scale'] != 0] 
        features_csv_filename = f"{safe_rna_name}.r_{test_corr_str}_{scale_info_str}_elasticnet_model_features.csv"
        features_df.to_csv(os.path.join(model_output_path, features_csv_filename), index=False)
        print(f"    Saved features and coefficients to: {features_csv_filename}")
        
        # Plot predictions vs actual for TRAIN set
        plot_predictions_vs_actual(
            y_true=regression_results['y_train_true'], 
            y_pred=regression_results['y_pred_train'],
            title_prefix=f"{safe_rna_name}_Train",
            output_dir=model_output_path,
            filename_suffix=f"r={test_corr_str}_{scale_info_str}"
        )
        # Plot predictions vs actual for TEST set
        plot_predictions_vs_actual(
            y_true=regression_results['y_test_true'], 
            y_pred=regression_results['y_pred_test'],
            title_prefix=f"{safe_rna_name}_Test",
            output_dir=model_output_path,
            filename_suffix=f"r={test_corr_str}_{scale_info_str}"
        )
        print(f"    Saved prediction plots to directory.")

        # Store metrics for overall summary
        current_metrics = {
            'signature': rna_signature_name,
            'output_directory': model_output_path,
            'scale_x': SCALE_X_FLAG,
            'scale_y': SCALE_Y_FLAG,
            'l1_ratio_final': regression_results['best_l1_ratio'],
            'alpha_final': regression_results['best_alpha'],
            'train_corr': regression_results['train_corr'],
            'test_corr': regression_results['test_corr'],
            'train_r2': regression_results['train_r2'],
            'test_r2': regression_results['test_r2'],
            'train_mse': regression_results['train_mse'],
            'test_mse': regression_results['test_mse'],
            'auc_train': regression_results['auc_train'],
            'auc_test': regression_results['auc_test'],
            'num_selected_features': regression_results['num_selected_features']
        }
        overall_metrics_list.append(current_metrics)
        
        sig_time = (time.time() - signature_process_start_time) / 60
        print(f"Time taken for signature {rna_signature_name}: {sig_time:.2f} minutes.")

    # Save overall metrics summary
    if overall_metrics_list:
        summary_df = pd.DataFrame(overall_metrics_list)
        summary_filename = os.path.join(base_output_dir, "regression_model_training_summary.csv")
        summary_df.to_csv(summary_filename, index=False)
        print(f"\nOverall regression model training summary saved to {summary_filename}")
        print("--- Overall Regression Summary ---")
        print(summary_df)
        print("---------------------------------")
    else:
        print("\nNo signatures were processed or no valid results obtained from regression training.")

    total_execution_time = (time.time() - overall_start_time) / 60
    print(f"\nTotal script execution time: {total_execution_time:.2f} minutes.")

if __name__ == '__main__':
    main() 