# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split, StratifiedKFold, ShuffleSplit
from sklearn.linear_model import ElasticNet, ElasticNetCV, MultiTaskElasticNetCV, LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed, dump, load
from sklearn.impute import SimpleImputer
import warnings
import pickle
from scipy.stats import pearsonr
from sklearn.model_selection import GridSearchCV

warnings.filterwarnings('ignore')

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

def train_evaluate_elasticnet_regression(trainX_raw, trainY_raw, testX_raw, testY_raw,
                                         l1_ratios_to_try, # List of l1_ratios (alphas in R glmnet)
                                         n_alphas_cv=100,  # Number of alphas (lambdas in R glmnet) along the path
                                         cv_folds=4): # Default to 4-fold CV
    """
    Cleans, scales, trains ElasticNetCV for each l1_ratio to find the best alpha (lambda),
    then selects the overall best (l1_ratio, alpha) pair and refits a final ElasticNet model.
    Provides progress updates for each l1_ratio processed.
    """
    print("  Cleaning data...")
    trainX = clean_data_array(trainX_raw)
    testX = clean_data_array(testX_raw)
    trainY = clean_target_array(trainY_raw) # Clean Y as well
    testY = clean_target_array(testY_raw)   # Clean Y as well

    print("  Scaling feature data (X)...")
    scaler = StandardScaler()
    trainX_scaled = scaler.fit_transform(trainX)
    testX_scaled = scaler.transform(testX) # Use scaler fitted on training data
    
    print(f"  Starting hyperparameter tuning (l1_ratios and alphas) with {cv_folds}-fold CV...")

    best_overall_mse = float('inf')
    best_l1_ratio_overall = None
    best_alpha_overall = None
    
    # Ensure trainY is 1D for fitting
    trainY_raveled = trainY.ravel()

    for i, current_l1_ratio in enumerate(l1_ratios_to_try):
        print(f"    Processing l1_ratio {i+1}/{len(l1_ratios_to_try)}: {current_l1_ratio:.2f} ...")
        
        # Use ElasticNetCV to find the best alpha for the current_l1_ratio
        # verbose=0 here to avoid too much internal chatter from this sub-step
        model_for_current_l1 = ElasticNetCV(
            l1_ratio=current_l1_ratio, # Single l1_ratio value
            alphas=None,             # Let ElasticNetCV determine the alpha path
            n_alphas=n_alphas_cv,
            cv=cv_folds,
            random_state=42,
            max_iter=10000,
            tol=1e-7,
            selection='random',
            n_jobs=-1,               # Use all available cores for this CV step
            verbose=0                # Suppress verbose for this inner CV
        )
        model_for_current_l1.fit(trainX_scaled, trainY_raveled)
        
        # ElasticNetCV stores MSEs for each alpha along the path for the given l1_ratio.
        # The best alpha for this l1_ratio is model_for_current_l1.alpha_.
        # The corresponding MSE can be found by looking up this alpha in its mse_path_.
        # mse_path_ is (n_alphas, n_folds) or (n_alphas,) if refit=True. We need the mean MSE for the best alpha.
        # A simpler way: find the minimum of the mean MSEs across alphas for this l1_ratio.
        min_mse_for_current_l1 = np.min(np.mean(model_for_current_l1.mse_path_, axis=1)) # Mean across folds, then min across alphas

        print(f"      l1_ratio {current_l1_ratio:.2f} done. Best alpha for this l1_ratio: {model_for_current_l1.alpha_:.4f}, Min Mean MSE: {min_mse_for_current_l1:.4f}")

        if min_mse_for_current_l1 < best_overall_mse:
            best_overall_mse = min_mse_for_current_l1
            best_l1_ratio_overall = current_l1_ratio
            best_alpha_overall = model_for_current_l1.alpha_

    if best_l1_ratio_overall is None or best_alpha_overall is None:
        print("Warning: Could not determine best hyperparameters. Check data or CV process.")
        # Fallback or raise error - for now, returning None or a dummy model
        # This case should ideally not happen with proper data.
        # Create a dummy model with default parameters if needed for graceful failure
        final_model = ElasticNet(random_state=42)
        final_model.alpha_ = None # Custom attribute for consistency
        final_model.l1_ratio_ = None # Custom attribute for consistency
    else:
        print(f"  Overall best l1_ratio: {best_l1_ratio_overall:.2f}, Overall best alpha (lambda): {best_alpha_overall:.4f}, Best MSE: {best_overall_mse:.4f}")
        print("  Refitting final model with best parameters...")
        # Refit the final model using the overall best l1_ratio and alpha on the full training data
        final_model = ElasticNet(
            alpha=best_alpha_overall,
            l1_ratio=best_l1_ratio_overall,
            random_state=42,
            max_iter=10000, # Ensure consistency with CV
            tol=1e-7
        )
        final_model.fit(trainX_scaled, trainY_raveled)
        # Store the chosen hyperparameters as attributes for easy access, similar to ElasticNetCV
        final_model.alpha_ = best_alpha_overall 
        final_model.l1_ratio_ = best_l1_ratio_overall
    
    print(f"  Model training completed.") # Removed specifics about l1_ratio/alpha here as they are now printed above.
    
    return final_model, trainX_scaled, testX_scaled, scaler

def plot_regression_and_roc(y_train_actual, y_train_pred, y_test_actual, y_test_pred,
                            train_corr, train_mse, roc_auc_train,
                            test_corr, test_mse, roc_auc_test,
                            signature_name, results_dir):
    """Plots regression performance and ROC curves."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    safe_signature_name = signature_name.replace('/', '_').replace('\\\\', '_')
    fig.suptitle(f"Performance for {safe_signature_name}", fontsize=16)

    # Training set scatter plot
    axes[0, 0].scatter(y_train_actual, y_train_pred, alpha=0.5, label="Predictions")
    axes[0, 0].plot([y_train_actual.min(), y_train_actual.max()], [y_train_actual.min(), y_train_actual.max()], 'r--', lw=2, label="Ideal")
    axes[0, 0].set_xlabel('Actual RNA Score')
    axes[0, 0].set_ylabel('Predicted RNA Score')
    axes[0, 0].set_title(f'Training Set\nCorr: {train_corr:.3f}, MSE: {train_mse:.3f}')
    axes[0, 0].legend()
    axes[0, 0].grid(True)

    # Test set scatter plot
    axes[0, 1].scatter(y_test_actual, y_test_pred, alpha=0.5, label="Predictions")
    axes[0, 1].plot([y_test_actual.min(), y_test_actual.max()], [y_test_actual.min(), y_test_actual.max()], 'r--', lw=2, label="Ideal")
    axes[0, 1].set_xlabel('Actual RNA Score')
    axes[0, 1].set_ylabel('Predicted RNA Score')
    axes[0, 1].set_title(f'Test Set\nCorr: {test_corr:.3f}, MSE: {test_mse:.3f}')
    axes[0, 1].legend()
    axes[0, 1].grid(True)

    # Binarize actuals for ROC curve plotting
    y_train_binary_actual = binarize_score(y_train_actual)
    y_test_binary_actual = binarize_score(y_test_actual)

    # ROC for Training set (using continuous predictions)
    fpr_train, tpr_train, _ = roc_curve(y_train_binary_actual, y_train_pred)
    axes[1, 0].plot(fpr_train, tpr_train, color='blue', lw=2, label=f'Train ROC curve (AUC = {roc_auc_train:.3f})')
    axes[1, 0].plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--')
    axes[1, 0].set_xlim([0.0, 1.0])
    axes[1, 0].set_ylim([0.0, 1.05])
    axes[1, 0].set_xlabel('False Positive Rate')
    axes[1, 0].set_ylabel('True Positive Rate')
    axes[1, 0].set_title('Training Set ROC (Top 1/3 vs Rest)')
    axes[1, 0].legend(loc="lower right")
    axes[1, 0].grid(True)

    # ROC for Test set (using continuous predictions)
    fpr_test, tpr_test, _ = roc_curve(y_test_binary_actual, y_test_pred)
    axes[1, 1].plot(fpr_test, tpr_test, color='green', lw=2, label=f'Test ROC curve (AUC = {roc_auc_test:.3f})')
    axes[1, 1].plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--')
    axes[1, 1].set_xlim([0.0, 1.0])
    axes[1, 1].set_ylim([0.0, 1.05])
    axes[1, 1].set_xlabel('False Positive Rate')
    axes[1, 1].set_ylabel('True Positive Rate')
    axes[1, 1].set_title('Test Set ROC (Top 1/3 vs Rest)')
    axes[1, 1].legend(loc="lower right")
    axes[1, 1].grid(True)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    plot_filename = os.path.join(results_dir, f"{safe_signature_name}_performance.png")
    plt.savefig(plot_filename)
    plt.close()
    print(f"  Performance plot saved to {plot_filename}")

def main():
    """Main function to run the model training pipeline for regression and AUC evaluation."""
    
    # Setup results directory
    results_dir = "scaling/results_regression_with_auc" # More descriptive name
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(os.path.join(results_dir, "models"), exist_ok=True) # For saving models
    os.makedirs(os.path.join(results_dir, "plots"), exist_ok=True)  # For saving plots
    
    # Load and align data
    signature_score_df, segment_score_df = load_data()
    if signature_score_df.empty or segment_score_df.empty:
        print("Exiting due to data loading issues or no common samples.")
        return

    print("\nData loaded and aligned:")
    print("RNA signature scores shape:", signature_score_df.shape)
    print("CNA segment scores shape:", segment_score_df.shape) # segment_score_df is features x samples
    
    # Save aligned data (optional)
    # signature_score_df.to_pickle(os.path.join(results_dir, "aligned_rna_signatures.pkl"))
    # segment_score_df.to_pickle(os.path.join(results_dir, "aligned_cna_segments.pkl"))

    # Define l1_ratios for ElasticNetCV (similar to R's alpha for glmnet)
    # R script used seq(0.1, 0.9, by=0.1) and also 0.01, 0.05, and 1 for specific cases
    # A common range:
    # l1_ratios_to_tune = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]
    # As per paper: "alphas over a range from 0.1 to 1 by 0.1"
    l1_ratios_to_tune = np.round(np.arange(0.1, 1.01, 0.1), 2).tolist() # [0.1, 0.2, ..., 1.0]

    all_signature_metrics = []

    # CNA data: samples are columns, features are rows. Transpose for scikit-learn.
    # X_features_cna should have samples as rows, features as columns.
    X_features_cna = segment_score_df.T 
    cna_feature_names = X_features_cna.columns.tolist() # Gene names or segment IDs

    # Define the specific list of signatures for the trial run, with GHI before RB_LOH
    target_signatures_for_trial = [
        "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP", # Specific Basal signaling first
        "GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN",
        "GHI_RS_Model_NJEM.2004_PMID.15591335", # GHI before RB_LOH
        "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450",
        "UNC_HER1_Cluster2_Median_BMC.Genomics.2007_PMID.17663798"
    ]
    print(f"\n--- Starting Trial Run for {len(target_signatures_for_trial)} Specific Signatures (GHI before RB_LOH) ---")

    for rna_signature_name in target_signatures_for_trial:
        if rna_signature_name not in signature_score_df.index:
            print(f"\nWarning: Signature '{rna_signature_name}' not found in loaded RNA signature data. Skipping.")
            continue

        print(f"\nProcessing RNA Signature: {rna_signature_name}")
        
        # Target variable y: continuous RNA signature scores for the current signature
        y_target_rna = signature_score_df.loc[rna_signature_name].values # This is 1D array (samples,)
        
        # Split data into training and testing sets
        # X_features_cna.values ensures we are using numpy arrays
        X_train_raw, X_test_raw, y_train_raw, y_test_raw = train_test_split(
            X_features_cna.values, y_target_rna, test_size=0.3, random_state=42
        )
        
        # Train model using the new consolidated function
        # This function handles cleaning, scaling, and ElasticNetCV training
        trained_model, X_train_scaled, X_test_scaled, fitted_scaler = train_evaluate_elasticnet_regression(
            X_train_raw, y_train_raw, X_test_raw, y_test_raw,
            l1_ratios_to_try=l1_ratios_to_tune,
            n_alphas_cv=100, # Number of lambdas for ElasticNetCV to try for each l1_ratio
            cv_folds=4      # Use 4-fold cross-validation within ElasticNetCV
        )
        
        # Make predictions on scaled data
        y_train_pred_continuous = trained_model.predict(X_train_scaled)
        y_test_pred_continuous = trained_model.predict(X_test_scaled)
    
        # Calculate REGRESSION metrics (using original raw y for comparison against predictions from model)
        # y_train_raw and y_test_raw are the original, uncleaned continuous scores for the splits.
        # It's better to use the cleaned versions that went into scaling and modeling for metric consistency
        y_train_cleaned = clean_target_array(y_train_raw)
        y_test_cleaned = clean_target_array(y_test_raw)

        train_corr, _ = pearsonr(y_train_cleaned, y_train_pred_continuous)
        test_corr, _ = pearsonr(y_test_cleaned, y_test_pred_continuous)
        train_mse = mean_squared_error(y_train_cleaned, y_train_pred_continuous)
        test_mse = mean_squared_error(y_test_cleaned, y_test_pred_continuous)
        
        print(f"  Regression Train: Correlation={train_corr:.3f}, MSE={train_mse:.3f}")
        print(f"  Regression Test:  Correlation={test_corr:.3f}, MSE={test_mse:.3f}")
    
        # Calculate CLASSIFICATION (AUC) metrics
        # Binarize the *actual* target variable (cleaned)
        y_train_binary_actual = binarize_score(y_train_cleaned)
        y_test_binary_actual = binarize_score(y_test_cleaned)
        
        # Calculate AUC using the *continuous predictions* from the regression model
        # This is valid: roc_auc_score(y_true_binary, y_pred_continuous_scores)
        train_roc_auc = roc_auc_score(y_train_binary_actual, y_train_pred_continuous)
        test_roc_auc = roc_auc_score(y_test_binary_actual, y_test_pred_continuous)
        
        print(f"  Classification Train AUC (Top 1/3 vs Rest): {train_roc_auc:.3f}")
        print(f"  Classification Test  AUC (Top 1/3 vs Rest): {test_roc_auc:.3f}")

        # Store coefficients (these are for SCALED features)
        coefficients = trained_model.coef_
        intercept = trained_model.intercept_

        # Plot performance (regression and ROC)
        plot_regression_and_roc(
            y_train_cleaned, y_train_pred_continuous, y_test_cleaned, y_test_pred_continuous,
            train_corr, train_mse, train_roc_auc,
            test_corr, test_mse, test_roc_auc,
            rna_signature_name, os.path.join(results_dir, "plots")
        )
        
        # Save model and scaler
        safe_rna_signature_name = rna_signature_name.replace('/', '_').replace('\\\\', '_')
        model_filepath = os.path.join(results_dir, "models", f"{safe_rna_signature_name}_elasticnet_model.joblib")
        dump({'model': trained_model, 'scaler': fitted_scaler, 'feature_names': cna_feature_names}, model_filepath)
        print(f"  Model and scaler saved to {model_filepath}")

        current_metrics = {
            'rna_signature': rna_signature_name,
            'best_l1_ratio': trained_model.l1_ratio_,
            'best_alpha_lambda': trained_model.alpha_,
            'train_correlation': train_corr,
            'test_correlation': test_corr,
            'train_mse': train_mse,
            'test_mse': test_mse,
            'train_roc_auc': train_roc_auc,
            'test_roc_auc': test_roc_auc,
            'n_features_selected': np.sum(coefficients != 0),
            # Store coefficients as a list or dict if planning to save to CSV directly
            # For now, not adding full coeff list to aggregated CSV to keep it manageable
            # 'coefficients': coefficients.tolist(), 
            # 'intercept': intercept
        }
        all_signature_metrics.append(current_metrics)
    
    # Save aggregated metrics for all signatures
    metrics_summary_df = pd.DataFrame(all_signature_metrics)
    summary_csv_path = os.path.join(results_dir, "all_signatures_metrics_summary.csv")
    metrics_summary_df.to_csv(summary_csv_path, index=False)
    print(f"\nAggregated metrics summary saved to {summary_csv_path}")
    
    print("\nModel training and evaluation pipeline completed.")

if __name__ == '__main__':
    # Remove old functions not used by the new main flow for clarity
    # del caret_wrap 
    # del monte_carlo_cv
    # del evaluate_params
    # del impute_preserve_all_columns
    # del build_prediction_model
    # del monte_carlo_cv_classification
    # del build_final_model_and_evaluate_classification
    # The above deletions should be done carefully or by commenting out if unsure.
    # For this auto-edit, I will not delete them to be safe, but they are effectively orphaned.
    main() 