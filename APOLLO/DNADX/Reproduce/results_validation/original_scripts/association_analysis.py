import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import ElasticNet
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from joblib import Parallel, delayed, dump, load
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("results_validation/score/signature_score.pkl")
    segment_score = pd.read_pickle("results_validation/score/segment_score.pkl")
    
    # Find common samples using exact names
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

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
#    results_df.to_csv('results_validation/association_results.csv', index=False)
#    print(f"\nSaved {len(results_df)} association results")
#    
#    # Print most significant correlations
#    print("\nTop 10 most significant correlations:")
#    print(results_df.head(10))
#    
#    return results_df

def get_lambda_range(X_df, y_series, current_l1_ratio, n_strengths=20, eps=0.01):
    """
    Determine a range of alpha (strength) values for ElasticNet for a given l1_ratio.
    The range is determined similarly to how scikit-learn's ElasticNetCV would find its alpha path.
    """
    X = X_df.values if isinstance(X_df, (pd.DataFrame, pd.Series)) else X_df
    y = y_series.values if isinstance(y_series, (pd.DataFrame, pd.Series)) else y_series
    n_samples = X.shape[0]

    if current_l1_ratio < 1e-8: # l1_ratio is practically zero (Ridge-like)
        # The formula for strength_max below is primarily for l1_ratio > 0.
        # For a pure Ridge or near-Ridge, provide a generic wide range of strengths.
        # print(f"Warning: l1_ratio ({current_l1_ratio}) is close to zero. Using a default wide strength range.")
        return np.logspace(-5, 3, n_strengths)

    # For ElasticNet (l1_ratio > 0):
    # Calculate strength_max = max(|X^T (y - y_mean)|) / (n_samples * l1_ratio)
    # This is the smallest strength (alpha in scikit-learn ElasticNet) for which all coefficients are zero.
    y_mean = np.mean(y)
    # Ensure X is 2D for dot product, even if it has only one feature
    X_2d = X.reshape(n_samples, -1) if X.ndim == 1 else X
    X_T_y_centered = np.dot(X_2d.T, y - y_mean)
    
    if X_T_y_centered.ndim == 0: # If X had only 1 feature, X_T_y_centered might be scalar
        strength_max_numerator = np.abs(X_T_y_centered)
    else:
        strength_max_numerator = np.max(np.abs(X_T_y_centered))

    if n_samples * current_l1_ratio == 0: # Avoid division by zero
         # print(f"Warning: n_samples * current_l1_ratio is zero. Using default strength range for l1_ratio {current_l1_ratio}.")
         return np.logspace(-4, 1, n_strengths)

    strength_max = strength_max_numerator / (n_samples * current_l1_ratio)

    if strength_max <= 1e-9: 
        # This might happen if y is constant or X has no correlation with y, or X is all zero.
        # print(f"Warning: Calculated strength_max is very small or zero for l1_ratio {current_l1_ratio}. Using default range.")
        return np.logspace(-4, 1, n_strengths) 

    strength_min = strength_max * eps
    
    # Ensure strength_min and strength_max are sensible and min < max
    if strength_min <= 1e-9 : strength_min = 1e-7 * strength_max 
    if strength_min == 0 : strength_min = 1e-9 
    if strength_min >= strength_max : strength_min = strength_max / (float(n_strengths * 10) if n_strengths > 0 else 10.0)
    if strength_min <= 1e-9 : strength_min = 1e-9

    if strength_max <= strength_min : 
         strength_max = strength_min * 100.0
         if strength_max <= 1e-9: strength_max = 0.01 
         if strength_min >= strength_max : strength_min = strength_max / 100.0 # re-adjust min if max was floored

    # Final check to prevent log10 of zero or negative
    strength_min = max(strength_min, 1e-9)
    strength_max = max(strength_max, strength_min * 1.1) # ensure max > min

    return np.logspace(np.log10(strength_min), np.log10(strength_max), n_strengths)

def monte_carlo_cv(X, y, l1_ratio_val, strength_val, n_splits=50, train_size=0.75):
    """Perform Monte Carlo cross validation"""
    accuracies = []
    # Ensure X and y are numpy arrays for scikit-learn
    X_np = X if isinstance(X, np.ndarray) else X.values
    y_np = y if isinstance(y, np.ndarray) else y.values

    for _ in range(n_splits):
        # Split data
        # Using fixed random_state=42 for reproducibility within monte_carlo_cv's splits
        # The outer split in build_prediction_model already has random_state=42
        X_train_mc, X_test_mc, y_train_mc, y_test_mc = train_test_split(
            X_np, y_np, train_size=train_size, random_state=42 # Using fixed random state for MCCV splits
        )
        
        # Train model - Corrected alpha and l1_ratio usage
        model = ElasticNet(alpha=strength_val, l1_ratio=l1_ratio_val, max_iter=10000, tol=1e-4) # Increased tol for faster CV
        model.fit(X_train_mc, y_train_mc)
        
        # Make predictions and calculate accuracy
        y_pred_mc = model.predict(X_test_mc)
        
        # Binarize based on medians. If y is scaled, median will be ~0.
        y_test_mc_median = np.median(y_test_mc)
        y_pred_mc_median = np.median(y_pred_mc)

        y_test_bin = (y_test_mc > y_test_mc_median).astype(int)
        y_pred_bin = (y_pred_mc > y_pred_mc_median).astype(int)
        
        accuracy = accuracy_score(y_test_bin, y_pred_bin)
        accuracies.append(accuracy)
    
    return np.mean(accuracies)

def evaluate_params(args):
    """Helper function for parallel parameter evaluation"""
    # Parameters are X_train_scaled, y_train_scaled, l1_ratio, strength
    X, y, l1_ratio, strength = args 
    return (l1_ratio, strength, monte_carlo_cv(X, y, l1_ratio, strength))

def impute_preserve_all_columns(X):
    """Impute each column: if all-NaN, fill with 0; else fill NaN with mean"""
    X_imputed = X.copy()
    for col in X_imputed.columns:
        if X_imputed[col].isna().all():
            X_imputed[col] = 0
        else:
            mean_val = X_imputed[col].mean()
            X_imputed[col] = X_imputed[col].fillna(mean_val)
    return X_imputed

def print_descriptive_stats(data, name, is_feature_matrix=False):
    """Prints descriptive statistics for a given dataset (features or target)."""
    print(f"Descriptive statistics for {name}:")
    if data is None:
        print("  Data is None.")
        return

    if is_feature_matrix:
        X_np = data if isinstance(data, np.ndarray) else data.values
        if X_np.ndim == 1: X_np = X_np.reshape(-1, 1) # Ensure 2D for feature matrix logic

        if X_np.shape[1] == 0:
            print(f"  {name} has 0 features.")
            return
            
        num_features_to_show = min(5, X_np.shape[1])
        print(f"  (Showing for first {num_features_to_show} of {X_np.shape[1]} features if more than 1 feature)")
        
        for i in range(num_features_to_show):
            col_data = X_np[:, i]
            col_finite = col_data[np.isfinite(col_data)]
            if col_finite.size == 0:
                print(f"    Feature {i+1}: All NaN/inf or empty.")
                continue
            
            mode_res = stats.mode(col_finite)
            actual_mode = mode_res.mode[0] if isinstance(mode_res.mode, np.ndarray) and mode_res.mode.size > 0 else mode_res.mode
            actual_count = mode_res.count[0] if isinstance(mode_res.count, np.ndarray) and mode_res.count.size > 0 else mode_res.count
            
            print(f"    Feature {i+1}: Min={np.min(col_finite):.3f}, Max={np.max(col_finite):.3f}, Mean={np.mean(col_finite):.3f}, Median={np.median(col_finite):.3f}, Mode={actual_mode:.3f} (Count: {actual_count})")
        if X_np.shape[1] > num_features_to_show:
            print(f"    ... (stats for {X_np.shape[1] - num_features_to_show} more features not shown)")
    else: # Target variable
        y_np = data if isinstance(data, np.ndarray) else data.values
        y_finite = y_np[np.isfinite(y_np)]

        if y_finite.size == 0:
            print(f"  {name} (target): All NaN/inf or empty.")
            return

        mode_res_y = stats.mode(y_finite)
        actual_mode_y = mode_res_y.mode[0] if isinstance(mode_res_y.mode, np.ndarray) and mode_res_y.mode.size > 0 else mode_res_y.mode
        actual_count_y = mode_res_y.count[0] if isinstance(mode_res_y.count, np.ndarray) and mode_res_y.count.size > 0 else mode_res_y.count
        print(f"  {name} (target): Min={np.min(y_finite):.3f}, Max={np.max(y_finite):.3f}, Mean={np.mean(y_finite):.3f}, Median={np.median(y_finite):.3f}, Mode={actual_mode_y:.3f} (Count: {actual_count_y})")

def build_prediction_model(signature_score, segment_score, target_signature):
    """Build prediction model for a target signature using segment scores"""
    print(f"\nBuilding prediction model for {target_signature}...")
    
    # Prepare data
    X = segment_score.T.astype(float)
    y = signature_score.loc[target_signature].astype(float)
    
    original_feature_names = X.columns.tolist() # Save original feature names
    original_X_median = X.median().to_dict() # Save original medians

    print(f"Initial data shape - X: {X.shape}, y: {y.shape}")
    
    # Check for NaN and infinite values
    print(f"NaN values in X: {np.isnan(X).sum().sum()}")
    print(f"NaN values in y: {np.isnan(y).sum()}")
    print(f"Infinite values in X: {np.isinf(X).sum().sum()}")
    print(f"Infinite values in y: {np.isinf(y).sum()}")
    
    # Impute missing values using the same strategy as training/validation
    X = impute_preserve_all_columns(X)
    
    # Keep only samples where y has finite values
    mask = np.isfinite(y)
    X = X[mask]
    y = y[mask]
    
    print(f"Data shape after filtering - X: {X.shape}, y: {y.shape}")
    
    if len(X) < 10:  # Require at least 10 samples
        print(f"Not enough valid samples after filtering (need at least 10, got {len(X)})")
        return None, None, None, None
    
    # Split data - using 70/30 split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    print(f"Training data shape (before scaling) - X_train: {X_train.shape}, y_train: {y_train.shape}")
    print(f"Test data shape (before scaling) - X_test: {X_test.shape}, y_test: {y_test.shape}")

    # --- Print stats before scaling ---
    print_descriptive_stats(X_train, "X_train (CNA scores) before scaling", is_feature_matrix=True)
    print_descriptive_stats(y_train, "y_train (RNA score) before scaling")

    # --- Scaling Data ---
    print("\nScaling data...")
    scaler_X = StandardScaler()
    X_train_scaled = scaler_X.fit_transform(X_train) # X_train is DataFrame, output is numpy
    X_test_scaled = scaler_X.transform(X_test)     # X_test is DataFrame, output is numpy

    scaler_y = StandardScaler()
    y_train_scaled = scaler_y.fit_transform(y_train.values.reshape(-1, 1)).ravel() # y_train is Series, output is numpy
    y_test_scaled = scaler_y.transform(y_test.values.reshape(-1, 1)).ravel()     # y_test is Series, output is numpy
    
    print("Data scaling complete.")
    print(f"Scaled training data shape - X_train_scaled: {X_train_scaled.shape}, y_train_scaled: {y_train_scaled.shape}")

    # --- Print stats after scaling ---
    print_descriptive_stats(X_train_scaled, "X_train_scaled (CNA scores) after scaling", is_feature_matrix=True)
    print_descriptive_stats(y_train_scaled, "y_train_scaled (RNA score) after scaling")

    # Parameter tuning (using scaled data)
    print("\nPerforming parameter tuning (using scaled data)...")
    # The 'alphas' here are L1_RATIOS
    l1_ratios_to_try = np.array([0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]) # Expanded l1_ratios, incl. near Ridge and near Lasso
    
    param_combinations = []
    # Pass scaled data to get_lambda_range and then to evaluate_params
    for l1_ratio_val in l1_ratios_to_try:
        # get_lambda_range expects DataFrame/Series for X,y if it uses .values, or numpy arrays
        # X_train_scaled and y_train_scaled are numpy arrays.
        # The updated get_lambda_range handles numpy arrays directly.
        strengths_for_l1 = get_lambda_range(X_train_scaled, y_train_scaled, l1_ratio_val)
        for strength_val in strengths_for_l1:
            param_combinations.append((X_train_scaled, y_train_scaled, l1_ratio_val, strength_val))
    
    if not param_combinations:
        print("Error: No parameter combinations generated for tuning. Skipping model for this signature.")
        return None, None, None, None

    # Run parameter tuning in parallel
    results = Parallel(n_jobs=-1)(delayed(evaluate_params)(params) for params in param_combinations)
    
    # Find best parameters
    best_accuracy = -np.inf
    best_params = None # Should store (best_l1_ratio, best_strength)
    for l1_ratio_res, strength_res, accuracy_res in results:
        if accuracy_res > best_accuracy:
            best_accuracy = accuracy_res
            best_params = (l1_ratio_res, strength_res)
    
    if best_params is None:
        print("Error: Could not find best parameters after tuning. Skipping model for this signature.")
        return None, None, None, None
        
    print(f"Best parameters found - L1 Ratio: {best_params[0]:.4f}, Strength (Alpha): {best_params[1]:.6f}, Best CV Accuracy: {best_accuracy:.4f}")
    
    # Train final model with best parameters on scaled data
    # Corrected: alpha is strength (best_params[1]), l1_ratio is mix (best_params[0])
    final_model = ElasticNet(alpha=best_params[1], l1_ratio=best_params[0], max_iter=10000, tol=1e-4) # Increased tol
    final_model.fit(X_train_scaled, y_train_scaled)
    
    # Save the model and preprocessing information
    # Feature names and median values should be from the *original* unscaled X
    model_info = {
        'model': final_model,
        'scaler_X': scaler_X, # Save scalers for potential future use
        'scaler_y': scaler_y,
        'feature_names': original_feature_names, 
        'median_values': original_X_median, 
        'best_params': {'l1_ratio': best_params[0], 'alpha_strength': best_params[1]},
        'cv_accuracy': best_accuracy
    }
    
    # Create directory if it doesn't exist
    import os
    os.makedirs('results_validation/models', exist_ok=True)
    
    # Save the model and preprocessing information
    dump(model_info, f'results_validation/models/{target_signature}_model.joblib')
    print(f"Saved model and preprocessing information to results_validation/models/{target_signature}_model.joblib")
    
    # Make predictions on scaled test data
    y_pred_scaled = final_model.predict(X_test_scaled)
    
    # Inverse transform predictions to original scale for metrics
    y_pred_orig_scale = scaler_y.inverse_transform(y_pred_scaled.reshape(-1, 1)).ravel()

    # Calculate performance metrics using original scale y_test and inverse_transformed y_pred
    # y_test is the original unscaled test set target values.
    correlation = np.corrcoef(y_test, y_pred_orig_scale)[0, 1] if len(y_test) > 1 and len(y_pred_orig_scale) > 1 else np.nan
    mse = np.mean((y_test - y_pred_orig_scale) ** 2) if len(y_test) > 0 and len(y_pred_orig_scale) > 0 else np.nan
    
    # For ROC/AUC, binarize y_test (original scale) using its median as threshold
    # Use y_pred_orig_scale for ROC curve as it's comparable to y_test's scale
    if len(y_test) > 1 and not np.all(np.isnan(y_test)) and not np.all(np.isnan(y_pred_orig_scale)):
        y_test_median_orig = np.median(y_test) 
        y_test_bin = (y_test > y_test_median_orig).astype(int)
        
        # Check if y_test_bin has both classes, otherwise AUC is not well-defined
        if len(np.unique(y_test_bin)) > 1:
            fpr, tpr, thresholds = roc_curve(y_test_bin, y_pred_orig_scale)
            auc_score = auc(fpr, tpr)
        else:
            print(f"Warning: y_test_bin for {target_signature} has only one class after binarization. AUC cannot be computed.")
            auc_score = np.nan
            fpr, tpr = np.array([]), np.array([]) # Empty arrays for plotting
    else:
        auc_score = np.nan
        fpr, tpr = np.array([]), np.array([])

    # Plot results (original scale)
    plt.figure(figsize=(10, 6))
    plt.scatter(y_test, y_pred_orig_scale, alpha=0.5)
    plt.plot([min(y_test.min(), y_pred_orig_scale.min()) if y_test.size > 0 and y_pred_orig_scale.size > 0 else 0, 
              max(y_test.max(), y_pred_orig_scale.max()) if y_test.size > 0 and y_pred_orig_scale.size > 0 else 1],
             [min(y_test.min(), y_pred_orig_scale.min()) if y_test.size > 0 and y_pred_orig_scale.size > 0 else 0, 
              max(y_test.max(), y_pred_orig_scale.max()) if y_test.size > 0 and y_pred_orig_scale.size > 0 else 1], 'r--')
    plt.xlabel('True Values (Original Scale)')
    plt.ylabel('Predictions (Original Scale)')
    plt.title(f'Prediction Performance for {target_signature}\nCorrelation: {correlation:.3f}, MSE: {mse:.3f}, AUC: {auc_score:.3f}')
    plt.savefig(f'results_validation/{target_signature}_prediction.png', dpi=300, bbox_inches='tight')
    plt.close() # Close plot to free memory

    # Return metrics (can be extended)
    return correlation, mse, auc_score, model_info

def main():
    print("Starting association analysis and model building...")
    signature_score, segment_score = load_data()

    if signature_score.empty or segment_score.empty:
        print("Data loading failed or resulted in empty dataframes. Exiting.")
        return

    # Example: Build model for a subset of signatures or all
    # target_signatures = signature_score.index[:3] # For testing with first 3 signatures
    target_signatures = signature_score.index

    results_summary = []
    # Enumerate for progress tracking
    for i, target_signature in enumerate(target_signatures):
        print(f"\n--- Processing signature {i+1}/{len(target_signatures)}: {target_signature} ---")
        try:
            correlation, mse, auc_score, model_info = build_prediction_model(signature_score, segment_score, target_signature)
            if correlation is not None: # Check if model building was successful
                print(f"Results for {target_signature}: Correlation={correlation:.4f}, MSE={mse:.4f}, AUC={auc_score:.4f}")
                results_summary.append({
                    'signature': target_signature,
                    'correlation': correlation,
                    'mse': mse,
                    'auc': auc_score,
                    'l1_ratio': model_info['best_params']['l1_ratio'] if model_info else np.nan,
                    'alpha_strength': model_info['best_params']['alpha_strength'] if model_info else np.nan,
                    'cv_accuracy': model_info['cv_accuracy'] if model_info else np.nan
                })
            else:
                print(f"Skipped model building for {target_signature} due to issues (e.g., not enough data).")
                results_summary.append({
                    'signature': target_signature,
                    'correlation': np.nan, 'mse': np.nan, 'auc': np.nan, 
                    'l1_ratio': np.nan, 'alpha_strength': np.nan, 'cv_accuracy': np.nan
                })
        except Exception as e:
            print(f"!!! Critical Error building model for signature {target_signature}: {e}")
            import traceback
            traceback.print_exc() # Print full traceback for debugging
            results_summary.append({
                'signature': target_signature,
                'correlation': np.nan, 'mse': np.nan, 'auc': np.nan,
                'l1_ratio': np.nan, 'alpha_strength': np.nan, 'cv_accuracy': np.nan,
                'error': str(e)
            })
            # Continue to the next signature

    # Save summary results
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        summary_df.to_csv("results_validation/model_metrics_summary.csv", index=False)
        print("\nSaved summary of model metrics to results_validation/model_metrics_summary.csv")
        print(summary_df.head())
    else:
        print("\nNo results to summarize.")
    
    print("\nAnalysis complete.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "elasticnet_only":
        main()
    else:
        main() 