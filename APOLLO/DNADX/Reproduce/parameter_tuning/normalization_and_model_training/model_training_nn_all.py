# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, roc_curve, auc
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shutil # Added for directory cleaning
from scipy.stats import pearsonr
import time
import joblib
import torch
from tqdm import tqdm # Added for progress bar in main

def sanitize_for_file_path(name):
    """Sanitize a string to be a valid file path component, preserving ., =, _, and spaces,
       but replacing characters problematic for file/directory names like / \ :"""
    name = name.replace('/', '__SLASH__')
    name = name.replace('\\', '__BACKSLASH__') # Note: single backslash is escape in Python string
    name = name.replace(':', '__COLON__')
    # Further characters to consider if issues arise on specific OS:
    # name = name.replace('*', '__STAR__')
    # name = name.replace('?', '__QUESTION__')
    # name = name.replace('"', '__QUOTE__')
    # name = name.replace('<', '__LT__')
    # name = name.replace('>', '__GT__')
    # name = name.replace('|', '__PIPE__')
    return name

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("training_data/rna_signature_score_median_no_norm.pkl")
    segment_score = pd.read_pickle("training_data/cna_segment_score_mean_no_norm.pkl")
    
    print("\nSample name formats:")
    print("RNA signature score sample names (first 5):", list(signature_score.columns[:5]))
    print("CNA segment score sample names (first 5):", list(segment_score.columns[:5]))
    
    # Convert CNA sample names to match RNA format (replace . with -)
    segment_score.columns = segment_score.columns.str.replace('.', '-', regex=False) # Added regex=False for clarity
    
    # Find common samples using exact names
    common_samples = sorted(set(signature_score.columns) & set(segment_score.columns))
    print("\nFound {} common samples".format(len(common_samples)))
    if len(common_samples) > 0:
        print("First 5 common samples:", common_samples[:5])
    
    # Align by common samples
    signature_score = signature_score[common_samples]
    segment_score = segment_score[common_samples]
    
    return signature_score, segment_score

def clean_data(X):
    """Clean data by handling NaN values"""
    if isinstance(X, pd.DataFrame):
        X = X.values
    
    # Convert to float type first
    X = X.astype(float)
    
    # Replace inf values with NaN
    X = np.where(np.isinf(X), np.nan, X)
    
    # Impute NaN values with column means
    imputer = SimpleImputer(strategy='mean')
    X_cleaned = imputer.fit_transform(X)
    
    return X_cleaned

def create_binary_labels(y_continuous, threshold_type='tercile'):
    """Convert continuous values to binary labels based on threshold.
       Handles potential NaNs in y_continuous.
    """
    if not isinstance(y_continuous, pd.Series):
        y_continuous = pd.Series(y_continuous)

    if y_continuous.empty or y_continuous.isnull().all():
        # print("    Warning: y_continuous is empty or all NaN. Cannot determine threshold.")
        return pd.Series([np.nan] * len(y_continuous), index=y_continuous.index), np.nan

    if threshold_type == 'tercile':
        valid_y = y_continuous.dropna()
        if len(valid_y) < 3: # Need at least 3 non-NaN points for meaningful terciles
            # print(f"    Warning: Not enough valid data points ({len(valid_y)}) in y_continuous for tercile. Using median if possible, else NaN.")
            if len(valid_y) > 0:
                threshold = np.median(valid_y)
            else: # No valid data points at all
                return pd.Series([np.nan] * len(y_continuous), index=y_continuous.index), np.nan
        else:
            threshold = np.percentile(valid_y, 100 * 2/3) # Top 1/3 vs bottom 2/3

        binary_labels = (y_continuous >= threshold).astype(float) # Use float to allow NaNs
        binary_labels[y_continuous.isnull()] = np.nan # Ensure original NaNs are preserved
        return binary_labels, threshold
    # elif threshold_type == 'median': # Keep other types if they might be used elsewhere, or remove
    #     threshold = np.median(y_continuous.dropna())
    #     binary_labels = (y_continuous >= threshold).astype(int)
    #     return binary_labels, threshold
    # elif threshold_type == 'quartile':
    #     threshold = np.percentile(y_continuous.dropna(), 75)
    #     binary_labels = (y_continuous >= threshold).astype(int)
    #     return binary_labels, quartile_75
    else:
        raise ValueError(f"Unsupported threshold_type: {threshold_type}. Only 'tercile' is actively supported for this script.")

def plot_performance(y_true, y_pred, signature_name, correlation, mse, set_name, output_dir):
    """Create simplified performance plots: actual vs predicted."""
    file_path_signature_name = sanitize_for_file_path(signature_name)
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    fig.suptitle(f'{set_name} Set Performance: {signature_name}', fontsize=14, fontweight='bold')
    
    # Ensure y_true and y_pred are numpy arrays for min/max operations if they are Series
    y_true_arr = y_true.values if isinstance(y_true, pd.Series) else y_true
    y_pred_arr = y_pred.values if isinstance(y_pred, pd.Series) else y_pred

    ax.scatter(y_true_arr, y_pred_arr, alpha=0.6, color='blue')
    min_val = min(np.nanmin(y_true_arr), np.nanmin(y_pred_arr)) # Use nanmin
    max_val = max(np.nanmax(y_true_arr), np.nanmax(y_pred_arr)) # Use nanmax
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
    ax.set_xlabel('Actual Values (Normalized)')
    ax.set_ylabel('Predicted Values (Normalized)')
    ax.set_title(f'Actual vs Predicted\nCorrelation: {correlation:.3f}, MSE: {mse:.3f}')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    os.makedirs(output_dir, exist_ok=True)
    plot_path = os.path.join(output_dir, f'{file_path_signature_name}_{set_name}_performance.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return plot_path

def plot_roc_analysis(
    y_train_true_continuous, y_train_pred_continuous_for_auc,
    y_val_true_continuous, y_val_pred_continuous_for_auc,
    y_test_true_continuous, y_test_pred_continuous_for_auc,
    signature_name,
    output_dir
):
    """Create ROC analysis plots for Top Tercile threshold for train, validation, and test sets."""
    file_path_signature_name = sanitize_for_file_path(signature_name)
    
    # Adjusted to 3 rows for ROC curves, 1 column. Or 1 row, 3 columns. Let's try 1 row, 3 columns.
    fig, axes = plt.subplots(1, 3, figsize=(21, 6)) # 1 row, 3 columns for ROCs
    fig.suptitle(f'Top Tercile ROC Analysis: {signature_name}', fontsize=16, fontweight='bold')
    
    datasets = [
        ('Train', y_train_true_continuous, y_train_pred_continuous_for_auc, axes[0]),
        ('Validation', y_val_true_continuous, y_val_pred_continuous_for_auc, axes[1]),
        ('Test', y_test_true_continuous, y_test_pred_continuous_for_auc, axes[2])
    ]
    
    auc_scores = {'train': np.nan, 'val': np.nan, 'test': np.nan}

    for i, (set_name, y_true_continuous, y_pred_continuous, roc_ax) in enumerate(datasets):
        auc_score_tercile_set = np.nan
        
        # Ensure y_true_continuous is a pd.Series for consistent NaN handling
        if not isinstance(y_true_continuous, pd.Series):
            y_true_continuous = pd.Series(y_true_continuous)
        if not isinstance(y_pred_continuous, pd.Series):
            # Assuming y_pred_continuous is a numpy array, it should align with y_true_continuous's index if made a Series
            y_pred_continuous = pd.Series(y_pred_continuous, index=y_true_continuous.index)


        # Create binary labels from y_true_continuous for Top Tercile
        y_binary, threshold_value = create_binary_labels(y_true_continuous, threshold_type='tercile')
        
        # Filter out NaNs from y_binary and corresponding y_pred_continuous
        valid_indices = y_binary.notna() & y_pred_continuous.notna()
        y_binary_clean = y_binary[valid_indices]
        y_pred_continuous_clean = y_pred_continuous[valid_indices]

        if len(np.unique(y_binary_clean)) < 2 or len(y_binary_clean) == 0:
            # print(f"    Skipping ROC for {set_name} set due to single class or no data in y_binary_clean for {signature_name}")
            roc_ax.text(0.5, 0.5, f'Single class or no data\nfor {set_name} ROC', ha='center', va='center', transform=roc_ax.transAxes)
        else:
            fpr, tpr, _ = roc_curve(y_binary_clean, y_pred_continuous_clean)
            auc_score_tercile_set = auc(fpr, tpr)
            
            roc_ax.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC (AUC = {auc_score_tercile_set:.3f})') # Simplified label
            roc_ax.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
            roc_ax.set_xlim([0.0, 1.0])
            roc_ax.set_ylim([0.0, 1.05])
            roc_ax.set_xlabel('False Positive Rate', fontsize=12)
            roc_ax.set_ylabel('True Positive Rate', fontsize=12)
            roc_ax.set_title(f'{set_name} Set - Tercile (Thresh: {threshold_value:.3f})', fontsize=12) # fontweight='bold'
            roc_ax.legend(loc="lower right", fontsize=10)
            roc_ax.grid(True, alpha=0.3)

        if set_name == 'Train': auc_scores['train'] = auc_score_tercile_set
        elif set_name == 'Validation': auc_scores['val'] = auc_score_tercile_set
        elif set_name == 'Test': auc_scores['test'] = auc_score_tercile_set
            
    plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust rect for suptitle
    
    os.makedirs(output_dir, exist_ok=True)
    plot_path = os.path.join(output_dir, f'{file_path_signature_name}_roc_train_val_test_tercile.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return plot_path, auc_scores['train'], auc_scores['val'], auc_scores['test']

def plot_boosting_performance(signature_name, eval_results, best_iteration, output_dir):
    """Plots the training and validation performance over boosting rounds."""
    file_path_signature_name = sanitize_for_file_path(signature_name)
    
    # XGBoost eval_results typically has keys like 'validation_0', 'validation_1' if multiple eval_sets
    # If only one eval_set (X_val, y_val) is passed, it's usually 'validation_0'
    # The metric name (e.g., 'rmse') is then a sub-key.
    
    metric_key_to_plot = None
    eval_scores = None

    if eval_results and isinstance(eval_results, dict):
        # Common case: eval_set = [(X_val, y_val)], so key is often 'validation_0'
        # Or if named: eval_set = [ (X_val, y_val, 'validation_set_name')] -> key is 'validation_set_name'
        # For this script, we pass one eval set, so it should be 'validation_0' if not named.
        # Let's assume `model.fit(..., eval_set=[(X_val_scaled, y_val_flat)], ...)`
        # then `eval_results` will be `{'validation_0': {'rmse': [values...]}}`
        
        # Find the first (and likely only) evaluation set key
        eval_set_label = next(iter(eval_results)) if eval_results else None
        if eval_set_label and eval_results[eval_set_label]:
            metric_key_to_plot = next(iter(eval_results[eval_set_label])) # e.g., 'rmse'
            eval_scores = eval_results[eval_set_label][metric_key_to_plot]
        
    if not eval_scores or not metric_key_to_plot:
        print(f"    Warning: Could not extract evaluation scores for boosting plot for {signature_name}. Eval results: {eval_results}")
        return None

    rounds = range(len(eval_scores))

    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(rounds, eval_scores, label=f'Validation {metric_key_to_plot.upper()}') # Using validation set for early stopping
    
    if best_iteration is not None: # best_iteration is 0-indexed
        ax.axvline(x=best_iteration, color='r', linestyle='--', label=f'Best Iteration ({best_iteration+1})') # Display as 1-indexed
    
    ax.set_xlabel("Boosting Round")
    ax.set_ylabel(f"{metric_key_to_plot.upper()} Value")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    fig.suptitle(f"XGBoost Boosting Performance: {signature_name}", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    os.makedirs(output_dir, exist_ok=True)
    plot_path = os.path.join(output_dir, f'{file_path_signature_name}_boosting_performance.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # print(f"    Boosting performance plot saved to: {plot_path}") # Reduced verbosity
    return plot_path

def train_single_signature(
    signature_name,
    X_train_df, y_train_series, # These are for the current signature
    X_val_df, y_val_series,
    X_test_df, y_test_series,
    feature_names, # Full list of feature names from original X_train_full
    output_base_dir
):
    """Train XGBoost model for a single signature using pre-split and pre-normalized data."""
    print(f"  Training for signature: {signature_name}")
    
    file_path_signature_name = sanitize_for_file_path(signature_name)
    model_dir = os.path.join(output_base_dir, "models")
    performance_plot_dir = os.path.join(output_base_dir, "performance_plots")
    roc_plot_dir = os.path.join(output_base_dir, "roc_plots")
    epoch_plot_dir = os.path.join(output_base_dir, "epoch_plots") # For boosting curves
    
    os.makedirs(model_dir, exist_ok=True)
    # Other plot dirs will be created by their respective functions

    try:
        # X data is already cleaned. Y data is already normalized.
        # X Scaling: Fit on X_train_df, transform all X sets for this signature
        x_scaler = StandardScaler()
        X_train_scaled = x_scaler.fit_transform(X_train_df)
        X_val_scaled = x_scaler.transform(X_val_df)
        X_test_scaled = x_scaler.transform(X_test_df)
        
        xgb_eval_metric_name = "rmse"

        model = xgb.XGBRegressor(
            objective='reg:squarederror',
            n_estimators=1000,
            learning_rate=0.01, # Example, consider tuning or using previous values
            max_depth=5,        # Example
            subsample=0.8,
            colsample_bytree=0.8,
            random_state=42,
            tree_method='gpu_hist', # For GPU
            gpu_id=0, # Ensure GPU ID is set
            eval_metric=xgb_eval_metric_name, # Explicitly set for early stopping
            early_stopping_rounds=20 # Increased early stopping rounds
        )
        
        y_train_flat = y_train_series.ravel()
        y_val_flat = y_val_series.ravel()
        y_test_flat = y_test_series.ravel() # For final evaluation

        print("    Training XGBoost model...")
        model.fit(X_train_scaled, y_train_flat,
                  eval_set=[(X_val_scaled, y_val_flat)], # Use validation set for early stopping
                  verbose=False)
        print("    XGBoost training completed.")

        best_iteration = model.best_iteration if hasattr(model, 'best_iteration') and model.best_iteration is not None else model.n_estimators_
        
        # Plot boosting performance (based on validation set)
        boosting_plot_path = None
        if hasattr(model, 'evals_result') and model.evals_result():
            eval_results = model.evals_result()
            boosting_plot_path = plot_boosting_performance(signature_name, eval_results, best_iteration, epoch_plot_dir)
        else:
            print(f"    Could not retrieve evaluation results for plotting boosting performance for {signature_name}.")

        # Predictions
        y_pred_train_continuous = model.predict(X_train_scaled)
        y_pred_val_continuous = model.predict(X_val_scaled)
        y_pred_test_continuous = model.predict(X_test_scaled)

        # Correlation and MSE for each set (using original non-inverted predictions for these)
        train_corr, _ = pearsonr(y_train_flat, y_pred_train_continuous)
        train_mse = mean_squared_error(y_train_flat, y_pred_train_continuous)
        
        val_corr, _ = pearsonr(y_val_flat, y_pred_val_continuous)
        val_mse = mean_squared_error(y_val_flat, y_pred_val_continuous)

        test_corr, _ = pearsonr(y_test_flat, y_pred_test_continuous)
        test_mse = mean_squared_error(y_test_flat, y_pred_test_continuous)

        # Plot performance for test set
        # Ensure y_test_series is passed as pd.Series for consistency if y_pred_test_continuous is an array
        performance_plot_path = plot_performance(y_test_series, pd.Series(y_pred_test_continuous, index=y_test_series.index),
                                                 signature_name, test_corr, test_mse, "Test", performance_plot_dir)
        
        # Inversion logic for AUC based on VALIDATION set correlation
        y_pred_train_for_auc = y_pred_train_continuous
        y_pred_val_for_auc = y_pred_val_continuous
        y_pred_test_for_auc = y_pred_test_continuous
        inverted_predictions_flag = False

        if pd.notna(val_corr) and val_corr < 0:
            print(f"    Negative validation correlation ({val_corr:.4f}). Inverting all predictions for AUC calculation.")
            y_pred_train_for_auc = -y_pred_train_continuous
            y_pred_val_for_auc = -y_pred_val_continuous
            y_pred_test_for_auc = -y_pred_test_continuous
            inverted_predictions_flag = True
        
        # ROC analysis & Tercile AUCs for all three sets
        # Ensure continuous prediction series align with true series indices for roc_plot_analysis
        y_pred_train_series_for_auc = pd.Series(y_pred_train_for_auc, index=y_train_series.index)
        y_pred_val_series_for_auc = pd.Series(y_pred_val_for_auc, index=y_val_series.index)
        y_pred_test_series_for_auc = pd.Series(y_pred_test_for_auc, index=y_test_series.index)

        roc_plot_path, train_auc_tercile, val_auc_tercile, test_auc_tercile = plot_roc_analysis(
            y_train_series, y_pred_train_series_for_auc,
            y_val_series, y_pred_val_series_for_auc,
            y_test_series, y_pred_test_series_for_auc,
            signature_name,
            roc_plot_dir
        )
        
        print(f"    Tercile AUCs - Train: {train_auc_tercile if pd.notna(train_auc_tercile) else 'N/A':.4f}, Val: {val_auc_tercile if pd.notna(val_auc_tercile) else 'N/A':.4f}, Test: {test_auc_tercile if pd.notna(test_auc_tercile) else 'N/A':.4f}")

        # Save model and X-scaler
        model_path_json = os.path.join(model_dir, f"{file_path_signature_name}_model.json") # Simplified name
        x_scaler_path_joblib = os.path.join(model_dir, f"{file_path_signature_name}_x_scaler.joblib")
        
        model.save_model(model_path_json)
        joblib.dump(x_scaler, x_scaler_path_joblib)
        # print(f"    Model saved to {model_path_json}") # Reduced verbosity
        # print(f"    X-Scaler saved to {x_scaler_path_joblib}")

        model_info = {
            'signature_name_original': signature_name,
            'feature_names': feature_names, # These are the columns of X_train_df etc.
            'x_scaler_path': x_scaler_path_joblib,
            'model_path': model_path_json,
            'best_iteration': best_iteration,
            'train_correlation': train_corr, 'train_mse': train_mse, 'train_auc_tercile': train_auc_tercile,
            'val_correlation': val_corr, 'val_mse': val_mse, 'val_auc_tercile': val_auc_tercile,
            'test_correlation': test_corr, 'test_mse': test_mse, 'test_auc_tercile': test_auc_tercile,
            'inverted_predictions_for_auc': inverted_predictions_flag,
            'eval_metric_for_early_stopping': xgb_eval_metric_name,
            'best_score_early_stopping': model.best_score if hasattr(model, 'best_score') else np.nan
        }
        metadata_pt_path = os.path.join(model_dir, f"{file_path_signature_name}_model_metadata.pt")
        torch.save(model_info, metadata_pt_path)
        # print(f"    Metadata saved to {metadata_pt_path}")
        
        return {
            "Signature": signature_name,
            "Train_AUC_Tercile": train_auc_tercile, "Val_AUC_Tercile": val_auc_tercile, "Test_AUC_Tercile": test_auc_tercile,
            "Train_Correlation": train_corr, "Val_Correlation": val_corr, "Test_Correlation": test_corr,
            "Train_MSE": train_mse, "Val_MSE": val_mse, "Test_MSE": test_mse,
            "Inverted_AUC_Preds": inverted_predictions_flag,
            "Best_Iteration": best_iteration,
            "Model_File": model_path_json, "X_Scaler_File": x_scaler_path_joblib, "Metadata_File": metadata_pt_path,
            "Performance_Plot_Test": performance_plot_path, "ROC_Plot": roc_plot_path, "Boosting_Plot": boosting_plot_path
        }
    
    except Exception as e:
        print(f"ERROR processing signature {signature_name}: {e}")
        import traceback
        traceback.print_exc()
        return {
            "Signature": signature_name, "Train_AUC_Tercile": np.nan, "Val_AUC_Tercile": np.nan, "Test_AUC_Tercile": np.nan,
            "Error": str(e)
        }

def main():
    start_time = time.time()
    
    output_base_dir = "results_xgb_retrained_60_20_20_split" # New output directory
    prepared_data_dir = "prepared_data_splits"

    print(f"Cleaning up output directory: {output_base_dir}")
    if os.path.exists(output_base_dir):
        shutil.rmtree(output_base_dir)
    # Create base and subdirectories
    subdirs = ["models", "performance_plots", "roc_plots", "epoch_plots", "summaries"]
    for subdir in subdirs:
        os.makedirs(os.path.join(output_base_dir, subdir), exist_ok=True)
    print(f"Output directory created/cleaned: {os.path.abspath(output_base_dir)}")

    # Load pre-split and pre-normalized data
    print(f"Loading data from {prepared_data_dir}...")
    try:
        X_train_full = pd.read_pickle(os.path.join(prepared_data_dir, "X_train.pkl"))
        X_val_full = pd.read_pickle(os.path.join(prepared_data_dir, "X_val.pkl"))
        X_test_full = pd.read_pickle(os.path.join(prepared_data_dir, "X_test.pkl"))

        Y_train_norm_df = pd.read_pickle(os.path.join(prepared_data_dir, "Y_train_norm.pkl"))
        Y_val_norm_df = pd.read_pickle(os.path.join(prepared_data_dir, "Y_val_norm.pkl"))
        Y_test_norm_df = pd.read_pickle(os.path.join(prepared_data_dir, "Y_test_norm.pkl"))
        print("  Data loaded successfully.")
        print(f"  X_train shape: {X_train_full.shape}, Y_train_norm shape: {Y_train_norm_df.shape}")
        print(f"  X_val shape: {X_val_full.shape}, Y_val_norm shape: {Y_val_norm_df.shape}")
        print(f"  X_test shape: {X_test_full.shape}, Y_test_norm shape: {Y_test_norm_df.shape}")

    except FileNotFoundError as e:
        print(f"ERROR: Could not load data from {prepared_data_dir}. {e}")
        print("Please ensure 'prepare_and_analyze_data_splits.py' was run successfully.")
        return

    all_feature_names = X_train_full.columns.tolist()
    rna_signature_names = Y_train_norm_df.columns.tolist()
    
    all_metrics = []
    
    print(f"\nStarting model retraining for {len(rna_signature_names)} RNA signatures...")
    for rna_signature_name in tqdm(rna_signature_names, desc="Retraining Models"):
        y_train_sig = Y_train_norm_df[rna_signature_name].copy()
        y_val_sig = Y_val_norm_df[rna_signature_name].copy()
        y_test_sig = Y_test_norm_df[rna_signature_name].copy()

        # X data is already aligned by sample IDs with Y data due to the preparation script
        # X_train_full, X_val_full, X_test_full have all features.
        # No further per-signature alignment of X needed here as X_scaler will be fit on X_train_full
        # and used for all signatures (assuming features are consistent).
        # If features differ per model, this needs more complex handling,
        # but current structure implies a global X set.
        # The train_single_signature expects the X data for the current signature.
        # Since X features are common, we pass the full X dataframes for the splits.

        # Basic check for sufficient data for the current signature
        if y_train_sig.isnull().all() or y_val_sig.isnull().all() or y_test_sig.isnull().all():
            print(f"  Skipping {rna_signature_name} due to all NaN Y values in one of the splits.")
            all_metrics.append({"Signature": rna_signature_name, "Error": "All NaN Y values in a split"})
            continue
        if len(y_train_sig.dropna()) < 10 or len(y_val_sig.dropna()) < 5 : # Min samples for train/val
             print(f"  Skipping {rna_signature_name} due to insufficient non-NaN samples: train={len(y_train_sig.dropna())}, val={len(y_val_sig.dropna())}")
             all_metrics.append({"Signature": rna_signature_name, "Error": f"Insufficient non-NaN samples"})
             continue

        metrics = train_single_signature(
            rna_signature_name,
            X_train_full.copy(), y_train_sig, # Pass copies to avoid modification issues
            X_val_full.copy(), y_val_sig,
            X_test_full.copy(), y_test_sig,
            all_feature_names, # Pass the global feature names
            output_base_dir
        )
        if metrics:
            all_metrics.append(metrics)
            
    summary_file_path = os.path.join(output_base_dir, "summaries", "all_signatures_retrained_results.csv")
    if all_metrics:
        metrics_df = pd.DataFrame(all_metrics)
        # Define preferred column order
        cols_order = [
            "Signature", "Train_AUC_Tercile", "Val_AUC_Tercile", "Test_AUC_Tercile",
            "Train_Correlation", "Val_Correlation", "Test_Correlation",
            "Train_MSE", "Val_MSE", "Test_MSE",
            "Inverted_AUC_Preds", "Best_Iteration", "Error",
            "Model_File", "X_Scaler_File", "Metadata_File",
            "Performance_Plot_Test", "ROC_Plot", "Boosting_Plot"
        ]
        # Filter to only include columns that actually exist in the DataFrame, maintaining order
        existing_cols_in_order = [col for col in cols_order if col in metrics_df.columns]
        metrics_df = metrics_df[existing_cols_in_order]
        metrics_df.to_csv(summary_file_path, index=False, float_format='%.4f')
        print(f"\nAll signature metrics saved to {summary_file_path}")
    else:
        print("\nNo metrics were collected. Summary CSV not saved.")
    
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")

if __name__ == '__main__':
    main() 