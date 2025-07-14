import pandas as pd
import numpy as np
import xgboost as xgb
import joblib
import os
import torch # Added for loading .pt metadata
from sklearn.metrics import roc_auc_score, roc_curve, auc # Added roc_curve, auc for completeness
from sklearn.preprocessing import StandardScaler # For type hinting Y_scalers
from scipy.stats import pearsonr # For correlation

# Define a file path sanitization function (consistent with training script)
def sanitize_for_file_path(name):
    name = str(name)
    name = name.replace('/', '__SLASH__').replace('\\', '__BACKSLASH__').replace(':', '__COLON__')
    name = "".join(c if c.isalnum() or c in ['_', '-', '.'] else '_' for c in name)
    return name

# Function to create binary labels (copied from model_training_nn_all.py)
def create_binary_labels(y_continuous, threshold_type='tercile'):
    if y_continuous.empty or len(y_continuous.dropna()) == 0:
        print("    Warning: y_continuous is empty or all NaN. Cannot determine threshold.")
        return pd.Series([np.nan] * len(y_continuous), index=y_continuous.index), np.nan
        
    if threshold_type == 'tercile':
        # Ensure there are enough unique, non-NaN values to compute percentile
        valid_y = y_continuous.dropna()
        if len(valid_y) < 3 : # Need at least 3 points for meaningful terciles
             print(f"    Warning: Not enough valid data points ({len(valid_y)}) in y_continuous for tercile. Using median as fallback or returning NaNs.")
             if len(valid_y) > 0:
                 threshold = np.median(valid_y)
             else: # No valid data points at all
                 return pd.Series([np.nan] * len(y_continuous), index=y_continuous.index), np.nan
        else:
            threshold = np.percentile(valid_y, 100 * 2/3) # Top 1/3 vs bottom 2/3

        binary_labels = (y_continuous >= threshold).astype(int)
        # For samples that were NaN in y_continuous, ensure they are NaN in binary_labels
        binary_labels[y_continuous.isnull()] = np.nan 
        return binary_labels, threshold
    else:
        raise ValueError("threshold_type must be 'tercile' for this script.")

def main():
    # Define paths
    base_path = "." 
    prepared_data_dir = os.path.join(base_path, "prepared_data_splits")
    models_base_dir = os.path.join(base_path, "results_xgb/models") # Base for .pt metadata which contains paths
    output_summary_dir = os.path.join(base_path, "results_xgb/summaries")
    os.makedirs(output_summary_dir, exist_ok=True)
    output_csv_path = os.path.join(output_summary_dir, "recalculated_train_val_test_auc_tercile.csv")

    # Load the prepared data splits
    print("Loading prepared data splits...")
    try:
        X_train_full = pd.read_pickle(os.path.join(prepared_data_dir, "X_train.pkl"))
        X_val_full = pd.read_pickle(os.path.join(prepared_data_dir, "X_val.pkl"))
        X_test_full = pd.read_pickle(os.path.join(prepared_data_dir, "X_test.pkl"))

        Y_train_norm_df = pd.read_pickle(os.path.join(prepared_data_dir, "Y_train_norm.pkl"))
        Y_val_norm_df = pd.read_pickle(os.path.join(prepared_data_dir, "Y_val_norm.pkl"))
        Y_test_norm_df = pd.read_pickle(os.path.join(prepared_data_dir, "Y_test_norm.pkl"))
        
        # Y_scalers = joblib.load(os.path.join(prepared_data_dir, "Y_scalers.joblib")) # Not directly used here, Y is already scaled
        print(f"  Loaded X_train_full: {X_train_full.shape}, Y_train_norm_df: {Y_train_norm_df.shape}")
        print(f"  Loaded X_val_full: {X_val_full.shape}, Y_val_norm_df: {Y_val_norm_df.shape}")
        print(f"  Loaded X_test_full: {X_test_full.shape}, Y_test_norm_df: {Y_test_norm_df.shape}")

    except FileNotFoundError as e:
        print(f"Error loading prepared data splits: {e}. Please run prepare_and_analyze_data_splits.py first.")
        print(f"Expected data in: {prepared_data_dir}")
        return

    RNA_SIGNATURE_NAMES = Y_train_norm_df.columns.tolist() # Get signatures from the loaded Y data
    results_list = []

    print(f"\nFound {len(RNA_SIGNATURE_NAMES)} signatures to process.")

    for rna_signature_name in RNA_SIGNATURE_NAMES:
        print(f"\nProcessing signature: {rna_signature_name}")
        
        # Sanitize signature name to find the .pt metadata file
        file_path_signature_name = sanitize_for_file_path(rna_signature_name)
        metadata_pt_path = os.path.join(models_base_dir, f"{file_path_signature_name}_final_model.pt")

        if not os.path.exists(metadata_pt_path):
            print(f"  Metadata file not found: {metadata_pt_path}. Skipping.")
            results_list.append({
                "Signature": rna_signature_name, "Train_AUC_Tercile": np.nan, "Val_AUC_Tercile": np.nan,
                "Test_AUC_Tercile": np.nan, "Test_Correlation": np.nan, 
                "Inverted_Predictions": False, "Error": "Metadata .pt file not found"
            })
            continue

        try:
            # Load metadata, which contains paths to model and scaler, and feature names
            model_info = torch.load(metadata_pt_path)
            model_path_json = model_info['model_path']
            scaler_path_joblib = model_info['scaler_path']
            feature_names_for_model = model_info['feature_names']
            
            # Ensure paths from .pt are absolute or relative to a known root if necessary
            # Assuming paths in .pt are directly usable or relative to project root
            if not os.path.isabs(model_path_json): model_path_json = os.path.join(base_path, model_path_json)
            if not os.path.isabs(scaler_path_joblib): scaler_path_joblib = os.path.join(base_path, scaler_path_joblib)


            if not (os.path.exists(model_path_json) and os.path.exists(scaler_path_joblib)):
                print(f"  Model or X-scaler not found based on paths in {metadata_pt_path}. Skipping.")
                print(f"    Expected model: {model_path_json}")
                print(f"    Expected X-scaler: {scaler_path_joblib}")
                results_list.append({
                    "Signature": rna_signature_name, "Train_AUC_Tercile": np.nan, "Val_AUC_Tercile": np.nan,
                    "Test_AUC_Tercile": np.nan, "Test_Correlation": np.nan,
                    "Inverted_Predictions": False, "Error": "Model or X-scaler file not found"
                })
                continue

            # Load X-scaler and XGBoost model
            x_scaler = joblib.load(scaler_path_joblib)
            model = xgb.XGBRegressor()
            model.load_model(model_path_json)
            print(f"  Loaded model: {model_path_json}")
            print(f"  Loaded X-scaler: {scaler_path_joblib}")
            print(f"  Model expects {len(feature_names_for_model)} features.")

            # Prepare X data for this model (feature selection and scaling)
            X_train_model = X_train_full[feature_names_for_model]
            X_val_model = X_val_full[feature_names_for_model]
            X_test_model = X_test_full[feature_names_for_model]

            X_train_scaled = x_scaler.transform(X_train_model)
            X_val_scaled = x_scaler.transform(X_val_model)
            X_test_scaled = x_scaler.transform(X_test_model)
            print(f"  X data prepared and scaled. Train shape: {X_train_scaled.shape}, Val shape: {X_val_scaled.shape}, Test shape: {X_test_scaled.shape}")

            # Get corresponding Y data (already normalized)
            y_train_norm_sig = Y_train_norm_df[rna_signature_name].copy()
            y_val_norm_sig = Y_val_norm_df[rna_signature_name].copy()
            y_test_norm_sig = Y_test_norm_df[rna_signature_name].copy()

            # Make continuous predictions
            y_pred_train_continuous = model.predict(X_train_scaled)
            y_pred_val_continuous = model.predict(X_val_scaled)
            y_pred_test_continuous = model.predict(X_test_scaled)

            # Determine if inversion is needed based on TEST set correlation
            # Ensure y_test_norm_sig and y_pred_test_continuous are aligned and have no NaNs for correlation
            valid_indices_test = y_test_norm_sig.notna() & ~np.isnan(y_pred_test_continuous)
            if valid_indices_test.sum() < 2: # Need at least 2 points for correlation
                test_corr = np.nan
                print(f"  Warning: Not enough valid data points to calculate test correlation for {rna_signature_name}.")
            else:
                test_corr, _ = pearsonr(y_test_norm_sig[valid_indices_test], y_pred_test_continuous[valid_indices_test])
            
            y_pred_train_for_auc = y_pred_train_continuous
            y_pred_val_for_auc = y_pred_val_continuous
            y_pred_test_for_auc = y_pred_test_continuous
            inverted_flag = False

            if pd.notna(test_corr) and test_corr < 0:
                # Invert predictions if correlation is negative (assuming higher score = higher probability of being in top tercile)
                # The original XGBoost model predicts continuous values. If these are negatively correlated with the true scores,
                # then a higher predicted value means a *lower* true score. For AUC where positive class is "high score",
                # we need to ensure that higher *prediction_for_auc* corresponds to higher likelihood of positive class.
                # If original predictions are anti-correlated, 1 - pred (if preds are ~0-1) or -pred (if preds are centered around 0) might be used.
                # Given XGBoost regressor output is unbounded, -pred is safer.
                y_pred_train_for_auc = -y_pred_train_continuous 
                y_pred_val_for_auc = -y_pred_val_continuous
                y_pred_test_for_auc = -y_pred_test_continuous
                inverted_flag = True
                print(f"  Negative test correlation ({test_corr:.3f}). Inverting predictions for AUC calculation.")

            # Create binary labels based on terciles of *normalized actual* scores for each set
            y_train_binary, train_thresh = create_binary_labels(y_train_norm_sig, threshold_type='tercile')
            y_val_binary, val_thresh = create_binary_labels(y_val_norm_sig, threshold_type='tercile')
            y_test_binary, test_thresh = create_binary_labels(y_test_norm_sig, threshold_type='tercile')
            
            print(f"  Thresholds for tercile - Train: {train_thresh:.3f}, Val: {val_thresh:.3f}, Test: {test_thresh:.3f}")

            # Calculate AUCs
            current_results = {"Signature": rna_signature_name, "Test_Correlation": test_corr, "Inverted_Predictions": inverted_flag}
            for set_name, y_binary, y_pred_for_auc, y_true_continuous in [
                ('Train', y_train_binary, y_pred_train_for_auc, y_train_norm_sig),
                ('Val', y_val_binary, y_pred_val_for_auc, y_val_norm_sig),
                ('Test', y_test_binary, y_pred_test_for_auc, y_test_norm_sig)
            ]:
                auc_val = np.nan
                # Ensure y_binary and y_pred_for_auc are aligned and have no NaNs for AUC
                valid_auc_indices = y_binary.notna() & ~np.isnan(y_pred_for_auc)
                y_binary_clean = y_binary[valid_auc_indices]
                y_pred_for_auc_clean = y_pred_for_auc[valid_auc_indices]

                if len(np.unique(y_binary_clean.dropna())) > 1 and len(y_binary_clean.dropna()) > 0:
                    try:
                        auc_val = roc_auc_score(y_binary_clean, y_pred_for_auc_clean)
                    except ValueError as e_auc:
                        print(f"    Warning: ROC AUC calculation error for {set_name} ({rna_signature_name}): {e_auc}")
                else:
                    print(f"    Warning: {set_name} labels for {rna_signature_name} (binary, cleaned) are monolithic or empty. Cannot calculate AUC.")
                current_results[f"{set_name}_AUC_Tercile"] = auc_val
                print(f"    {set_name} AUC Tercile: {auc_val if pd.notna(auc_val) else 'N/A':<6}")
            
            results_list.append(current_results)

        except Exception as e:
            print(f"  Error processing {rna_signature_name}: {e}")
            import traceback
            traceback.print_exc()
            results_list.append({
                "Signature": rna_signature_name, "Train_AUC_Tercile": np.nan, "Val_AUC_Tercile": np.nan,
                "Test_AUC_Tercile": np.nan, "Test_Correlation": np.nan, 
                "Inverted_Predictions": False, "Error": str(e)
            })

    # Save results to CSV
    results_df = pd.DataFrame(results_list)
    # Reorder columns for clarity
    cols_order = ["Signature", "Train_AUC_Tercile", "Val_AUC_Tercile", "Test_AUC_Tercile", 
                  "Test_Correlation", "Inverted_Predictions", "Error"]
    results_df = results_df[[col for col in cols_order if col in results_df.columns]] # Ensure only existing columns

    results_df.to_csv(output_csv_path, index=False, float_format='%.4f')
    print(f"\nResults saved to {output_csv_path}")

if __name__ == "__main__":
    main() 