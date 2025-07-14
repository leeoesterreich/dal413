#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import xgboost as xgb
import joblib
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from tqdm import tqdm
import torch

def create_binary_labels(y_continuous, threshold_type='median'):
    """Convert continuous values to binary labels based on threshold.
    Copied from model_training_nn_all.py (XGBoost version)."""
    if threshold_type == 'median':
        threshold = np.median(y_continuous)
        binary_labels = (y_continuous >= threshold).astype(int)
        return binary_labels, threshold
    elif threshold_type == 'tercile':
        # Top tercile means values >= 66.67th percentile are considered "high"
        tercile_66 = np.percentile(y_continuous, 100 * (2/3))
        binary_labels = (y_continuous >= tercile_66).astype(int)
        return binary_labels, tercile_66
    elif threshold_type == 'quartile':
        # Top quartile means values >= 75th percentile are considered "high"
        quartile_75 = np.percentile(y_continuous, 75)
        binary_labels = (y_continuous >= quartile_75).astype(int)
        return binary_labels, quartile_75
    else:
        raise ValueError("threshold_type must be 'median', 'tercile', or 'quartile'")

def normalize_name(name):
    # Remove file extension if present
    name = os.path.splitext(name)[0]
    # Remove '_final_model' if present
    name = name.replace('_final_model', '')
    # Replace special characters with underscores
    name = name.replace('.', '_').replace('-', '_').replace(' ', '_')
    return name.lower()  # Convert to lowercase for case-insensitive matching

def load_data():
    print("Loading data...")
    # Load CNA data
    cna_data = pd.read_pickle("validation_data/metabric_scores/cna_signature_score_cleaned.pkl")
    print("CNA data shape:", cna_data.shape)
    print("CNA NaN values before filling:", cna_data.isna().sum().sum())
    cna_data = cna_data.fillna(0)  # Fill NaN with 0
    print("CNA NaN values after filling:", cna_data.isna().sum().sum())
    
    # Load RNA data
    rna_data = pd.read_pickle("validation_data/metabric_scores/rna_signature_score_cleaned.pkl")
    if 'Entrez_Gene_Id' in rna_data.columns:
        rna_data = rna_data.drop('Entrez_Gene_Id', axis=1)
    print("RNA data shape:", rna_data.shape)
    print("RNA NaN values before filling:", rna_data.isna().sum().sum())
    rna_data = rna_data.fillna(0)  # Fill NaN with 0
    print("RNA NaN values after filling:", rna_data.isna().sum().sum())
    
    # Align samples between RNA and CNA data
    common_samples = sorted(list(set(cna_data.columns) & set(rna_data.columns)))
    print(f"\nFound {len(common_samples)} common samples")
    
    cna_data = cna_data[common_samples]
    rna_data = rna_data[common_samples]
    print("After alignment:")
    print("CNA data shape:", cna_data.shape)
    print("RNA data shape:", rna_data.shape)
    
    # Create normalized index mapping for RNA data
    rna_name_map = {normalize_name(idx): idx for idx in rna_data.index}
    print(f"\nFound {len(rna_name_map)} unique normalized RNA signatures")
    
    return cna_data, rna_data, rna_name_map

def validate_models():
    # Load data
    cna_data, rna_data, rna_name_map = load_data()
    
    # Get all model metadata files
    model_dir = "results_xgb/models" 
    # These .pt files are expected to be dictionaries saved by torch.save,
    # containing paths to the actual XGBoost model (.json) and scaler (.joblib), plus feature_names.
    metadata_files = [f for f in os.listdir(model_dir) if f.endswith('_final_model.pt')]
    print(f"\nFound {len(metadata_files)} model metadata files to process from {model_dir}")
    
    results = []
    
    print("\nValidating all models...")
    for metadata_filename in tqdm(metadata_files):
        model_name_for_display = metadata_filename.replace('_final_model.pt', '') # For display and output
        safe_model_name = normalize_name(model_name_for_display) # For matching RNA sig and saving plots
        
        try:
            # Get RNA signature name from model file and normalize it
            if safe_model_name not in rna_name_map:
                print(f"\nSkipping {metadata_filename} - RNA signature '{safe_model_name}' not found in METABRIC RNA data")
                continue
            
            rna_sig_metabric = rna_name_map[safe_model_name]
            
            # Load metadata (which was saved as a .pt file)
            metadata_path = os.path.join(model_dir, metadata_filename)
            if not os.path.exists(metadata_path):
                print(f"\nSkipping {metadata_filename} - metadata file not found at {metadata_path}")
                continue
                
            # Load with map_location to handle potential CPU/GPU differences during saving/loading
            device = 'cpu' # XGBoost runs on CPU for prediction unless explicitly configured for GPU prediction
            model_info = torch.load(metadata_path, map_location=device) 

            required_keys = ['feature_names', 'scaler_path', 'model_path']
            missing_keys = [key for key in required_keys if key not in model_info]
            if missing_keys:
                print(f"\nSkipping {metadata_filename} - metadata is missing keys: {missing_keys}")
                continue
            
            expected_features = model_info['feature_names']
            actual_model_path = model_info['model_path'] # This should be the .json path
            scaler_path = model_info['scaler_path']

            if not os.path.isabs(actual_model_path):
                actual_model_path = os.path.join(model_dir, os.path.basename(actual_model_path))
            if not os.path.isabs(scaler_path):
                scaler_path = os.path.join(model_dir, os.path.basename(scaler_path))

            if not os.path.exists(actual_model_path):
                print(f"\nSkipping {metadata_filename} - XGBoost model file not found: {actual_model_path}")
                continue
            if not os.path.exists(scaler_path):
                print(f"\nSkipping {metadata_filename} - Scaler file not found: {scaler_path}")
                continue
            
            # Feature alignment with validation CNA data
            cna_features_for_model = []
            missing_in_cna_data = []
            cna_data_normalized_index = {normalize_name(idx): idx for idx in cna_data.index}

            for feat_name_from_training in expected_features:
                normalized_feat_name = normalize_name(feat_name_from_training)
                if normalized_feat_name in cna_data_normalized_index:
                    cna_features_for_model.append(cna_data_normalized_index[normalized_feat_name])
                else:
                    missing_in_cna_data.append(feat_name_from_training)
            
            if missing_in_cna_data:
                # If many features are missing, it could be a systematic issue.
                # For now, we will proceed if *some* features are found, but the model might underperform.
                # Crucially, the order and number must match what the scaler and model expect.
                # The current logic requires all expected_features to be found.
                print(f"\nSkipping {metadata_filename} - {len(missing_in_cna_data)} features used in training were not found (after normalization) in validation CNA data.")
                if len(missing_in_cna_data) < 10: print(f"Missing: {missing_in_cna_data[:20]}")
                continue
            
            if len(cna_features_for_model) != len(expected_features):
                print(f"\nSkipping {metadata_filename} - Mismatch in feature counts after alignment. Expected {len(expected_features)}, found {len(cna_features_for_model)}.")
                continue

            # Load Scaler and XGBoost Model
            scaler = joblib.load(scaler_path)
            model = xgb.XGBRegressor() # Initialize new model instance
            model.load_model(actual_model_path) # Load the trained state
            
            print(f"\nProcessing {model_name_for_display}")
            print(f"  Matched METABRIC RNA signature: {rna_sig_metabric}")
            print(f"  Using {len(cna_features_for_model)} CNA features from validation set.")
            
            # Prepare input data: select features in the order they were trained, then transpose for samples as rows
            # The cna_features_for_model list is already in the order of expected_features.
            X_val_metabric = cna_data.loc[cna_features_for_model].values.T 
            X_val_metabric_scaled = scaler.transform(X_val_metabric) # Apply scaling
            
            # Get predictions
            y_pred_continuous = model.predict(X_val_metabric_scaled)
            
            # Get ground truth from METABRIC RNA data
            y_true_continuous = rna_data.loc[rna_sig_metabric].values
            
            # Filter out y_true > 40 as in original script (and corresponding y_pred)
            valid_indices = y_true_continuous <= 40
            y_true_filtered = y_true_continuous[valid_indices]
            y_pred_filtered = y_pred_continuous[valid_indices]

            if len(y_true_filtered) < 10: # Require at least 10 samples after filtering
                print(f"  Skipping {model_name_for_display} - insufficient valid points after filtering y_true<=40 ({len(y_true_filtered)})")
                continue
            
            # Convert to binary classification using TOP TERCILE for y_true_filtered
            y_true_binary, tercile_threshold = create_binary_labels(y_true_filtered, threshold_type='tercile')
            print(f"  Applied top tercile threshold ({tercile_threshold:.3f}) to binarize METABRIC RNA scores.")

            if len(np.unique(y_true_binary)) < 2:
                print(f"  Skipping {model_name_for_display} - Not enough class diversity after top tercile binarization. Only one class present.")
                continue
            
            # Calculate ROC AUC using the continuous predictions (y_pred_filtered)
            fpr, tpr, _ = roc_curve(y_true_binary, y_pred_filtered)
            roc_auc_tercile = auc(fpr, tpr)
                
            results.append((model_name_for_display, rna_sig_metabric, roc_auc_tercile))
            print(f"  Top Tercile AUC = {roc_auc_tercile:.4f} for {model_name_for_display}")
                
            # Save ROC plot
            plt.figure(figsize=(8, 6))
            plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (Top Tercile AUC = {roc_auc_tercile:.3f})')
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'METABRIC Validation ROC for {model_name_for_display}\n(Top Tercile Threshold on True RNA)')
            plt.legend(loc="lower right")
            
            os.makedirs('results_xgb/validation_plots', exist_ok=True)
            # Use safe_model_name which is normalized and suitable for filenames
            plt.savefig(f'results_xgb/validation_plots/{safe_model_name}_roc_tercile_auc.png')
            plt.close()
            
        except Exception as e:
            print(f"\nError processing {metadata_filename}: {str(e)}")
            import traceback
            traceback.print_exc() # Print full traceback for debugging
            continue
    
    if results:
        results_df = pd.DataFrame(results, columns=['ModelName', 'MatchedRNASignature', 'TopTercileAUC'])
        results_df = results_df.sort_values(by='TopTercileAUC', ascending=False)
        output_csv_path = "results_xgb/validation_summary_tercile_auc_xgb.csv"
        results_df.to_csv(output_csv_path, index=False)
        print(f"\nValidation summary saved to {output_csv_path}")
        
        if not results_df.empty:
            print(f"\n--- Summary of METABRIC Validation (Top Tercile AUC) ---")
            print(f"Total models successfully validated: {len(results_df)}")
            print(f"Mean Top Tercile AUC: {results_df['TopTercileAUC'].mean():.4f}")
            print(f"Median Top Tercile AUC: {results_df['TopTercileAUC'].median():.4f}")
            print(f"Std Dev Top Tercile AUC: {results_df['TopTercileAUC'].std():.4f}")
            print(f"Min Top Tercile AUC: {results_df['TopTercileAUC'].min():.4f} (Model: {results_df.iloc[results_df['TopTercileAUC'].idxmin()]['ModelName']})")
            print(f"Max Top Tercile AUC: {results_df['TopTercileAUC'].max():.4f} (Model: {results_df.iloc[results_df['TopTercileAUC'].idxmax()]['ModelName']})")
        
        # Plot AUC distribution
        if len(results_df['TopTercileAUC']) > 1:
            plt.figure(figsize=(10,6))
            plt.hist(results_df['TopTercileAUC'], bins=20, edgecolor='black')
            plt.title('Distribution of Top Tercile AUCs on METABRIC Validation (XGBoost Models)')
            plt.xlabel('Top Tercile AUC')
            plt.ylabel('Number of Models')
            os.makedirs('results_xgb/validation_plots', exist_ok=True)
            plt.savefig('results_xgb/validation_plots/auc_distribution_tercile_xgb.png')
        plt.close()
    else:
        print("\nNo models were successfully validated.")

if __name__ == "__main__":
    validate_models() 