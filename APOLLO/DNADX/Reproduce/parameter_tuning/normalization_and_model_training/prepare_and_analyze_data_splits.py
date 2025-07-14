import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt
import os
import joblib

# --- Data Loading and Initial Processing (adapted from model_training_nn_all.py) ---
def load_and_align_data(cna_path, rna_path):
    print("Loading data...")
    signature_score_df = pd.read_pickle(rna_path)  # Signatures as rows, samples as columns
    segment_score_df = pd.read_pickle(cna_path)  # Segments as rows, samples as columns

    print("\\nOriginal Sample name formats:")
    print("RNA signature score sample names (first 5):", list(signature_score_df.columns[:5]))
    print("CNA segment score sample names (first 5):", list(segment_score_df.columns[:5]))

    # Convert CNA sample names to match RNA format (replace . with -)
    # segment_score_df.columns = segment_score_df.columns.str.replace('.', '-', regex=False) # Removed as per user feedback

    common_samples = sorted(list(set(signature_score_df.columns) & set(segment_score_df.columns)))
    print(f"\\nFound {len(common_samples)} common samples")
    if not common_samples:
        raise ValueError("No common samples found between CNA and RNA data. Check sample naming.")
    print("First 5 common samples:", common_samples[:5])

    # Align by common samples
    signature_score_aligned = signature_score_df[common_samples]
    segment_score_aligned = segment_score_df[common_samples]

    # Transpose so samples are rows
    # X: samples x segments
    # Y: samples x signatures
    X_full = segment_score_aligned.T
    Y_full_df = signature_score_aligned.T

    print(f"X_full shape (samples x segments): {X_full.shape}")
    print(f"Y_full_df shape (samples x signatures): {Y_full_df.shape}")
    
    return X_full, Y_full_df, common_samples

def clean_cna_data(X_df):
    print("\\nCleaning CNA data (X)...")
    X_values = X_df.values.astype(float)
    X_values = np.where(np.isinf(X_values), np.nan, X_values)
    
    # Impute NaN values with column means
    imputer = SimpleImputer(strategy='mean')
    X_cleaned_values = imputer.fit_transform(X_values)
    
    X_cleaned_df = pd.DataFrame(X_cleaned_values, index=X_df.index, columns=X_df.columns)
    print(f"NaNs in X after cleaning: {X_cleaned_df.isnull().sum().sum()}")
    return X_cleaned_df

# --- Distribution Analysis (adapted from analyze_dataset_distributions.py) ---
def analyze_y_split_distributions(y_dict, output_dir, base_filename="Y_normalized_distributions"):
    os.makedirs(output_dir, exist_ok=True)
    num_signatures_to_plot = min(5, y_dict['train'].shape[1] if isinstance(y_dict['train'], pd.DataFrame) else 1)
    
    if isinstance(y_dict['train'], pd.DataFrame):
        signatures_to_plot = y_dict['train'].columns[:num_signatures_to_plot]
    else: # Assuming single series if not a DataFrame (though current plan is DataFrame)
        signatures_to_plot = [y_dict['train'].name if hasattr(y_dict['train'], 'name') else "Signature_0"]
        # Temporarily convert to DataFrame for consistent plotting logic if single series
        if not isinstance(y_dict['train'], pd.DataFrame):
             for split_name in y_dict:
                y_dict[split_name] = pd.DataFrame(y_dict[split_name])


    for sig_idx, signature_name in enumerate(signatures_to_plot):
        plt.figure(figsize=(15, 5))
        plot_idx = 1
        for split_name, y_data_df in y_dict.items():
            plt.subplot(1, len(y_dict), plot_idx)
            current_y_series = y_data_df.iloc[:, sig_idx] if isinstance(y_data_df, pd.DataFrame) else y_data_df # Use sig_idx for multi-column Y
            
            finite_values = current_y_series.dropna()
            finite_values = finite_values[np.isfinite(finite_values)]

            if finite_values.empty:
                print(f"No finite values for {signature_name} in {split_name} set for plotting.")
                plt.text(0.5, 0.5, "No finite data", ha='center', va='center')
            else:
                plt.hist(finite_values, bins=30, edgecolor='black', alpha=0.7, density=True)
                mean_val = finite_values.mean()
                median_val = finite_values.median()
                std_val = finite_values.std()
                plt.axvline(mean_val, color='red', linestyle='dashed', linewidth=1, label=f'Mean: {mean_val:.2f}')
                plt.axvline(median_val, color='green', linestyle='dashed', linewidth=1, label=f'Median: {median_val:.2f}')
                print(f"  Stats for {split_name} - {signature_name}: Mean={mean_val:.2f}, Median={median_val:.2f}, Std={std_val:.2f}")

            plt.title(f'{split_name} Set - {signature_name.replace("_", " ")}')
            plt.xlabel('Normalized Y Value')
            plt.ylabel('Density')
            plt.legend()
            plt.grid(True, alpha=0.5)
            plot_idx += 1
        
        plt.suptitle(f'Distribution of Normalized Y for Signature: {signature_name.replace("_", " ")}', fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plot_filename = os.path.join(output_dir, f"{base_filename}_signature_{sig_idx}_{sanitize_for_file_path(signature_name)}.png")
        plt.savefig(plot_filename)
        print(f"Distribution plot saved to: {plot_filename}")
        plt.close()

def sanitize_for_file_path(name):
    name = str(name) # Ensure it's a string
    name = name.replace('/', '__SLASH__').replace('\\\\', '__BACKSLASH__').replace(':', '__COLON__')
    name = "".join(c if c.isalnum() or c in ['_', '-', '.'] else '_' for c in name)
    return name

# --- Main Data Preparation Function ---
def prepare_data(cna_path, rna_path, output_data_dir="prepared_data_splits", output_plot_dir="prepared_data_plots", random_seed=42):
    os.makedirs(output_data_dir, exist_ok=True)
    os.makedirs(output_plot_dir, exist_ok=True)

    # 1. Load and align data
    X_full_raw, Y_full_raw_df, common_samples = load_and_align_data(cna_path, rna_path)

    # 2. Clean CNA data (X)
    X_full_cleaned = clean_cna_data(X_full_raw)

    # 3. Split sample IDs (60% train, 20% validation, 20% test)
    # Ensure Y_full_raw_df is indexed by common_samples if not already (it should be from load_and_align)
    if not X_full_cleaned.index.equals(Y_full_raw_df.index):
        print("Warning: X_full_cleaned and Y_full_raw_df indices do not match after alignment and cleaning. Re-aligning.")
        common_idx = X_full_cleaned.index.intersection(Y_full_raw_df.index)
        X_full_cleaned = X_full_cleaned.loc[common_idx]
        Y_full_raw_df = Y_full_raw_df.loc[common_idx]
        print(f"Re-aligned to {len(common_idx)} samples.")


    sample_ids = X_full_cleaned.index
    train_ids, temp_ids = train_test_split(sample_ids, test_size=0.4, random_state=random_seed) # 60% train, 40% temp
    val_ids, test_ids = train_test_split(temp_ids, test_size=0.5, random_state=random_seed)    # 20% val, 20% test (50% of 40%)

    print(f"\\nData Splitting:")
    print(f"Total samples: {len(sample_ids)}")
    print(f"Training samples: {len(train_ids)}")
    print(f"Validation samples: {len(val_ids)}")
    print(f"Test samples: {len(test_ids)}")

    # 4. Create data splits based on IDs
    X_train = X_full_cleaned.loc[train_ids]
    X_val = X_full_cleaned.loc[val_ids]
    X_test = X_full_cleaned.loc[test_ids]

    Y_train_raw = Y_full_raw_df.loc[train_ids]
    Y_val_raw = Y_full_raw_df.loc[val_ids]
    Y_test_raw = Y_full_raw_df.loc[test_ids]

    # 5. Normalize Y (RNA scores) - Fit on train, transform all
    # Each signature (column in Y) is scaled independently
    Y_scalers = {}
    Y_train_norm = pd.DataFrame(index=Y_train_raw.index, columns=Y_train_raw.columns)
    Y_val_norm = pd.DataFrame(index=Y_val_raw.index, columns=Y_val_raw.columns)
    Y_test_norm = pd.DataFrame(index=Y_test_raw.index, columns=Y_test_raw.columns)

    print("\\nNormalizing Y data (per signature, fit on training set)...")
    for signature in Y_full_raw_df.columns:
        scaler = StandardScaler()
        Y_train_norm[signature] = scaler.fit_transform(Y_train_raw[[signature]])
        Y_val_norm[signature] = scaler.transform(Y_val_raw[[signature]])
        Y_test_norm[signature] = scaler.transform(Y_test_raw[[signature]])
        Y_scalers[signature] = scaler
    
    print(f"NaNs in Y_train_norm: {Y_train_norm.isnull().sum().sum()}")
    print(f"NaNs in Y_val_norm: {Y_val_norm.isnull().sum().sum()}")
    print(f"NaNs in Y_test_norm: {Y_test_norm.isnull().sum().sum()}")


    # 6. Save data splits and Y scalers
    print("\\nSaving data splits and Y scalers...")
    pd.to_pickle(X_train, os.path.join(output_data_dir, "X_train.pkl"))
    pd.to_pickle(X_val, os.path.join(output_data_dir, "X_val.pkl"))
    pd.to_pickle(X_test, os.path.join(output_data_dir, "X_test.pkl"))

    pd.to_pickle(Y_train_norm, os.path.join(output_data_dir, "Y_train_norm.pkl"))
    pd.to_pickle(Y_val_norm, os.path.join(output_data_dir, "Y_val_norm.pkl"))
    pd.to_pickle(Y_test_norm, os.path.join(output_data_dir, "Y_test_norm.pkl"))
    
    joblib.dump(Y_scalers, os.path.join(output_data_dir, "Y_scalers.joblib"))
    joblib.dump(train_ids, os.path.join(output_data_dir, "train_ids.joblib"))
    joblib.dump(val_ids, os.path.join(output_data_dir, "val_ids.joblib"))
    joblib.dump(test_ids, os.path.join(output_data_dir, "test_ids.joblib"))

    print(f"Data saved to: {os.path.abspath(output_data_dir)}")

    # 7. Analyze and plot Y distributions
    print("\\nAnalyzing Y distributions for splits...")
    y_distributions_for_plot = {
        'train': Y_train_norm,
        'validation': Y_val_norm,
        'test': Y_test_norm
    }
    analyze_y_split_distributions(y_distributions_for_plot, output_plot_dir)

    print("\\nData preparation and analysis complete.")

if __name__ == "__main__":
    CNA_DATA_PATH = "training_data/cna_segment_score_mean_no_norm_renamed.pkl"
    RNA_DATA_PATH = "training_data/rna_signature_score_median_no_norm.pkl"
    
    # Check if input files exist
    if not os.path.exists(CNA_DATA_PATH):
        print(f"ERROR: CNA data file not found at {CNA_DATA_PATH}")
        exit(1)
    if not os.path.exists(RNA_DATA_PATH):
        print(f"ERROR: RNA data file not found at {RNA_DATA_PATH}")
        exit(1)

    prepare_data(CNA_DATA_PATH, RNA_DATA_PATH) 