import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from scipy.stats import pearsonr
from sklearn.metrics import roc_curve, auc, mean_squared_error
import matplotlib.pyplot as plt
import os
import pickle

def clean_data(X):
    """Clean data by handling NaN and infinite values"""
    # Convert to numpy array if not already and ensure float type
    X = np.array(X, dtype=float)
    
    # Handle 1D arrays
    if len(X.shape) == 1:
        X = X.reshape(-1, 1)
    
    # Replace inf and NaN with 0
    X = np.where(np.isinf(X) | np.isnan(X), 0, X)
    
    # Return to original shape if 1D
    if len(X.shape) == 2 and X.shape[1] == 1:
        X = X.ravel()
    
    return X

def print_data_ranges(X, y, name=""):
    """Print the range of values in X and y"""
    print(f"\n{name} data ranges:")
    print(f"X range: [{np.min(X):.2f}, {np.max(X):.2f}]")
    print(f"y range: [{np.min(y):.2f}, {np.max(y):.2f}]")
    print(f"X > 40: {np.sum(np.abs(X) > 40)} values")
    print(f"y > 40: {np.sum(np.abs(y) > 40)} values")

def remove_outliers(X, y, threshold=100):
    """Remove samples where either X or y values exceed the threshold"""
    print_data_ranges(X, y, "Before outlier removal")
    
    # Find samples to keep
    valid_samples = np.all(np.abs(X) <= threshold, axis=1) & (np.abs(y) <= threshold)
    
    # Print number of removed samples
    removed_count = len(y) - np.sum(valid_samples)
    if removed_count > 0:
        print(f"Removed {removed_count} samples with values > {threshold}")
    
    X_clean = X[valid_samples]
    y_clean = y[valid_samples]
    
    print_data_ranges(X_clean, y_clean, "After outlier removal")
    return X_clean, y_clean

def load_data():
    """Load and prepare the training and validation data"""
    print("Loading training data...")
    rna_data = pd.read_pickle("training_data/rna_signature_score_median_no_norm.pkl")
    cna_data = pd.read_pickle("training_data/cna_segment_score_mean_no_norm.pkl")
    
    # Standardize TCGA sample names
    cna_data.columns = cna_data.columns.str.replace('.', '-')
    
    print("\nLoading validation data...")
    val_rna_data = pd.read_pickle("validation_data/metabric_scores/rna_signature_score_cleaned.pkl")
    val_cna_data = pd.read_pickle("validation_data/metabric_scores/cna_signature_score_cleaned.pkl")
    
    # Find common features between training and validation CNA data
    common_features = sorted(set(cna_data.index) & set(val_cna_data.index))
    print(f"Found {len(common_features)} common CNA features")
    
    # Keep only common features
    cna_data = cna_data.loc[common_features]
    val_cna_data = val_cna_data.loc[common_features]
    
    # Find common samples within each dataset
    train_common_samples = sorted(set(rna_data.columns) & set(cna_data.columns))
    val_common_samples = sorted(set(val_rna_data.columns) & set(val_cna_data.columns))
    print(f"Found {len(train_common_samples)} common training samples")
    print(f"Found {len(val_common_samples)} common validation samples")
    
    # Align data
    rna_aligned = rna_data[train_common_samples]
    cna_aligned = cna_data[train_common_samples]
    val_rna_aligned = val_rna_data[val_common_samples]
    val_cna_aligned = val_cna_data[val_common_samples]
    
    return rna_aligned, cna_aligned, val_rna_aligned, val_cna_aligned

def train_and_evaluate_model(X, y, X_val, y_val, model_name):
    """Train an ElasticNet model and evaluate on test and validation sets"""
    print(f"\nTraining model: {model_name}")
    print(f"Input data shapes: X={X.shape}, y={y.shape}")
    print(f"Validation data shapes: X_val={X_val.shape}, y_val={y_val.shape}")
    
    # Clean all datasets (replace NaN with 0)
    X = clean_data(X)
    y = clean_data(y)
    X_val = clean_data(X_val)
    y_val = clean_data(y_val)
    
    # Remove outliers from training data
    X, y = remove_outliers(X, y, threshold=40)
    
    # Remove outliers from validation data
    X_val, y_val = remove_outliers(X_val, y_val, threshold=40)
    
    # Split training data
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.3,
        random_state=42
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    X_val_scaled = scaler.transform(X_val)
    
    # Train model
    model = ElasticNetCV(
        l1_ratio=0.5,
        alphas=np.logspace(-3, 1, 100),
        cv=5,
        random_state=42,
        max_iter=10000,
        selection='random'
    )
    
    # Fit model
    model.fit(X_train_scaled, y_train)
    print(f"Best alpha: {model.alpha_}")
    
    # Make predictions
    y_test_pred = model.predict(X_test_scaled)
    y_val_pred = model.predict(X_val_scaled)
    
    # Calculate correlations
    test_corr, _ = pearsonr(y_test, y_test_pred)
    val_corr, _ = pearsonr(y_val, y_val_pred)
    
    print(f"Test correlation: {test_corr:.3f}")
    print(f"Validation correlation: {val_corr:.3f}")
    
    # Calculate ROC curves using upper 1/3 threshold
    test_thresh = np.percentile(y_test, [33.33, 66.67])
    y_test_binary = np.where(y_test >= test_thresh[1], 1, 0)
    fpr_test, tpr_test, _ = roc_curve(y_test_binary, y_test_pred)
    roc_auc_test = auc(fpr_test, tpr_test)
    
    val_thresh = np.percentile(y_val, [33.33, 66.67])
    y_val_binary = np.where(y_val >= val_thresh[1], 1, 0)
    fpr_val, tpr_val, _ = roc_curve(y_val_binary, y_val_pred)
    roc_auc_val = auc(fpr_val, tpr_val)
    
    print(f"Test ROC AUC: {roc_auc_test:.3f}")
    print(f"Validation ROC AUC: {roc_auc_val:.3f}")
    
    # Create ROC plot
    plt.figure(figsize=(8, 8))
    
    # Plot ROC curves
    plt.plot(fpr_test, tpr_test, 'r-', label=f'Test (AUC = {roc_auc_test:.3f})')
    plt.plot(fpr_val, tpr_val, 'b-', label=f'Validation (AUC = {roc_auc_val:.3f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    
    # Set title based on model name
    if model_name == 'UNC RB LOH':
        title = 'RB-LOH'
    elif model_name == 'Basal Signaling':
        title = 'Basal signaling'
    elif model_name == 'Estrogen Signaling':
        title = 'Estrogen signaling'
    else:
        title = model_name
    
    plt.title(title, fontsize=14, pad=10)
    plt.legend(loc="lower right", fontsize=14, frameon=True, framealpha=1)
    plt.grid(True, alpha=0.3)
    
    # Create ROC_plots directory if it doesn't exist
    os.makedirs('ROC_plots', exist_ok=True)
    plt.savefig(f'ROC_plots/{model_name.lower().replace(" ", "_")}_roc.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return model, roc_auc_test, roc_auc_val

def save_model(model, model_name):
    """Save model in Python 3 compatible format"""
    os.makedirs('ROC_plots/models', exist_ok=True)
    model_path = f'ROC_plots/models/{model_name.lower().replace(" ", "_")}_model.pkl'
    
    with open(model_path, 'wb') as f:
        pickle.dump({'model': model}, f, protocol=3)  # Use protocol 3 for Python 3 compatibility
    print(f"Model saved to {model_path}")

def main():
    # Load data
    rna_data, cna_data, val_rna_data, val_cna_data = load_data()
    
    # Models to train
    models_to_train = {
        'Basal Signaling': 'SMID_BREAST_CANCER_BASAL_UP',
        'Estrogen Signaling': 'SMID_BREAST_CANCER_BASAL_DN',
        'UNC RB LOH': 'UNC_RB_LOH',
        'Oncotype DX': 'GHI_RS_Model_NJEM.2004_PMID.15591335'
    }
    
    # Train and save models
    results = {}
    for model_name, signature in models_to_train.items():
        # Find the RNA signature
        matching_signatures = [col for col in rna_data.index if signature in col]
        if not matching_signatures:
            print(f"Warning: No matching signature found for {signature}")
            continue
            
        signature_name = matching_signatures[0]
        print(f"\nProcessing {model_name} using signature: {signature_name}")
        
        # Get the RNA signature values
        y = rna_data.loc[signature_name].values
        X = cna_data.values.T
        
        # Get validation data
        matching_val_signatures = [col for col in val_rna_data.index if signature in col]
        if not matching_val_signatures:
            print(f"Warning: No matching validation signature found for {signature}")
            continue
            
        val_signature_name = matching_val_signatures[0]
        y_val = val_rna_data.loc[val_signature_name].values
        X_val = val_cna_data.values.T
        
        # Train and evaluate model
        model, test_auc, val_auc = train_and_evaluate_model(X, y, X_val, y_val, model_name)
        
        # Save model
        save_model(model, model_name)
        
        results[model_name] = (test_auc, val_auc)
    
    # Print summary
    print("\nResults Summary:")
    for model_name, (test_auc, val_auc) in results.items():
        print(f"{model_name}:")
        print(f"  Test AUC = {test_auc:.3f}")
        print(f"  Validation AUC = {val_auc:.3f}")

if __name__ == '__main__':
    main() 