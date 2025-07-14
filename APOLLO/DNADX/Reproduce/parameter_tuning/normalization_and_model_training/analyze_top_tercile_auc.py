import os
import numpy as np
import pandas as pd
import torch
import re
from sklearn.metrics import roc_auc_score
from model_training_nn import DeepCNANet
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("training_data/rna_signature_score_median_no_norm.pkl")
    segment_score = pd.read_pickle("training_data/cna_segment_score_mean_no_norm.pkl")
    
    print("\nFirst 10 signature names in data:")
    print('\n'.join(signature_score.index[:10]))
    print("\nTotal signatures in data:", len(signature_score.index))
    
    # Convert CNA sample names to match RNA format (replace . with -)
    segment_score.columns = segment_score.columns.str.replace('.', '-')
    
    # Find common samples using exact names
    common_samples = sorted(set(signature_score.columns) & set(segment_score.columns))
    print("\nFound {} common samples".format(len(common_samples)))
    
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

def create_binary_labels(y_continuous):
    """Create binary labels for top tercile vs bottom two terciles."""
    tercile_threshold = np.percentile(y_continuous, 66.67)
    return (y_continuous >= tercile_threshold).astype(int)

def calculate_top_tercile_auc(y_true, y_pred):
    """Calculate AUC for top tercile classification."""
    # Convert continuous predictions to binary labels (top tercile vs rest)
    y_true_binary = create_binary_labels(y_true)
    
    try:
        # Use continuous predictions against binary labels for ROC
        auc = roc_auc_score(y_true_binary, y_pred)
        if auc < 0.5:  # If correlation is negative, invert predictions
            auc = roc_auc_score(y_true_binary, -y_pred)
    except ValueError:
        auc = np.nan
    
    return auc

def analyze_model_performance(model_info, X_train, y_train, X_test, y_test, device='cpu'):
    """Analyze model's top tercile classification performance."""
    # Load model architecture from saved info
    input_size = model_info['architecture']['input_size']
    hidden_sizes = model_info['architecture']['hidden_sizes']
    
    # Initialize model with correct architecture
    model = DeepCNANet(input_size=input_size, hidden_sizes=hidden_sizes).to(device)
    model.load_state_dict(model_info['model'])
    model.eval()
    
    # Apply the same scaling as used during training
    scaler = model_info['scaler']
    X_train_scaled = scaler.transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Get predictions
    with torch.no_grad():
        y_train_pred = model(torch.FloatTensor(X_train_scaled).to(device)).cpu().numpy().flatten()
        y_test_pred = model(torch.FloatTensor(X_test_scaled).to(device)).cpu().numpy().flatten()
    
    # Calculate top tercile AUC
    train_top_tercile_auc = calculate_top_tercile_auc(y_train, y_train_pred)
    test_top_tercile_auc = calculate_top_tercile_auc(y_test, y_test_pred)
    
    return {
        'train_top_tercile_auc': train_top_tercile_auc,
        'test_top_tercile_auc': test_top_tercile_auc
    }

def normalize_name(name):
    """Strip all non-alphanumeric characters and convert to lowercase for matching"""
    # Remove common suffixes first
    name = name.replace('_final_model.pt', '').replace('.pt', '')
    
    # Keep only letters and numbers, convert to lowercase
    normalized = re.sub(r'[^a-zA-Z0-9]', '', name).lower()
    return normalized

def find_best_signature_match(model_file, signature_names):
    """Find the best matching signature for a model file"""
    # Get normalized model name
    model_normalized = normalize_name(model_file)
    
    # Create normalized signature names
    signature_normalized = {normalize_name(sig): sig for sig in signature_names}
    
    # Try exact match first
    if model_normalized in signature_normalized:
        return signature_normalized[model_normalized]
    
    # If no exact match, try partial matches (signature contains model or vice versa)
    for norm_sig, orig_sig in signature_normalized.items():
        if model_normalized in norm_sig or norm_sig in model_normalized:
            # Check if it's a reasonable match (not too different in length)
            if abs(len(model_normalized) - len(norm_sig)) / max(len(model_normalized), len(norm_sig)) < 0.3:
                return orig_sig
    
    return None

def main():
    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # Load data
    print("Loading data...")
    signature_score, segment_score = load_data()
    print("Data loaded successfully!")
    print(f"Signature score shape: {signature_score.shape}")
    print(f"Segment score shape: {segment_score.shape}")
    
    # Print all signature names for analysis
    print(f"\nAll {len(signature_score.index)} RNA signature names:")
    for i, sig_name in enumerate(signature_score.index):
        print(f"{i+1:3d}: {sig_name}")
    
    # Get list of all model files
    models_dir = 'results_nn/models'
    model_files = [f for f in os.listdir(models_dir) if f.endswith('_final_model.pt')]
    total_models = len(model_files)
    print(f"\nAll {total_models} model files:")
    for i, model_file in enumerate(model_files):
        model_normalized = normalize_name(model_file)
        print(f"{i+1:3d}: {model_file} -> {model_normalized}")
    
    # Create a mapping of model names to signature names using robust matching
    print("\nMatching models to signatures...")
    model_to_signature = {}
    unmatched_models = []
    
    for model_file in model_files:
        best_match = find_best_signature_match(model_file, signature_score.index)
        if best_match:
            model_to_signature[model_file] = best_match
            print(f"✓ {model_file} -> {best_match}")
        else:
            unmatched_models.append(model_file)
            print(f"✗ {model_file} -> NO MATCH")
    
    print(f"\nMatching summary:")
    print(f"Successfully matched: {len(model_to_signature)} models")
    print(f"Unmatched: {len(unmatched_models)} models")
    
    if unmatched_models:
        print(f"\nUnmatched models:")
        for model in unmatched_models[:10]:  # Show first 10
            print(f"  {model} -> normalized: {normalize_name(model)}")
        if len(unmatched_models) > 10:
            print(f"  ... and {len(unmatched_models) - 10} more")
    
    # Prepare shared data
    X_data = segment_score.T.values  # Transpose so samples are rows, features are columns
    feature_names = list(segment_score.index)
    
    # Initialize results list
    results = []
    
    print(f"\nProcessing {len(model_to_signature)} matched models...")
    
    for i, (model_file, signature_name) in enumerate(model_to_signature.items(), 1):
        print(f"\nProcessing {i}/{len(model_to_signature)}: {model_file}")
        print(f"  -> Signature: {signature_name}")
        
        try:
            model_path = os.path.join(models_dir, model_file)
            
            # Get signature data
            y_data = signature_score.loc[signature_name].values
            
            # Clean and split data
            X_data_clean = clean_data(X_data)
            X_train, X_test, y_train, y_test = train_test_split(
                X_data_clean, y_data, test_size=0.2, random_state=42, stratify=None
            )
            
            # Load model info
            model_info = torch.load(model_path)
            
            # Analyze model performance
            performance = analyze_model_performance(
                model_info, X_train, y_train, X_test, y_test, device
            )
            
            # Add signature name and original correlation
            performance['signature'] = signature_name
            performance['original_correlation'] = model_info.get('test_correlation', np.nan)
            results.append(performance)
            
            print(f"  Train AUC: {performance['train_top_tercile_auc']:.3f}")
            print(f"  Test AUC: {performance['test_top_tercile_auc']:.3f}")
            print(f"  Original correlation: {performance['original_correlation']:.3f}")
            
        except Exception as e:
            print(f"Error processing {model_file}: {str(e)}")
    
    if results:
        # Convert results to DataFrame
        results_df = pd.DataFrame(results)
        
        # Sort by test AUC
        results_df = results_df.sort_values('test_top_tercile_auc', ascending=False)
        
        # Save results
        output_path = os.path.join('results_nn', 'top_tercile_auc_summary.csv')
        print(f"\nAttempting to save results to: {output_path}")
        try:
            results_df.to_csv(output_path, index=False)
            print(f"Results successfully saved to {output_path}")
        except Exception as e:
            print(f"Error saving results: {str(e)}")
            # Try saving to current directory as fallback
            fallback_path = 'top_tercile_auc_summary.csv'
            print(f"Attempting to save to current directory: {fallback_path}")
            results_df.to_csv(fallback_path, index=False)
        
        # Print summary statistics
        print("\nSummary Statistics for Top Tercile AUC:")
        print("\nTraining Set:")
        print(f"Mean: {results_df['train_top_tercile_auc'].mean():.3f}")
        print(f"Std: {results_df['train_top_tercile_auc'].std():.3f}")
        print("\nTest Set:")
        print(f"Mean: {results_df['test_top_tercile_auc'].mean():.3f}")
        print(f"Std: {results_df['test_top_tercile_auc'].std():.3f}")
        
        # Print top 5 performing signatures
        print("\nTop 5 signatures by test set performance:")
        top_5 = results_df.head()
        for _, row in top_5.iterrows():
            print(f"{row['signature']}: Test AUC = {row['test_top_tercile_auc']:.3f} "
                  f"(Original correlation: {row['original_correlation']:.3f})")
    
    else:
        print("No results were generated. Check for errors above.")

if __name__ == "__main__":
    main() 