import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import os
import sys

# Print Python version for debugging
print("Python version: {}".format(sys.version))

def load_model(model_path):
    """Load a saved scikit-learn model"""
    try:
        # Try loading with protocol 2 for Python 2 compatibility
        with open(model_path, 'rb') as f:
            if sys.version_info[0] >= 3:
                # For Python 3, try different encodings
                try:
                    model_dict = pickle.load(f, encoding='bytes')
                    if isinstance(model_dict, dict) and b'model' in model_dict:
                        return model_dict[b'model']
                    return model_dict
                except:
                    f.seek(0)
                    model_dict = pickle.load(f, encoding='latin1')
                    if isinstance(model_dict, dict) and 'model' in model_dict:
                        return model_dict['model']
                    return model_dict
            else:
                # For Python 2, load directly
                model_dict = pickle.load(f)
                if isinstance(model_dict, dict) and 'model' in model_dict:
                    return model_dict['model']
                return model_dict
    except Exception as e:
        print("Error loading model: {}".format(str(e)))
        raise

def load_data(file_path):
    """Load data with appropriate protocol"""
    try:
        return pd.read_pickle(file_path)
    except Exception as e:
        print("Error loading data: {}".format(str(e)))
        # If the error is due to protocol version, try loading with protocol 2
        try:
            with open(file_path, 'rb') as f:
                return pickle.load(f, encoding='bytes')
        except Exception as e2:
            print("Error with alternative loading method: {}".format(str(e2)))
            raise

def standardize_tcga_id(sample_id):
    """Standardize TCGA sample IDs by converting between dot and hyphen format"""
    if isinstance(sample_id, str):
        if 'TCGA-' in sample_id:
            return sample_id
        elif 'TCGA.' in sample_id:
            return sample_id.replace('.', '-')
    return sample_id

def check_data_structure(data, name):
    """Check and print data structure information"""
    print("\n{} data structure:".format(name))
    print("Shape: {}".format(data.shape))
    print("Index type: {}".format(type(data.index)))
    print("First few index entries: {}".format(list(data.index[:5])))
    if isinstance(data, pd.DataFrame):
        print("Columns type: {}".format(type(data.columns)))
        print("First few column names: {}".format(list(data.columns[:5])))

def prepare_data(rna_data, cna_data):
    """Prepare and align RNA and CNA data"""
    # Check initial data structure
    check_data_structure(rna_data, "RNA")
    check_data_structure(cna_data, "CNA")
    
    # Standardize sample names
    if isinstance(rna_data.columns, pd.Index):
        rna_data.columns = [standardize_tcga_id(col) for col in rna_data.columns]
    if isinstance(cna_data.columns, pd.Index):
        cna_data.columns = [standardize_tcga_id(col) for col in cna_data.columns]
    
    # Ensure samples are in columns (features in rows)
    if rna_data.shape[0] > rna_data.shape[1]:  # More rows than columns
        print("Transposing RNA data to have samples as columns")
        rna_data = rna_data.T
    if cna_data.shape[0] > cna_data.shape[1]:  # More rows than columns
        print("Transposing CNA data to have samples as columns")
        cna_data = cna_data.T
    
    # Get common samples
    common_samples = sorted(set(rna_data.columns) & set(cna_data.columns))
    print("Found {} common samples".format(len(common_samples)))
    print("First few common samples: {}".format(common_samples[:5]))
    
    if len(common_samples) == 0:
        # Try with index/columns swapped
        print("\nNo common samples found. Trying with transposed data...")
        rna_data = rna_data.T
        cna_data = cna_data.T
        common_samples = sorted(set(rna_data.columns) & set(cna_data.columns))
        print("After transposition: Found {} common samples".format(len(common_samples)))
        if len(common_samples) > 0:
            print("First few common samples after transposition: {}".format(common_samples[:5]))
    
    if len(common_samples) == 0:
        raise ValueError("No common samples found between RNA and CNA data")
    
    # Align data
    rna_aligned = rna_data[common_samples]
    cna_aligned = cna_data[common_samples]
    
    print("\nFinal aligned shapes:")
    print("RNA aligned: {}".format(rna_aligned.shape))
    print("CNA aligned: {}".format(cna_aligned.shape))
    
    return rna_aligned, cna_aligned

def get_predictions(model, X_data):
    """Get model predictions"""
    try:
        return model.predict(X_data)
    except Exception as e:
        print("Error during prediction: {}".format(str(e)))
        if hasattr(model, 'predict_proba'):
            print("Trying predict_proba instead...")
            return model.predict_proba(X_data)[:, 1]
        raise

def calculate_roc(y_true, y_pred):
    """Calculate ROC curve and AUC"""
    # Create binary labels using median as threshold
    threshold = np.median(y_true)
    y_binary = (y_true >= threshold).astype(int)
    
    fpr, tpr, _ = roc_curve(y_binary, y_pred)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc

def create_roc_plot(test_fpr, test_tpr, test_auc, val_fpr, val_tpr, val_auc, title, output_path):
    """Create and save ROC plot"""
    plt.figure(figsize=(8, 8))
    
    # Plot testing ROC curve
    plt.plot(test_fpr, test_tpr, color='red', lw=2,
             label='TCGA testing AUC: {:.3f}'.format(test_auc))
    
    # Plot validation ROC curve
    plt.plot(val_fpr, val_tpr, color='blue', lw=2,
             label='METABRIC validation AUC: {:.3f}'.format(val_auc))
    
    # Add diagonal line
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', alpha=0.8)
    
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title(title, fontsize=14, pad=20)
    plt.legend(loc="lower right", fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Load training/testing data
    print("Loading training/testing data...")
    try:
        rna_train = load_data("training_data/rna_signature_score_median_no_norm.pkl")
        cna_train = load_data("training_data/cna_segment_score_mean_no_norm.pkl")
        
        # Load validation data
        print("Loading validation data...")
        rna_val = load_data("validation_data/metabric_scores/rna_signature_score_cleaned.pkl")
        cna_val = load_data("validation_data/metabric_scores/cna_signature_score_cleaned.pkl")
    except Exception as e:
        print("Error loading data files: {}".format(str(e)))
        raise
    
    # Prepare training data
    print("\nPreparing training data...")
    rna_train_aligned, cna_train_aligned = prepare_data(rna_train, cna_train)
    
    # Prepare validation data
    print("\nPreparing validation data...")
    rna_val_aligned, cna_val_aligned = prepare_data(rna_val, cna_val)
    
    # Split training data into train/test (70/30 split)
    X_train, X_test, y_train, y_test = train_test_split(
        rna_train_aligned.values.T,  # Transpose to have samples as rows
        cna_train_aligned.values.T,  # Transpose to have samples as rows
        test_size=0.3,  # 70/30 split
        random_state=42
    )
    
    print("\nAfter train/test split:")
    print("X_train shape: {}".format(X_train.shape))
    print("X_test shape: {}".format(X_test.shape))
    print("y_train shape: {}".format(y_train.shape))
    print("y_test shape: {}".format(y_test.shape))
    
    # Scale the data
    scaler = StandardScaler()
    X_test = scaler.fit_transform(X_test)
    X_val = scaler.fit_transform(rna_val_aligned.values.T)  # Transpose to have samples as rows
    
    # Models to evaluate
    models = {
        'Basal Signaling': 'ROC_plots/models/basal_signaling_model.pkl',
        'Estrogen Signaling': 'ROC_plots/models/estrogen_signaling_model.pkl',
        'UNC RB LOH': 'ROC_plots/models/unc_rb_loh_model.pkl'
    }
    
    # Process each model
    for model_name, model_path in models.items():
        print("\nProcessing {}...".format(model_name))
        
        try:
            # Load model
            model = load_model(model_path)
            print("Model type: {}".format(type(model)))
            print("Model attributes: {}".format(dir(model)))
            
            # Get predictions
            test_pred = get_predictions(model, X_test)
            val_pred = get_predictions(model, X_val)
            
            print("Generated predictions:")
            print("Test predictions shape: {}".format(test_pred.shape))
            print("Validation predictions shape: {}".format(val_pred.shape))
            
            # Calculate ROC curves
            test_fpr, test_tpr, test_auc = calculate_roc(y_test[:, 0], test_pred)  # Using first CNA feature
            val_fpr, val_tpr, val_auc = calculate_roc(cna_val_aligned.values.T[:, 0], val_pred)  # Using first CNA feature
            
            # Create and save plot
            output_path = 'ROC_plots/{}_roc_comparison.png'.format(model_name.lower().replace(" ", "_"))
            create_roc_plot(
                test_fpr, test_tpr, test_auc,
                val_fpr, val_tpr, val_auc,
                '{} ROC Curves'.format(model_name),
                output_path
            )
            print("Plot saved as {}".format(output_path))
            
        except Exception as e:
            print("Error processing {}: {}".format(model_name, str(e)))
            continue

if __name__ == '__main__':
    main() 