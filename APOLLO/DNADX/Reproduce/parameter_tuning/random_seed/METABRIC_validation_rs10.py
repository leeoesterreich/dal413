import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from joblib import load
from sklearn.metrics import roc_curve, auc
from model_training_rs10 import load_data, build_prediction_model
from validate_model import clean_feature_names, impute_preserve_all_columns

def load_all_signatures():
    """Load all signatures from the signature_score.pkl file"""
    print("Loading signature scores...")
    signature_score = pd.read_pickle("/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/signature_score.pkl")
    all_signatures = signature_score.index.tolist()
    print("Found {} signatures in total".format(len(all_signatures)))
    return pd.DataFrame({'signature': all_signatures})

def train_and_validate_signatures():
    """Train and validate models for all signatures"""
    # Load all signatures
    all_signatures = load_all_signatures()
    
    # Load training data
    signature_score, segment_score = load_data()
    
    # Load METABRIC validation data
    metabric_signature_score, metabric_segment_score = load_metabric_data()
    
    # Create results directory
    os.makedirs('results_rs10/validation_summary', exist_ok=True)
    os.makedirs('results_rs10/models', exist_ok=True)
    os.makedirs('results_rs10/plots', exist_ok=True)
    os.makedirs('results_rs10/validation_plots', exist_ok=True)
    
    # Store results for visualization
    results = []
    
    # Process each signature
    for _, row in all_signatures.iterrows():
        signature = row['signature']
        print("\nProcessing signature: {}".format(signature))
        
        try:
            # Train model
            print("\nTraining model...")
            model, _, _, _ = build_prediction_model(signature_score, segment_score, signature)
            
            # Load the trained model
            safe_name = signature.replace('/', '_').replace('\\', '_')
            model_info = load('results_rs10/models/{}_model.joblib'.format(safe_name))
            
            # Validate model with outlier removal
            print("\nValidating model...")
            correlation, auc = validate_model_with_outlier_removal(
                model_info, metabric_signature_score, metabric_segment_score, signature
            )
            
            # Store results
            results.append({
                'signature': signature,
                'train_auc': row['train_auc'],
                'test_auc': row['test_auc'],
                'validation_auc': auc,
                'validation_correlation': correlation
            })
            
        except Exception as e:
            print("Error processing signature {}: {}".format(signature, str(e)))
            continue
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results
    results_df.to_csv('results_rs10/validation_summary/auc_comparison.csv', index=False)
    
    # Create visualization
    plot_auc_comparison(results_df)

if __name__ == "__main__":
    train_and_validate_signatures() 