#!/usr/bin/env python3

import joblib
import pickle
import numpy as np

def test_load_model(model_file):
    """Test loading the model with different methods"""
    
    print(f"Testing load of: {model_file}")
    
    # Try method 1: joblib
    try:
        print("Trying joblib...")
        model_data = joblib.load(model_file)
        print(f"✓ Joblib success: {type(model_data)}")
        if isinstance(model_data, dict):
            print(f"  Keys: {list(model_data.keys())}")
            if 'model' in model_data:
                model = model_data['model']
                print(f"  Model type: {type(model)}")
                if hasattr(model, 'coef_'):
                    coef = model.coef_
                    print(f"  Coefficients shape: {coef.shape}")
                    print(f"  Non-zero coefficients: {np.sum(coef != 0)}")
            if 'feature_names' in model_data:
                print(f"  Feature names: {len(model_data['feature_names'])}")
        return model_data
    except Exception as e:
        print(f"✗ Joblib failed: {e}")
    
    # Try method 2: pickle
    try:
        print("Trying pickle...")
        with open(model_file, 'rb') as f:
            model_data = pickle.load(f)
        print(f"✓ Pickle success: {type(model_data)}")
        return model_data
    except Exception as e:
        print(f"✗ Pickle failed: {e}")
    
    return None

if __name__ == "__main__":
    model_file = "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP_model.pkl"
    model_data = test_load_model(model_file) 