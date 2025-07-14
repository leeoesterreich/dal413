from joblib import load
import numpy as np

# Load the model
model_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results_validation/results/models/UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450_model.joblib"
model_info = load(model_path)

# Print model information
print("\nModel Information:")
print("-----------------")
print(f"Model type: {type(model_info['model'])}")
print(f"Number of features: {len(model_info['feature_names'])}")
print("\nFeature names (first 10):")
print(model_info['feature_names'][:10])
print("\nModel parameters:")
print(model_info['model'].get_params())

# If it's an ElasticNet model, print non-zero coefficients and their feature names
if hasattr(model_info['model'], 'coef_'):
    coefs = model_info['model'].coef_
    feature_names = model_info['feature_names']
    nonzero_indices = np.where(coefs != 0)[0]
    print(f"\nNumber of non-zero coefficients: {len(nonzero_indices)}")
    print("\nNon-zero coefficients and their feature names:")
    for idx in nonzero_indices:
        print(f"{feature_names[idx]}: {coefs[idx]}")
else:
    print("Model does not have coefficients.") 