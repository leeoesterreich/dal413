import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def safe_stats(data):
    """Calculate statistics safely handling empty arrays"""
    if len(data) == 0:
        return {
            'mean': 0,
            'std': 0,
            'median': 0,
            'q1': 0,
            'q3': 0,
            'iqr': 0
        }
    
    return {
        'mean': np.mean(data),
        'std': np.std(data),
        'median': np.median(data),
        'q1': np.percentile(data, 25),
        'q3': np.percentile(data, 75),
        'iqr': np.percentile(data, 75) - np.percentile(data, 25)
    }

def clean_and_normalize_data():
    print("Loading data...")
    # Load METABRIC data
    metabric_cna = pd.read_pickle("validation_data/metabric_scores/cna_signature_score.pkl")
    metabric_rna = pd.read_pickle("validation_data/metabric_scores/rna_signature_score.pkl")
    
    # Convert to float and handle NaN
    metabric_cna = metabric_cna.astype(float)
    metabric_rna = metabric_rna.astype(float)
    
    # Clean CNA data
    print("\nCleaning CNA data...")
    clean_cna = clean_dataset(metabric_cna, "CNA")
    
    # Clean RNA data
    print("\nCleaning RNA data...")
    clean_rna = clean_dataset(metabric_rna, "RNA")
    
    # Save cleaned data
    print("\nSaving cleaned data...")
    clean_cna.to_pickle("validation_data/metabric_scores/cna_signature_score_cleaned.pkl")
    clean_rna.to_pickle("validation_data/metabric_scores/rna_signature_score_cleaned.pkl")

def clean_dataset(data, data_type):
    """Clean a dataset by removing outliers using IQR method"""
    # Calculate statistics per feature
    print(f"\nCalculating {data_type} feature-wise statistics...")
    feature_stats = {}
    for feature in tqdm(data.index, desc=f"Computing {data_type} feature statistics"):
        feature_data = data.loc[feature]
        valid_data = feature_data[~np.isnan(feature_data) & ~np.isinf(feature_data)]
        feature_stats[feature] = safe_stats(valid_data)
    
    # Create clean dataset
    print(f"\nCleaning {data_type} outliers...")
    clean_data = data.copy()
    
    # Track modifications
    total_outliers = 0
    modified_features = []
    
    for feature in tqdm(data.index, desc=f"Processing {data_type} features"):
        stats = feature_stats[feature]
        lower = stats['q1'] - 3 * stats['iqr']
        upper = stats['q3'] + 3 * stats['iqr']
        
        # Count outliers
        feature_data = data.loc[feature]
        valid_mask = ~np.isnan(feature_data) & ~np.isinf(feature_data)
        if not np.any(valid_mask):
            continue
            
        outliers = np.sum((feature_data[valid_mask] < lower) | (feature_data[valid_mask] > upper))
        
        if outliers > 0:
            total_outliers += outliers
            modified_features.append((feature, outliers))
            
            # Cap values
            clean_data.loc[feature][valid_mask] = np.clip(feature_data[valid_mask], lower, upper)
    
    print(f"\nTotal {data_type} outliers capped: {total_outliers:,}")
    print(f"{data_type} features modified: {len(modified_features)}")
    
    # Print top 10 most affected features
    print(f"\nTop 10 {data_type} features with most outliers:")
    modified_features.sort(key=lambda x: x[1], reverse=True)
    for feature, count in modified_features[:10]:
        valid_mask = ~np.isnan(data.loc[feature]) & ~np.isinf(data.loc[feature])
        clean_valid_mask = ~np.isnan(clean_data.loc[feature]) & ~np.isinf(clean_data.loc[feature])
        
        if not np.any(valid_mask) or not np.any(clean_valid_mask):
            continue
            
        print(f"{feature}: {count:,} outliers")
        print(f"  Original range: [{np.min(data.loc[feature][valid_mask]):.3f}, {np.max(data.loc[feature][valid_mask]):.3f}]")
        print(f"  Cleaned range: [{np.min(clean_data.loc[feature][clean_valid_mask]):.3f}, {np.max(clean_data.loc[feature][clean_valid_mask]):.3f}]")
        print(f"  Stats: mean={feature_stats[feature]['mean']:.3f}, std={feature_stats[feature]['std']:.3f}")
    
    # Plot distributions before and after cleaning
    print(f"\nGenerating {data_type} visualization...")
    plt.figure(figsize=(15, 5))
    
    # Get valid values for plotting
    orig_values = data.values.flatten()
    orig_valid = orig_values[~np.isnan(orig_values) & ~np.isinf(orig_values)]
    
    clean_values = clean_data.values.flatten()
    clean_valid = clean_values[~np.isnan(clean_values) & ~np.isinf(clean_values)]
    
    plt.subplot(1, 2, 1)
    plt.hist(orig_valid, bins=50, range=(-2, 2), density=True)
    plt.title(f'Original {data_type} Distribution\n(limited to [-2, 2] range)')
    plt.xlabel('Value')
    plt.ylabel('Density')
    
    plt.subplot(1, 2, 2)
    plt.hist(clean_valid, bins=50, range=(-2, 2), density=True)
    plt.title(f'Cleaned {data_type} Distribution\n(limited to [-2, 2] range)')
    plt.xlabel('Value')
    plt.ylabel('Density')
    
    plt.tight_layout()
    plt.savefig(f'metabric_{data_type.lower()}_cleaning_results.png')
    plt.close()
    
    # Calculate and print statistics before and after cleaning
    print(f"\n{data_type} Data Statistics:")
    print("\nBefore cleaning:")
    print(f"Mean: {np.nanmean(orig_values):.3f}")
    print(f"Std: {np.nanstd(orig_values):.3f}")
    print(f"Min: {np.nanmin(orig_values):.3f}")
    print(f"Max: {np.nanmax(orig_values):.3f}")
    print(f"NaN %: {100 * np.isnan(orig_values).mean():.2f}%")
    
    print("\nAfter cleaning:")
    print(f"Mean: {np.nanmean(clean_values):.3f}")
    print(f"Std: {np.nanstd(clean_values):.3f}")
    print(f"Min: {np.nanmin(clean_values):.3f}")
    print(f"Max: {np.nanmax(clean_values):.3f}")
    print(f"NaN %: {100 * np.isnan(clean_values).mean():.2f}%")
    
    return clean_data

if __name__ == "__main__":
    clean_and_normalize_data() 