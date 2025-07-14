import pandas as pd
import numpy as np
from scipy import stats

def print_signature_statistics(data_obj, signatures):
    print("\nAnalyzing specific signatures:")
    if data_obj is None:
        print("  Data object is None.")
        return

    if isinstance(data_obj, pd.DataFrame):
        print(f"  DataFrame Shape: {data_obj.shape}")
        if data_obj.empty:
            print("  DataFrame is empty.")
            return
        
        for signature in signatures:
            if signature in data_obj.index:
                col_data = data_obj.loc[signature].dropna()
                if col_data.empty:
                    print(f"    Signature '{signature}': All NaN or empty after dropping NaNs.")
                    continue
                
                col_finite = col_data[np.isfinite(col_data)]
                if col_finite.empty:
                    print(f"    Signature '{signature}': All non-finite (inf/NaN) or empty after filtering non-finite.")
                    continue

                try:
                    mode_res = stats.mode(col_finite.values)
                    actual_mode = mode_res.mode[0] if isinstance(mode_res.mode, np.ndarray) and mode_res.mode.size > 0 else mode_res.mode
                    actual_count = mode_res.count[0] if isinstance(mode_res.count, np.ndarray) and mode_res.count.size > 0 else mode_res.count
                    
                    print(f"\n    Signature: {signature}")
                    print(f"    Minimum: {col_finite.min():.4f}")
                    print(f"    Maximum: {col_finite.max():.4f}")
                    print(f"    Mean: {col_finite.mean():.4f}")
                    print(f"    Median: {col_finite.median():.4f}")
                    print(f"    Mode: {actual_mode:.4f} (Count: {actual_count})")
                    print(f"    Standard Deviation: {col_finite.std():.4f}")
                    print(f"    Number of non-null values: {len(col_finite)}")
                    
                    # Calculate quartiles
                    q1, q3 = np.percentile(col_finite, [25, 75])
                    print(f"    Q1 (25th percentile): {q1:.4f}")
                    print(f"    Q3 (75th percentile): {q3:.4f}")
                    print(f"    IQR (Interquartile Range): {(q3 - q1):.4f}")
                    
                except Exception as e:
                    print(f"    Error calculating statistics for signature '{signature}'. Error: {e}")
            else:
                print(f"\n    Signature '{signature}' not found in the dataset.")

def main():
    rna_file_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/training_data/rna_signature_score_median_no_norm.pkl'
    
    target_signatures = [
        'GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP',
        'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN'
    ]

    print(f"Loading RNA data from: {rna_file_path}")
    try:
        rna_data = pd.read_pickle(rna_file_path)
        if any(sig in rna_data.columns for sig in target_signatures):
            print("Signatures found in columns, using as is")
        else:
            print("Transposing data to check for signatures in rows")
            rna_data = rna_data.T
        print_signature_statistics(rna_data, target_signatures)
    except FileNotFoundError:
        print(f"Error: RNA file not found at {rna_file_path}")
    except Exception as e:
        print(f"Error loading or processing RNA file: {e}")

if __name__ == '__main__':
    main() 