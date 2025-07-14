import pandas as pd
import numpy as np
from scipy import stats

def print_detailed_statistics(data_obj, name):
    print(f"\nDescriptive statistics for: {name}")
    if data_obj is None:
        print("  Data object is None.")
        return

    if isinstance(data_obj, pd.DataFrame):
        print(f"  DataFrame Shape: {data_obj.shape}")
        if data_obj.empty:
            print("  DataFrame is empty.")
            return
        # Try to infer if features/signatures are rows or columns
        # If number of rows > number of columns by a large margin, assume columns are features/signatures
        # This is a heuristic, might need adjustment based on actual data structure
        data_to_iterate = data_obj
        if data_obj.shape[0] > data_obj.shape[1] * 2 and data_obj.shape[1] < 1000: # Heuristic for features as columns
             print(f"  Iterating over {data_obj.shape[1]} columns (assumed to be features/signatures).")
        elif data_obj.shape[1] > data_obj.shape[0] * 2 and data_obj.shape[0] < 1000: # Heuristic for features as rows
            print(f"  Data has more columns than rows. Transposing to iterate over {data_obj.shape[0]} rows (assumed to be features/signatures).")
            data_to_iterate = data_obj.T
        else:
            print(f"  Could not clearly determine feature/signature orientation. Assuming columns are features/signatures (Iterating over {data_obj.shape[1]} columns).")

        for col_name in data_to_iterate.columns:
            col_data = data_to_iterate[col_name].dropna()
            if col_data.empty:
                print(f"    Column '{col_name}': All NaN or empty after dropping NaNs.")
                continue
            if not pd.api.types.is_numeric_dtype(col_data):
                print(f"    Column '{col_name}': Non-numeric data. Skipping stats.")
                continue
            
            col_finite = col_data[np.isfinite(col_data)]
            if col_finite.empty:
                print(f"    Column '{col_name}': All non-finite (inf/NaN) or empty after filtering non-finite.")
                continue

            try:
                mode_res = stats.mode(col_finite.values) # Pass numpy array to stats.mode
                actual_mode = mode_res.mode[0] if isinstance(mode_res.mode, np.ndarray) and mode_res.mode.size > 0 else mode_res.mode
                actual_count = mode_res.count[0] if isinstance(mode_res.count, np.ndarray) and mode_res.count.size > 0 else mode_res.count
                print(f"    Column '{col_name}': Min={col_finite.min():.4f}, Max={col_finite.max():.4f}, Mean={col_finite.mean():.4f}, Median={col_finite.median():.4f}, Mode={actual_mode:.4f} (Count: {actual_count})")
            except Exception as e:
                print(f"    Could not calculate mode for column '{col_name}'. Error: {e}. Stats without mode: Min={col_finite.min():.4f}, Max={col_finite.max():.4f}, Mean={col_finite.mean():.4f}, Median={col_finite.median():.4f}")

    elif isinstance(data_obj, pd.Series):
        print(f"  Series Name: {data_obj.name}, Shape: {data_obj.shape}")
        if data_obj.empty:
            print("  Series is empty.")
            return
        
        col_data = data_obj.dropna()
        if col_data.empty:
            print(f"    Series '{data_obj.name}': All NaN or empty after dropping NaNs.")
            return
        if not pd.api.types.is_numeric_dtype(col_data):
            print(f"    Series '{data_obj.name}': Non-numeric data. Skipping stats.")
            return

        col_finite = col_data[np.isfinite(col_data)]
        if col_finite.empty:
            print(f"    Series '{data_obj.name}': All non-finite (inf/NaN) or empty after filtering non-finite.")
            return
        
        try:
            mode_res = stats.mode(col_finite.values) # Pass numpy array
            actual_mode = mode_res.mode[0] if isinstance(mode_res.mode, np.ndarray) and mode_res.mode.size > 0 else mode_res.mode
            actual_count = mode_res.count[0] if isinstance(mode_res.count, np.ndarray) and mode_res.count.size > 0 else mode_res.count
            print(f"    Stats for '{data_obj.name}': Min={col_finite.min():.4f}, Max={col_finite.max():.4f}, Mean={col_finite.mean():.4f}, Median={col_finite.median():.4f}, Mode={actual_mode:.4f} (Count: {actual_count})")
        except Exception as e:
            print(f"    Could not calculate mode for series '{data_obj.name}'. Error: {e}. Stats without mode: Min={col_finite.min():.4f}, Max={col_finite.max():.4f}, Mean={col_finite.mean():.4f}, Median={col_finite.median():.4f}")
    else:
        print(f"  Data is of an unexpected type: {type(data_obj)}. Cannot calculate detailed stats.")


def main():
    cna_file_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/training_data/cna_segment_score_mean_no_norm.pkl'
    rna_file_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/training_data/rna_signature_score_median_no_norm.pkl'

    print(f"Loading CNA data from: {cna_file_path}")
    try:
        cna_data = pd.read_pickle(cna_file_path)
        print_detailed_statistics(cna_data, "CNA Scores (cna_segment_score_mean_no_norm.pkl)")
    except FileNotFoundError:
        print(f"Error: CNA file not found at {cna_file_path}")
    except Exception as e:
        print(f"Error loading or processing CNA file: {e}")

    print(f"\nLoading RNA data from: {rna_file_path}")
    try:
        rna_data = pd.read_pickle(rna_file_path)
        print_detailed_statistics(rna_data, "RNA Scores (rna_signature_score_median_no_norm.pkl)")
    except FileNotFoundError:
        print(f"Error: RNA file not found at {rna_file_path}")
    except Exception as e:
        print(f"Error loading or processing RNA file: {e}")

if __name__ == '__main__':
    main() 