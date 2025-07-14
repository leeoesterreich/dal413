# Placeholder for analyze_dataset_distributions.py 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mode as scipy_mode
import os
import traceback

def analyze_data(data_obj, obj_name, output_dir):
    """
    Analyzes a pandas DataFrame, Series, or NumPy array to print summary statistics 
    and plot value distribution of all numeric data it contains.
    """
    print(f"--- Analyzing: {obj_name} ---")

    all_values_list = []
    original_nan_count = 0
    object_shape = "N/A"
    object_size = "N/A"

    if isinstance(data_obj, pd.DataFrame):
        object_shape = data_obj.shape
        object_size = data_obj.size
        if data_obj.empty:
            print("DataFrame is empty. Skipping.")
            print("\n")
            return
        original_nan_count = data_obj.isnull().sum().sum()
        for col_name in data_obj.columns:
            try:
                numeric_series = pd.to_numeric(data_obj[col_name], errors='coerce')
                all_values_list.extend(numeric_series.dropna().tolist())
            except Exception as e:
                print(f"  Warning: Could not process column {col_name} in DataFrame {obj_name} to numeric: {e}")
    elif isinstance(data_obj, pd.Series):
        object_shape = data_obj.shape
        object_size = data_obj.size
        if data_obj.empty:
            print("Series is empty. Skipping.")
            print("\n")
            return
        original_nan_count = data_obj.isnull().sum()
        try:
            numeric_series = pd.to_numeric(data_obj, errors='coerce')
            all_values_list.extend(numeric_series.dropna().tolist())
        except Exception as e:
            print(f"  Warning: Could not process Series {obj_name} to numeric: {e}")
    elif isinstance(data_obj, np.ndarray):
        object_shape = data_obj.shape
        object_size = data_obj.size
        if data_obj.size == 0:
            print("NumPy array is empty. Skipping.")
            print("\n")
            return
        original_nan_count = np.isnan(data_obj).sum()
        try:
            # Ensure data is float before pd.to_numeric if it's a numpy array of objects or mixed types
            numeric_flat_array = pd.to_numeric(data_obj.astype(float).flatten(), errors='coerce')
            all_values_list.extend(pd.Series(numeric_flat_array).dropna().tolist())
        except Exception as e:
            print(f"  Warning: Could not process NumPy array {obj_name} to numeric: {e}")
    else:
        print(f"Object {obj_name} is of unhandled type: {type(data_obj)}. Skipping detailed analysis.")
        print("\n")
        return

    print(f"Original Object Shape: {object_shape}")
    print(f"Original Object Size: {object_size}")
    print(f"Total NaN values in original loaded data: {original_nan_count}")

    if not all_values_list:
        print("No numeric data found after attempting conversion and removing NaNs. Skipping statistics and plotting.\n")
        return
        
    all_values_np = np.array(all_values_list, dtype=float)
    print(f"Total non-NaN numeric data points for stats: {all_values_np.size}")

    if all_values_np.size == 0:
        print("No numeric data points available for statistics calculation.\n")
        return

    stats_mean = np.mean(all_values_np)
    stats_median = np.median(all_values_np)
    
    try:
        mode_result = scipy_mode(all_values_np, keepdims=False) 
    except TypeError: 
        mode_result = scipy_mode(all_values_np)

    stats_mode_val = np.nan
    if hasattr(mode_result, 'mode') and isinstance(mode_result.mode, np.ndarray) and mode_result.mode.size > 0:
        stats_mode_val = mode_result.mode[0] 
    elif hasattr(mode_result, 'mode') and not isinstance(mode_result.mode, np.ndarray): 
        stats_mode_val = mode_result.mode
    elif isinstance(mode_result, (int, float)): # Handle cases where mode_result itself is the mode
        stats_mode_val = mode_result

    
    stats_std = np.std(all_values_np)
    stats_min = np.min(all_values_np)
    stats_max = np.max(all_values_np)

    print(f"Mean of all numeric values: {stats_mean:.4f}")
    print(f"Median of all numeric values: {stats_median:.4f}")
    print(f"Mode of all numeric values: {stats_mode_val if pd.notna(stats_mode_val) else 'N/A'}")
    print(f"Standard Deviation of all numeric values: {stats_std:.4f}")
    print(f"Min of all numeric values: {stats_min:.4f}")
    print(f"Max of all numeric values: {stats_max:.4f}")

    plt.figure(figsize=(10, 6))
    # Filter out potential NaNs or Infs that might have slipped through, for plotting
    finite_values_for_plot = all_values_np[np.isfinite(all_values_np)]
    if finite_values_for_plot.size > 0:
        plt.hist(finite_values_for_plot, bins=50, edgecolor='black')
    else:
        print("No finite values to plot for histogram.")
        
    safe_obj_name_for_title = obj_name.replace('_', ' ').title()
    safe_obj_name_for_file = obj_name.replace(' ', '_').replace('/', '_').replace(':', '_')

    plt.title(f'Distribution of All Numeric Values in {safe_obj_name_for_title}')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(True, alpha=0.5)
    
    plot_filename = os.path.join(output_dir, f"{safe_obj_name_for_file}_overall_distribution.png")
    try:
        plt.savefig(plot_filename)
        print(f"Overall distribution plot saved to: {plot_filename}")
    except Exception as e:
        print(f"Error saving plot {plot_filename}: {e}")
    plt.close()
    print("\n")


def main():
    data_paths = {
        "Recalculated_Training_CNA_Scores_V2_Fixed": "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/recalculated_scores_raw_inputs_v2_fixed_entrez/training_cna_segment_scores_recalculated.pkl",
        "Recalculated_Training_RNA_Scores_V2_Fixed": "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/recalculated_scores_raw_inputs_v2_fixed_entrez/training_rna_signature_scores_recalculated.pkl",
        "Recalculated_Validation_CNA_Scores_V2_Fixed": "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/recalculated_scores_raw_inputs_v2_fixed_entrez/validation_cna_segment_scores_recalculated.pkl",
        "Recalculated_Validation_RNA_Scores_V2_Fixed": "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/recalculated_scores_raw_inputs_v2_fixed_entrez/validation_rna_signature_scores_recalculated.pkl"
    }

    output_plot_dir = "recalculated_scores_distributions_v2_fixed_entrez"
    os.makedirs(output_plot_dir, exist_ok=True)
    print(f"Output plots will be saved to: {os.path.abspath(output_plot_dir)}\n")

    for name, path in data_paths.items():
        print(f"Loading dataset: {name} from {path}")
        try:
            data_object = pd.read_pickle(path)
            analyze_data(data_object, name, output_plot_dir)
        except FileNotFoundError:
            print(f"File not found: {path}. Skipping.\n")
        except Exception as e:
            print(f"General error loading or processing {name} from {path}: {e}\n")
            traceback.print_exc()
            
if __name__ == "__main__":
    main() 