#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re
import argparse
import os
import xgboost
import torch # For loading .pt metadata

# hg19 chromosome lengths and p-arm lengths (approximations)
# Adapted from user's bar_plot_trained_model.py and plot_nn_feature_landscape.py
CYTOBAND_HG19 = { # p-arm lengths, used to define centromere start for q-arm (not directly used in this script yet)
    '1': 125000000, '2': 93300000, '3': 91000000, '4': 50400000, '5': 48400000,
    '6': 61000000, '7': 59900000, '8': 45600000, '9': 49000000, '10': 40200000,
    '11': 53700000, '12': 35800000, '13': 17900000, '14': 17600000, '15': 19000000,
    '16': 36600000, '17': 24000000, '18': 17200000, '19': 26500000, '20': 27500000,
    '21': 13200000, '22': 14700000, 'X': 60600000
}

CHR_LENGTHS_HG19 = {
    '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260,
    '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747,
    '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
    '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
    '21': 48129895, '22': 51304566, 'X': 155270560
}
ORDERED_CHROMS = [str(i) for i in range(1, 23)] + ['X']

# Summary CSV and output directory constants
SUMMARY_CSV_PATH = 'results_xgb_retrained_60_20_20_split/summaries/all_signatures_retrained_results.csv'
DEFAULT_OUTPUT_DIR = 'results_xgb_retrained_60_20_20_split/feature_landscape_plots/'


def sanitize_filename(name):
    """Sanitizes a string to be a valid filename part."""
    name = name.replace('.r=', '_r_') # Specific replacements first
    name = re.sub(r'[^a-zA-Z0-9_.-]', '_', name)
    return name

def parse_feature_to_coords(feature_name_str):
    """
    Parses genomic coordinates from a feature name string.
    Prioritizes 'chrN:start-end' format.
    Returns a dictionary like {'chr': chr_num, 'start_genomic': start, 'end_genomic': end} or None.
    """
    # More robust: try to find chrN:start-end pattern anywhere in the string first
    search_match = re.search(r'(chr(\d+|X|Y)):(\d+)-(\d+)', feature_name_str, re.IGNORECASE)
    if search_match:
        # group(1) is like 'chr1', group(2) is '1', group(3) is start, group(4) is end
        chrom = search_match.group(2)
        start = int(search_match.group(3))
        end = int(search_match.group(4))
    else:
        # Fallback: try matching if the string STARTS with (optional chr)N:start-end
        segment_match_strict = re.match(r'^(?:chr)?(\d+|X|Y):(\d+)-(\d+)', feature_name_str, re.IGNORECASE)
        if segment_match_strict:
            chrom = segment_match_strict.group(1)
            start = int(segment_match_strict.group(2))
            end = int(segment_match_strict.group(3))
        else:
            # print(f"Warning: Could not parse genomic coordinates from feature: {feature_name_str}.")
            return None

    if chrom.upper() == 'Y' or chrom not in CHR_LENGTHS_HG19:
        # print(f"Warning: Chromosome {chrom} from feature {feature_name_str} is not recognized or is Y. Skipping.")
        return None
    return {'chr': chrom, 'start_genomic': start, 'end_genomic': end}

def get_model_data_paths(signature_name, summary_csv_path):
    """Reads the summary CSV to find model and metadata file paths for a signature."""
    try:
        summary_df = pd.read_csv(summary_csv_path)
    except FileNotFoundError:
        print(f"Error: Summary CSV file not found at {summary_csv_path}")
        return None, None

    signature_row = summary_df[summary_df['Signature'] == signature_name]
    if signature_row.empty:
        print(f"Error: Signature '{signature_name}' not found in {summary_csv_path}")
        return None, None
    
    model_file = signature_row['Model_File'].iloc[0]
    metadata_file = signature_row['Metadata_File'].iloc[0]
    
    if not model_file or pd.isna(model_file) or not os.path.exists(model_file):
        print(f"Error: Model file path not found, invalid, or file does not exist for {signature_name}: {model_file}")
        return None, None
    if not metadata_file or pd.isna(metadata_file) or not os.path.exists(metadata_file):
        print(f"Error: Metadata file path not found, invalid, or file does not exist for {signature_name}: {metadata_file}")
        return None, None
        
    return model_file, metadata_file

def plot_landscape(signature_name, model_json_path, metadata_pt_path, output_dir, importance_type, log_scale_y):
    """Generates and saves the feature importance landscape plot for an XGBoost model."""
    processing_message = f"Processing signature: {signature_name} (Importance type: {importance_type}"
    if log_scale_y:
        processing_message += ", Y-axis: Log10 scale"
    processing_message += ")"
    print(processing_message)
    try:
        xgb_model = xgboost.Booster()
        xgb_model.load_model(model_json_path)
        print(f"  Loaded model from: {model_json_path}")
    except Exception as e:
        print(f"Error loading XGBoost model {model_json_path}: {e}")
        return

    try:
        metadata = torch.load(metadata_pt_path, map_location='cpu')
        feature_names = metadata.get('feature_names')
        if feature_names is None:
            print(f"Error: 'feature_names' not found in metadata file {metadata_pt_path}")
            return
        print(f"  Loaded metadata from: {metadata_pt_path}, found {len(feature_names)} feature names.")
    except Exception as e:
        print(f"Error loading metadata file {metadata_pt_path}: {e}")
        return

    importances_raw = xgb_model.get_score(importance_type=importance_type)
    
    # Map generic 'f_idx' to actual feature names
    feature_importance_map = {}
    for f_idx_str, score in importances_raw.items():
        try:
            actual_idx = int(f_idx_str[1:]) # remove 'f' and convert to int
            if actual_idx < len(feature_names):
                feature_importance_map[feature_names[actual_idx]] = score
        except ValueError:
            print(f"Warning: Could not parse feature index {f_idx_str}")
            continue
    
    if not feature_importance_map:
        print(f"No feature importances ({importance_type}) found in the model or could not map to feature names.")
        #return # Continue to plot an empty graph if desired, or return

    parsed_plot_features = []
    unparsed_count = 0
    for i, name in enumerate(feature_names):
        coords = parse_feature_to_coords(name)
        if coords:
            importance = feature_importance_map.get(name, 0.0) # Default to 0 if not in map (i.e., gain was 0)
            parsed_plot_features.append({
                'chr': coords['chr'],
                'start_genomic': coords['start_genomic'],
                'end_genomic': coords['end_genomic'],
                'importance': importance,
                'name': name
            })
        else:
            unparsed_count +=1
            if unparsed_count < 10: # Print only a few warnings
                 print(f"Warning: Could not parse feature '{name}' to coordinates. Skipping this feature.")
            elif unparsed_count == 10:
                 print("Further parsing warnings will be suppressed...")


    if unparsed_count > 0:
        print(f"Total unparsed features: {unparsed_count} out of {len(feature_names)}")
    
    if not parsed_plot_features:
        print(f"No features could be parsed and mapped for signature {signature_name}. Plot will be empty or not generated.")
        # Optionally, create an empty plot or return
        # For now, let's allow an attempt to plot even if empty, matplotlib might handle it or error.

    # Calculate absolute genomic positions for plotting
    chr_offsets = {}
    current_offset = 0
    for chrom in ORDERED_CHROMS:
        chr_offsets[chrom] = current_offset
        current_offset += CHR_LENGTHS_HG19.get(chrom, 0) # Use .get for safety
    genome_length_total = current_offset

    plot_data_for_df = []
    for feat in parsed_plot_features:
        if feat['chr'] not in chr_offsets or feat['chr'] not in CHR_LENGTHS_HG19:
            print(f"Warning: Chromosome {feat['chr']} for feature {feat['name']} not in offset/length map. Skipping.")
            continue
        
        abs_start = chr_offsets[feat['chr']] + feat['start_genomic']
        abs_end = chr_offsets[feat['chr']] + feat['end_genomic']
        abs_midpoint = (abs_start + abs_end) / 2
        plot_data_for_df.append({
            'abs_midpoint': abs_midpoint,
            'importance': feat['importance'],
            'chr': feat['chr']
        })
    
    df_plot = pd.DataFrame(plot_data_for_df)

    fig, ax = plt.subplots(figsize=(25, 7)) # Wider figure

    y_plot_column = 'importance'
    y_axis_label = f"XGBoost Feature Importance ({importance_type.capitalize()})"

    if log_scale_y:
        if not df_plot.empty and df_plot['importance'].gt(0).any():
            # Filter out non-positive values before log transform
            df_plot_positive = df_plot[df_plot['importance'] > 0].copy()
            df_plot_positive['importance_log10'] = np.log10(df_plot_positive['importance'])
            y_plot_column = 'importance_log10'
            y_axis_label = f"XGBoost Feature Importance ({importance_type.capitalize()}) - Log10 Scale"
            df_plot = df_plot_positive # Use the dataframe with log values
        else:
            print("Warning: No positive importance values to plot on log scale. Using linear scale or plot will be empty.")
            # Keep y_plot_column = 'importance' and linear scale label if all are zero or less
            # Or, if df_plot becomes empty, plotting logic below handles it.

    # Plotting vertical lines for each segment's importance
    if not df_plot.empty and y_plot_column in df_plot.columns:
        ax.vlines(df_plot['abs_midpoint'], ymin=0, ymax=df_plot[y_plot_column], color='steelblue', linewidth=0.6, alpha=0.7)
    else:
        print("DataFrame for plotting is empty or lacks y_plot_column. No lines will be drawn.")

    # Set up chromosome labels and boundary lines
    ax.set_xlim(0, genome_length_total)
    chr_label_positions = []
    chr_boundary_positions = [0]
    current_pos = 0
    for chrom in ORDERED_CHROMS:
        length = CHR_LENGTHS_HG19.get(chrom, 0)
        if length == 0: continue
        chr_label_positions.append(current_pos + length / 2)
        current_pos += length
        chr_boundary_positions.append(current_pos)

    ax.set_xticks(chr_label_positions)
    ax.set_xticklabels(ORDERED_CHROMS, fontsize=9)
    ax.tick_params(axis='x', which='major', length=0) 

    for boundary in chr_boundary_positions:
        ax.axvline(boundary, color='darkgray', linestyle='-', linewidth=0.8)
    
    # Y-axis limits and label
    if not df_plot.empty and y_plot_column in df_plot.columns and df_plot[y_plot_column].notna().any():
        max_val = df_plot[y_plot_column].max()
        min_val = df_plot[y_plot_column].min()
        if log_scale_y and y_plot_column == 'importance_log10':
            # For log scale, ymin=0 might not make sense if all log values are far from 0.
            # Let's set ymin to be slightly below the smallest log value, or 0 if smallest is positive.
            plot_ymin = min(0, min_val - 0.1 * abs(min_val)) if min_val < 0 else 0
            ax.set_ylim(plot_ymin, max_val + 0.1 * abs(max_val) if max_val != 0 else 1) 
        else: # Linear scale
            ax.set_ylim(0, max_val * 1.1 if max_val > 0 else 1) 
    else:
        ax.set_ylim(0, 1) # Default if no data

    ax.set_ylabel(y_axis_label, fontsize=12)
    ax.axhline(0, color='black', linewidth=0.75)

    # Title
    signature_name_display = signature_name[:70] + '...' if len(signature_name) > 70 else signature_name
    ax.set_title(f"XGBoost Feature Importance ({importance_type.capitalize()}) Landscape for\n{signature_name_display}", fontsize=14, pad=10)
    
    plt.tight_layout()
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    sanitized_sig_name = sanitize_filename(signature_name)
    filename_suffix = f"_xgb_{importance_type}"
    if log_scale_y:
        filename_suffix += "_log10"
    filename_suffix += "_landscape.png"
    output_png_path = os.path.join(output_dir, f"{sanitized_sig_name}{filename_suffix}")
    
    try:
        plt.savefig(output_png_path, dpi=300)
        print(f"Successfully saved plot to: {output_png_path}")
    except Exception as e:
        print(f"Error saving plot: {e}")
    plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot XGBoost feature importance (gain) landscape.")
    parser.add_argument("signature_name", type=str, help="Full signature name as it appears in the summary CSV.")
    parser.add_argument("--summary_csv", type=str, default=SUMMARY_CSV_PATH, help="Path to the summary CSV file.")
    parser.add_argument("--output_dir", type=str, default=DEFAULT_OUTPUT_DIR, help="Directory to save the output plot.")
    parser.add_argument("--importance_type", type=str, default='gain', choices=['gain', 'weight', 'cover'], help="Type of importance to plot (gain, weight, or cover).")
    parser.add_argument("--log_scale_y", action='store_true', help="Apply log10 scale to the y-axis for importance values.")
    
    args = parser.parse_args()

    model_file, metadata_file = get_model_data_paths(args.signature_name, args.summary_csv)
    
    if model_file and metadata_file:
        plot_landscape(args.signature_name, model_file, metadata_file, args.output_dir, args.importance_type, args.log_scale_y)
    else:
        print(f"Could not proceed with plotting for signature: {args.signature_name} due to missing file paths.") 