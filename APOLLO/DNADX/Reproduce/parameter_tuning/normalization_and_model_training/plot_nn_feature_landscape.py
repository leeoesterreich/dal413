#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse
import os

# hg19 chromosome lengths and p-arm lengths (approximations)
# Adapted from user's bar_plot_trained_model.py
CYTOBAND_HG19 = {
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

def parse_genomic_coordinates(feature_name_str, importance_net_val):
    """
    Parses genomic coordinates from a feature name string.
    Adapted from user's script. Uses importance_net_val as the 'weight'.
    Returns a dictionary with chr, start, end, absolute_start, absolute_end, and importance_net.
    Focuses on parsing chrN:start-end format, as neural network features are likely segments.
    """
    # Primary expected format: chrN:start-end (potentially with suffixes)
    segment_match = re.match(r'^(chr(\d+|X|Y)):(\d+)-(\d+)', feature_name_str, re.IGNORECASE)
    if segment_match:
        chrom_prefix = segment_match.group(1) # e.g. chr1, chrX
        chrom = segment_match.group(2)      # e.g. 1, X
        start = int(segment_match.group(3))
        end = int(segment_match.group(4))

        if chrom.upper() == 'Y' or chrom not in CHR_LENGTHS_HG19:
            return None

        return {
            'chr': chrom,
            'start_genomic': start,
            'end_genomic': end,
            'importance_net': importance_net_val,
            'feature_name': feature_name_str
        }
    # Fallback for other potential simple formats if needed, but prioritize above.
    # E.g., if some features are just "1p" or similar, they won't be parsed by this simplified version.
    # The NN features from the previous script were assumed to be direct segment names.
    print(f"Warning: Could not parse genomic coordinates from feature: {feature_name_str}")
    return None

def plot_nn_feature_landscape(csv_file_path):
    """Plots the neural network feature landscape from the importance CSV file."""
    if not os.path.exists(csv_file_path):
        print(f"Error: CSV file not found at {csv_file_path}")
        return

    df_importances = pd.read_csv(csv_file_path)
    if 'Segment' not in df_importances.columns or 'Importance_Net' not in df_importances.columns:
        print("Error: CSV must contain 'Segment' and 'Importance_Net' columns.")
        return

    parsed_features = []
    for _, row in df_importances.iterrows():
        parsed = parse_genomic_coordinates(row['Segment'], row['Importance_Net'])
        if parsed:
            parsed_features.append(parsed)
    
    if not parsed_features:
        print("No features could be parsed for plotting.")
        return

    # Calculate absolute genomic positions for plotting
    chr_offsets = {}
    current_offset = 0
    for chrom in ORDERED_CHROMS:
        chr_offsets[chrom] = current_offset
        current_offset += CHR_LENGTHS_HG19[chrom]
    genome_length_total = current_offset

    plot_data = []
    for feat in parsed_features:
        if feat['chr'] not in chr_offsets: # Should not happen if ORDERED_CHROMS is correct
            continue
        abs_start = chr_offsets[feat['chr']] + feat['start_genomic']
        abs_end = chr_offsets[feat['chr']] + feat['end_genomic']
        # For vlines, we typically plot at the midpoint of the segment
        abs_midpoint = (abs_start + abs_end) / 2
        plot_data.append({
            'abs_midpoint': abs_midpoint,
            'importance_net': feat['importance_net'],
            'chr': feat['chr']
        })
    
    df_plot = pd.DataFrame(plot_data)

    fig, ax = plt.subplots(figsize=(20, 6))

    # Plotting the bars (vertical lines for each segment)
    for _, row in df_plot.iterrows():
        color = 'red' if row['importance_net'] > 0 else 'blue'
        ax.vlines(x=row['abs_midpoint'], ymin=0, ymax=row['importance_net'], color=color, linewidth=0.5) # Thin lines for many segments

    # Set up chromosome labels and grid lines
    ax.set_xlim(0, genome_length_total)
    chr_label_positions = []
    chr_boundary_positions = [0]
    current_pos = 0
    for chrom in ORDERED_CHROMS:
        length = CHR_LENGTHS_HG19[chrom]
        chr_label_positions.append(current_pos + length / 2)
        current_pos += length
        chr_boundary_positions.append(current_pos)

    ax.set_xticks(chr_label_positions)
    ax.set_xticklabels(ORDERED_CHROMS)
    ax.tick_params(axis='x', which='major', length=0) # Remove major ticks on x-axis for cleaner look

    for boundary in chr_boundary_positions:
        ax.axvline(boundary, color='gray', linestyle='-', linewidth=0.75)
    
    # Y-axis limits and label
    max_abs_importance = df_plot['importance_net'].abs().max()
    if pd.isna(max_abs_importance) or max_abs_importance == 0:
        max_abs_importance = 0.5 # Default if no data or all zero
    ax.set_ylim(-max_abs_importance * 1.1, max_abs_importance * 1.1)
    ax.set_ylabel("Net Importance (1st Layer Weights)")
    ax.axhline(0, color='black', linewidth=0.75) # Line at y=0

    # Title
    plot_title = os.path.basename(csv_file_path).replace('_feature_weights.csv', '')
    ax.set_title(f"Feature Landscape for {plot_title} (NN Model)", fontsize=16)
    
    plt.tight_layout()
    
    # Save the plot
    output_png_path = csv_file_path.replace('_feature_weights.csv', '_landscape_plot.png')
    try:
        plt.savefig(output_png_path, dpi=300)
        print(f"Successfully saved plot to: {output_png_path}")
    except Exception as e:
        print(f"Error saving plot: {e}")
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot feature landscape from NN importance CSV.")
    parser.add_argument("csv_file_path", type=str, help="Path to the feature weights CSV file.")
    args = parser.parse_args()
    plot_nn_feature_landscape(args.csv_file_path) 