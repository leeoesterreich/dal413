#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.patches import Rectangle

def parse_genomic_coordinate(feature_name):
    """Parse genomic coordinates from feature names"""
    # Handle whole chromosome arms
    if '.p.del' in feature_name or '.p.amp' in feature_name:
        chr_match = re.search(r'chr(\d+|X|Y)\.p\.', feature_name)
        if chr_match:
            chr_num = chr_match.group(1)
            return chr_num, 'p', 0.25  # Position for p-arm
    elif '.q.del' in feature_name or '.q.amp' in feature_name:
        chr_match = re.search(r'chr(\d+|X|Y)\.q\.', feature_name)
        if chr_match:
            chr_num = chr_match.group(1)
            return chr_num, 'q', 0.75  # Position for q-arm
    
    # Handle specific segments with coordinates
    coord_match = re.search(r'chr(\d+|X|Y):(\d+)-(\d+)', feature_name)
    if coord_match:
        chr_num = coord_match.group(1)
        start_pos = int(coord_match.group(2))
        end_pos = int(coord_match.group(3))
        
        # Estimate position on chromosome (simplified)
        # Chromosome lengths (approximate, in base pairs)
        chr_lengths = {
            '1': 249250621, '2': 242193529, '3': 198295559, '4': 190214555,
            '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
            '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
            '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
            '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
            '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
        }
        
        if chr_num in chr_lengths:
            chr_length = chr_lengths[chr_num]
            center_pos = (start_pos + end_pos) / 2
            relative_pos = center_pos / chr_length
            
            # Determine if it's on p-arm (< 0.5) or q-arm (>= 0.5)
            if relative_pos < 0.5:
                arm = 'p'
                arm_pos = relative_pos * 2  # Scale to 0-1 within p-arm
            else:
                arm = 'q'
                arm_pos = (relative_pos - 0.5) * 2  # Scale to 0-1 within q-arm
            
            return chr_num, arm, arm_pos
    
    return None, None, None

def plot_feature_landscape_from_excel(excel_file, sheet_name='GSEA Basal signaling'):
    """Plot feature landscape from Excel file weights"""
    
    # Read the Excel file
    print(f"Reading Excel file: {excel_file}")
    try:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        print(f"Successfully read sheet '{sheet_name}' with {len(df)} rows")
        print("Columns:", df.columns.tolist())
        print("\nFirst few rows:")
        print(df.head())
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        # Try to see what sheets are available
        try:
            xl_file = pd.ExcelFile(excel_file)
            print("Available sheets:", xl_file.sheet_names)
        except:
            pass
        return
    
    # Look for weight/coefficient column
    weight_cols = [col for col in df.columns if 'weight' in col.lower() or 'coeff' in col.lower() or 'value' in col.lower()]
    feature_cols = [col for col in df.columns if 'feature' in col.lower() or 'name' in col.lower() or col == df.columns[0]]
    
    print(f"Potential weight columns: {weight_cols}")
    print(f"Potential feature columns: {feature_cols}")
    
    if not weight_cols or not feature_cols:
        print("Available columns:")
        for i, col in enumerate(df.columns):
            print(f"{i}: {col}")
        return
    
    # Use the first weight and feature columns found
    weight_col = weight_cols[0]
    feature_col = feature_cols[0]
    
    print(f"Using weight column: {weight_col}")
    print(f"Using feature column: {feature_col}")
    
    # Filter non-zero weights
    non_zero_mask = df[weight_col] != 0
    features_with_weights = df[non_zero_mask]
    
    print(f"Found {len(features_with_weights)} non-zero weights")
    
    if len(features_with_weights) == 0:
        print("No non-zero weights found!")
        return
    
    # Parse genomic coordinates
    chromosomes = []
    arms = []
    positions = []
    weights = []
    feature_names = []
    
    for _, row in features_with_weights.iterrows():
        feature_name = str(row[feature_col])
        weight = float(row[weight_col])
        
        chr_num, arm, pos = parse_genomic_coordinate(feature_name)
        
        if chr_num is not None:
            chromosomes.append(chr_num)
            arms.append(arm)
            positions.append(pos)
            weights.append(weight)
            feature_names.append(feature_name)
    
    print(f"Successfully parsed coordinates for {len(chromosomes)} features")
    
    if len(chromosomes) == 0:
        print("No features with parseable genomic coordinates found!")
        print("Sample feature names:")
        for feature in features_with_weights[feature_col].head(10):
            print(f"  - {feature}")
        return
    
    # Create chromosome order
    chr_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    # Set up the plot
    fig, ax = plt.subplots(figsize=(20, 8))
    
    # Plot parameters
    bar_width = 0.35
    chr_positions = {}
    x_pos = 0
    
    # Create x-axis positions for chromosomes
    for chr_num in chr_order:
        chr_positions[chr_num] = x_pos
        x_pos += 1
    
    # Add background shading for chromosome arms
    for chr_num in chr_order:
        if chr_num in chr_positions:
            x_center = chr_positions[chr_num]
            # P-arm background (light gray)
            ax.add_patch(Rectangle((x_center - 0.4, -1), 0.35, 2, 
                                 facecolor='lightgray', alpha=0.3, zorder=0))
            # Q-arm background (slightly darker gray)
            ax.add_patch(Rectangle((x_center + 0.05, -1), 0.35, 2, 
                                 facecolor='gray', alpha=0.3, zorder=0))
    
    # Plot the weights
    colors = []
    x_positions = []
    y_values = []
    
    for i, (chr_num, arm, pos, weight) in enumerate(zip(chromosomes, arms, positions, weights)):
        if chr_num in chr_positions:
            # Calculate x position based on chromosome and arm
            chr_x = chr_positions[chr_num]
            if arm == 'p':
                x = chr_x - 0.4 + (pos * 0.35)  # Position within p-arm
            else:  # q-arm
                x = chr_x + 0.05 + (pos * 0.35)  # Position within q-arm
            
            x_positions.append(x)
            y_values.append(weight)
            
            # Color based on weight (positive = red, negative = blue)
            if weight > 0:
                colors.append('red')
            else:
                colors.append('blue')
    
    # Create scatter plot
    scatter = ax.scatter(x_positions, y_values, c=colors, alpha=0.7, s=50, zorder=2)
    
    # Customize the plot
    ax.set_xlabel('Chromosome', fontsize=14)
    ax.set_ylabel('Weight', fontsize=14)
    ax.set_title(f'Feature Landscape - {sheet_name} (from Excel)', fontsize=16, fontweight='bold')
    
    # Set x-axis
    ax.set_xticks([chr_positions[chr_num] for chr_num in chr_order if chr_num in chr_positions])
    ax.set_xticklabels([chr_num for chr_num in chr_order if chr_num in chr_positions])
    ax.set_xlim(-0.5, max(chr_positions.values()) + 0.5)
    
    # Add horizontal line at y=0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
                             markersize=8, label='Positive weight'),
                      Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
                             markersize=8, label='Negative weight')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Tight layout
    plt.tight_layout()
    
    # Save the plot
    output_file = f'{sheet_name.replace(" ", "_")}_excel_weights.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")
    
    # Show statistics
    print(f"\nWeight statistics:")
    print(f"Total features plotted: {len(y_values)}")
    print(f"Positive weights: {sum(1 for w in y_values if w > 0)}")
    print(f"Negative weights: {sum(1 for w in y_values if w < 0)}")
    print(f"Weight range: {min(y_values):.4f} to {max(y_values):.4f}")
    
    plt.show()

if __name__ == "__main__":
    excel_file = "41467_2019_13588_MOESM4_ESM_Elastic_Net_gene_signatures.xlsx"
    plot_feature_landscape_from_excel(excel_file, sheet_name='Sheet1') 