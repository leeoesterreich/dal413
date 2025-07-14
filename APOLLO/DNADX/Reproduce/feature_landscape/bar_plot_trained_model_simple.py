#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import joblib
from matplotlib.patches import Rectangle

def plot_trained_model_bars(model_file):
    """Plot feature landscape for trained model using bar plot"""
    
    print(f"Loading model file: {model_file}")
    
    # Load the model with joblib
    try:
        model_data = joblib.load(model_file)
        print(f"Model data type: {type(model_data)}")
        print(f"Model data keys: {model_data.keys() if isinstance(model_data, dict) else 'Not a dict'}")
        
        # Extract model and feature names
        model = model_data.get('model')
        feature_names = model_data.get('feature_names', [])
        
        print(f"Model type: {type(model)}")
        print(f"Number of feature names: {len(feature_names)}")
        
        # Get coefficients
        coefficients = model.coef_
        if len(coefficients.shape) > 1:
            coefficients = coefficients[0]  # Take first row if 2D
        
        print(f"Number of coefficients: {len(coefficients)}")
        print(f"Non-zero coefficients: {np.sum(coefficients != 0)}")
        
    except Exception as e:
        print(f"Error loading model: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Find non-zero weights
    non_zero_mask = coefficients != 0
    non_zero_features = [feature_names[i] for i in range(len(feature_names)) if non_zero_mask[i]]
    non_zero_weights = coefficients[non_zero_mask]
    
    print(f"Non-zero features: {len(non_zero_weights)}")
    
    if len(non_zero_weights) == 0:
        print("No non-zero weights found!")
        return
    
    # Parse genomic coordinates and plot
    chromosomes = []
    arms = []
    positions = []
    final_weights = []
    final_feature_names = []
    
    chr_lengths = {
        '1': 249250621, '2': 242193529, '3': 198295559, '4': 190214555,
        '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
        '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
        '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
        '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
        '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
    }
    
    for feature, weight in zip(non_zero_features, non_zero_weights):
        feature_str = str(feature)
        
        # Handle whole chromosome arms
        if '.p' in feature_str and 'wholearm' in feature_str:
            chr_match = re.search(r'(\d+|X|Y)\.p', feature_str)
            if chr_match:
                chr_num = chr_match.group(1)
                chromosomes.append(chr_num)
                arms.append('p')
                positions.append(0.25)
                final_weights.append(weight)
                final_feature_names.append(feature_str)
        elif '.q' in feature_str and 'wholearm' in feature_str:
            chr_match = re.search(r'(\d+|X|Y)\.q', feature_str)
            if chr_match:
                chr_num = chr_match.group(1)
                chromosomes.append(chr_num)
                arms.append('q')
                positions.append(0.75)
                final_weights.append(weight)
                final_feature_names.append(feature_str)
        else:
            # Handle specific segments
            coord_match = re.search(r'chr(\d+|X|Y):(\d+)-(\d+)', feature_str)
            if coord_match:
                chr_num = coord_match.group(1)
                start_pos = int(coord_match.group(2))
                end_pos = int(coord_match.group(3))
                
                if chr_num in chr_lengths:
                    chr_length = chr_lengths[chr_num]
                    center_pos = (start_pos + end_pos) / 2
                    relative_pos = center_pos / chr_length
                    
                    if relative_pos < 0.5:
                        arm = 'p'
                        arm_pos = relative_pos * 2
                    else:
                        arm = 'q'
                        arm_pos = (relative_pos - 0.5) * 2
                    
                    chromosomes.append(chr_num)
                    arms.append(arm)
                    positions.append(arm_pos)
                    final_weights.append(weight)
                    final_feature_names.append(feature_str)
    
    print(f"Successfully parsed {len(chromosomes)} features with coordinates")
    
    if len(chromosomes) == 0:
        print("No features with parseable coordinates!")
        print("Sample feature names:")
        for i, feature in enumerate(non_zero_features[:10]):
            print(f"  {i}: {feature}")
        return
    
    # Sort features by chromosome and position for better visualization
    chr_order_map = {str(i): i for i in range(1, 23)}
    chr_order_map.update({'X': 23, 'Y': 24})
    
    # Create sorting key: chromosome number * 100 + arm (p=0, q=50) + position
    sort_keys = []
    for chr_num, arm, pos in zip(chromosomes, arms, positions):
        chr_idx = chr_order_map.get(chr_num, 999)
        arm_offset = 0 if arm == 'p' else 50
        sort_key = chr_idx * 100 + arm_offset + pos
        sort_keys.append(sort_key)
    
    # Sort all arrays by the sort key
    sorted_indices = np.argsort(sort_keys)
    chromosomes = [chromosomes[i] for i in sorted_indices]
    arms = [arms[i] for i in sorted_indices]
    positions = [positions[i] for i in sorted_indices]
    final_weights = [final_weights[i] for i in sorted_indices]
    final_feature_names = [final_feature_names[i] for i in sorted_indices]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(25, 10))
    
    # Create x positions for bars
    x_positions = np.arange(len(final_weights))
    
    # Create colors for bars
    colors = ['red' if w > 0 else 'blue' for w in final_weights]
    
    # Create bar plot
    bars = ax.bar(x_positions, final_weights, color=colors, alpha=0.7, width=0.8)
    
    # Customize plot
    ax.set_xlabel('Features (ordered by genomic position)', fontsize=14)
    ax.set_ylabel('Weight', fontsize=14)
    ax.set_title('GSEA Basal Signaling Feature Landscape (Bar Plot from Trained Model)', fontsize=16, fontweight='bold')
    
    # Add chromosome labels
    chr_boundaries = []
    current_chr = None
    for i, chr_num in enumerate(chromosomes):
        if chr_num != current_chr:
            chr_boundaries.append(i)
            current_chr = chr_num
    chr_boundaries.append(len(chromosomes))
    
    # Add chromosome labels at boundaries
    for i in range(len(chr_boundaries) - 1):
        start_idx = chr_boundaries[i]
        end_idx = chr_boundaries[i + 1]
        middle_x = (start_idx + end_idx - 1) / 2
        chr_name = chromosomes[start_idx]
        ax.text(middle_x, ax.get_ylim()[1] * 0.95, f'Chr{chr_name}', 
                ha='center', va='top', fontsize=10, fontweight='bold')
    
    # Add vertical lines between chromosomes
    for boundary in chr_boundaries[1:-1]:
        ax.axvline(x=boundary - 0.5, color='gray', linestyle='--', alpha=0.5)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Remove x-axis tick labels (too many features to show)
    ax.set_xticks([])
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', alpha=0.7, label='Positive weight'),
                      Patch(facecolor='blue', alpha=0.7, label='Negative weight')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Generate output filename based on input
    output_name = model_file.split('/')[-1].replace('.pkl', '').replace('.joblib', '')
    output_file = f'{output_name}_trained_bars.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")
    
    # Print statistics
    print(f"\nStatistics:")
    print(f"Total features plotted: {len(final_weights)}")
    print(f"Positive weights: {sum(1 for w in final_weights if w > 0)}")
    print(f"Negative weights: {sum(1 for w in final_weights if w < 0)}")
    print(f"Weight range: {min(final_weights):.6f} to {max(final_weights):.6f}")
    
    # Print top 10 features by absolute weight
    abs_weights = [abs(w) for w in final_weights]
    top_indices = np.argsort(abs_weights)[-10:][::-1]
    print(f"\nTop 10 features by absolute weight:")
    for i, idx in enumerate(top_indices):
        print(f"{i+1:2d}. {final_feature_names[idx][:60]}... : {final_weights[idx]:.6f}")

if __name__ == "__main__":
    model_file = "GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP_model.pkl"
    plot_trained_model_bars(model_file) 