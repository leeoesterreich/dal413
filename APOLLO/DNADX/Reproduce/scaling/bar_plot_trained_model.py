#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import joblib
import pickle
from matplotlib.patches import Rectangle, Patch
import argparse
import os

# hg19 chromosome lengths (adjust as needed, e.g., from UCSC or NCBI)
# p-arm lengths (approximate, used to define centromere start for q-arm)
# For acrocentric chromosomes (13, 14, 15, 21, 22), p-arm is very short, effectively 0 for this plotting.
# X and Y p-arm lengths are also included.
CYTOBAND_HG19 = {
    '1': 125000000, '2': 93300000, '3': 91000000, '4': 50400000, '5': 48400000,
    '6': 61000000, '7': 59900000, '8': 45600000, '9': 49000000, '10': 40200000,
    '11': 53700000, '12': 35800000, '13': 17900000, '14': 17600000, '15': 19000000,
    '16': 36600000, '17': 24000000, '18': 17200000, '19': 26500000, '20': 27500000,
    '21': 13200000, '22': 14700000, 'X': 60600000
    # 'Y': 12500000 # Y chromosome excluded
}

CHR_LENGTHS_HG19 = {
    '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260,
    '6': 171115067, '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747,
    '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
    '16': 90354753, '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
    '21': 48129895, '22': 51304566, 'X': 155270560
    # 'Y': 59373566 # Y chromosome excluded
}
ORDERED_CHROMS = [str(i) for i in range(1, 23)] + ['X'] # Y chromosome excluded

def parse_genomic_coordinates(feature_name_str, weight_val):
    """
    Parses genomic coordinates from a feature name string.
    Returns a dictionary with chr, start, end, absolute_start, absolute_end, and weight.
    Handles segment features (e.g., chr1:1000-2000) and whole arm features (e.g., chr1p_amp, chr1q_del).
    """
    feature_name_str_lower = feature_name_str.lower() # For case-insensitive checks on keywords

    # Whole arm features (e.g., "chr1p_amp", "chr1q", "chr1_whole_arm_gain")
    # Regex tries to capture chrom, arm (p/q), and if _whole_arm suffix is present
    # Reverting to a previously working style for whole-arm, with minor checks
    arm_match = re.match(r'chr(\d+|X|Y)([pq])?(_whole_arm)?([._]?(amp|del|gain|loss))?$', feature_name_str, re.IGNORECASE)
    if arm_match:
        chrom = arm_match.group(1)
        if chrom.upper() == 'Y': return None # Exclude Y chromosome
        arm = arm_match.group(2) # p or q, might be None
        has_whole_arm_suffix = bool(arm_match.group(3)) # True if _whole_arm is present

        if chrom not in CHR_LENGTHS_HG19:
            return None

        p_len = CYTOBAND_HG19.get(chrom, 0)
        total_len = CHR_LENGTHS_HG19[chrom]
        
        is_w_arm = False
        if has_whole_arm_suffix:
            is_w_arm = True
        elif arm and not (':' in feature_name_str and '-' in feature_name_str): # arm ('p' or 'q') present, no explicit range
            is_w_arm = True
        elif not arm and not (':' in feature_name_str and '-' in feature_name_str): # No arm, no explicit range -> whole chromosome
            is_w_arm = True


        if arm == 'p':
            start = 1
            end = p_len
        elif arm == 'q':
            start = p_len + 1
            end = total_len
        elif arm is None: # Whole chromosome if no arm specified or implied
            start = 1
            end = total_len
            is_w_arm = True # Explicitly whole chromosome
        else: 
            return None
        
        # Final check for is_w_arm for cases like "chr1p" without suffix, should be whole arm
        if arm and not has_whole_arm_suffix and not (':' in feature_name_str and '-' in feature_name_str) :
             is_w_arm = True


        return {
            'chr': chrom, 'arm': arm if arm else 'pq', 
            'start_genomic': start, 'end_genomic': end,
            'weight': weight_val, 'feature_name': feature_name_str,
            'is_whole_arm': is_w_arm
        }

    # Segment features (e.g., "chr1:10000-20000")
    # Trying a new regex to capture chr:start-end, then anything, then amp/del.
    # Removed $ from end and will pass re.IGNORECASE to re.match()
    segment_match = re.match(r'^(chr(\d+|X|Y):(\d+)-(\d+))\S*?(amp|del|gain|loss)', feature_name_str, re.IGNORECASE)
    # Group 1: full chr:start-end part, Group 2: chr, Group 3: start, Group 4: end, Group 5: amp/del
    if segment_match:
        # captured_segment_base = segment_match.group(1) # e.g. chr1:10000-20000
        chrom = segment_match.group(2)
        if chrom.upper() == 'Y': return None # Exclude Y chromosome
        start = int(segment_match.group(3))
        end = int(segment_match.group(4))
        # alteration_type = segment_match.group(5) # amp or del

        if chrom not in CHR_LENGTHS_HG19:
            return None
        
        p_len = CYTOBAND_HG19.get(chrom, 0)
        mid_point = (start + end) / 2
        arm_val = 'p' if mid_point <= p_len else 'q'
        if chrom in ['13', '14', '15', '21', '22'] and arm_val == 'p':
             arm_val = 'q'

        return {
            'chr': chrom, 'arm': arm_val,
            'start_genomic': start, 'end_genomic': end,
            'weight': weight_val, 'feature_name': feature_name_str,
            'is_whole_arm': False # Segments are not whole arm
        }
    
    # Legacy or other formats (e.g., "1.p.wholearm") - Attempt to parse
    legacy_arm_match = re.match(r'(\\d+|X|Y)\\.([pq])\\.wholearm', feature_name_str, re.IGNORECASE)
    if legacy_arm_match:
        chrom = legacy_arm_match.group(1)
        if chrom.upper() == 'Y': return None # Exclude Y chromosome
        arm = legacy_arm_match.group(2)
        if chrom not in CHR_LENGTHS_HG19:
            return None
        p_len = CYTOBAND_HG19.get(chrom, 0)
        total_len = CHR_LENGTHS_HG19[chrom]
        if arm == 'p':
            start, end = 1, p_len
        else: # q arm
            start, end = p_len + 1, total_len
        return {
            'chr': chrom, 'arm': arm,
            'start_genomic': start, 'end_genomic': end,
            'weight': weight_val, 'feature_name': feature_name_str,
            'is_whole_arm': True # These are explicitly wholearm
        }

    # New specific format: e.g., "3.q   wholearm.chr3:91700000-199501827" or "X.qwholearm.chrX:..."
    # Made the space between arm and "wholearm" optional (\s*)
    specific_format_match = re.match(r'^(\d+|X|Y)\.([pq])\s*wholearm\.chr\1:(\d+)-(\d+)', feature_name_str, re.IGNORECASE)
    if specific_format_match:
        chrom = specific_format_match.group(1)
        if chrom.upper() == 'Y': return None # Exclude Y chromosome
        arm = specific_format_match.group(2)
        # start_coord = int(specific_format_match.group(3)) # These coords can be used for validation if needed
        # end_coord = int(specific_format_match.group(4))

        if chrom not in CHR_LENGTHS_HG19:
            return None
        
        p_len = CYTOBAND_HG19.get(chrom, 0)
        total_len = CHR_LENGTHS_HG19[chrom]
        start_genomic, end_genomic = (1, p_len) if arm == 'p' else (p_len + 1, total_len)

        return {
            'chr': chrom, 'arm': arm,
            'start_genomic': start_genomic, 'end_genomic': end_genomic,
            'weight': weight_val, 'feature_name': feature_name_str,
            'is_whole_arm': True
        }

    return None


def plot_trained_model_bars(model_file, excel_data=None):
    """Plot feature landscape using bar plot, for trained models or Excel data."""
    
    # Determine plot title based on model_file content
    if "Basal_signaling" in model_file:
        plot_title = "Feature landscape of basal signaling"
    elif "Estrogen_signaling" in model_file:
        plot_title = "Feature landscape of estrogen signaling"
    elif "RB_LOH" in model_file:
        plot_title = "Feature landscape of RB_LOH"
    else:
        plot_title = f'Feature Landscape from {model_file}' # Default if no keywords match

    if excel_data is not None:
        print(f"Processing data from Excel file: {model_file}") # model_file is used as a label here
        feature_names = excel_data['Genomic Region'].tolist()
        coefficients = excel_data['Weight'].tolist()
        # Title for Excel data might need separate logic or use the filename as before
        if "Basal_signaling" in model_file or "Estrogen_signaling" in model_file or "RB_LOH" in model_file:
            pass # Title already set
        else:
            plot_title = f'Feature Landscape from {model_file}' # Override for excel if keywords not in label
    else:
        print(f"Loading model file: {model_file}")
        # Title for model files is already set above based on path
        try:
            try:
                model_data = joblib.load(model_file)
                print("Loaded with joblib")
            except:
                with open(model_file, 'rb') as f:
                    try:
                        model_data = pickle.load(f)
                        print("Loaded with pickle")
                    except:
                        f.seek(0)
                        import pickle5
                        model_data = pickle5.load(f)
                        print("Loaded with pickle5")
            
            if isinstance(model_data, dict):
                # Use the keys as saved by model_training.py
                coefficients = model_data.get('coefficients_on_original_scale')
                feature_names = model_data.get('selected_feature_names', [])
                # The raw model object might be needed if other attributes were to be checked,
                # but for coefficients and features, use the pre-processed ones.
                # model = model_data.get('model_object') 
                
                if coefficients is None:
                    print("Error: 'coefficients_on_original_scale' not found in model file.")
                    return
                if not feature_names: # check if list is empty
                    print("Error: 'selected_feature_names' not found or empty in model file.")
                    return

            else: # Assuming model_data is the model itself (legacy, not expected from model_training.py)
                print("Warning: Model file is not a dictionary. Attempting to load as a raw model object.")
                model = model_data
                if hasattr(model, 'coef_'):
                    coefficients = model.coef_
                    if len(coefficients.shape) > 1: # For models like LogisticRegressionCV
                        if coefficients.shape[0] == 1:
                            coefficients = coefficients[0]
                        else:
                            print(f"Warning: Raw model has {coefficients.shape[0]} sets of coefficients. Using the first set.")
                            coefficients = coefficients[0]
                else:
                    print("Error: Raw model object does not have coef_ attribute.")
                    return

                if hasattr(model, 'feature_names_in_'):
                    feature_names = model.feature_names_in_
                elif hasattr(model, 'features_'): # Some older models might use this
                    feature_names = model.features_
                else: # Fallback if no feature names found on raw model
                    print(f"Warning: No feature names found on raw model object. Using dummy names.")
                    if coefficients is not None:
                         feature_names = [f"feature_{i}" for i in range(len(coefficients))]
                    else:
                         feature_names = []
            
            if not list(feature_names): # If feature_names is still empty or not list-like
                print(f"Critical Error: No feature names could be determined. Aborting.")
                return

            if len(feature_names) != len(coefficients):
                print(f"Critical Error: Feature names ({len(feature_names)}) and coefficients ({len(coefficients)}) length mismatch after loading. Aborting.")
                print(f"Feature names sample: {list(feature_names)[:5]}")
                return
        
        except Exception as e:
            print(f"Error loading model or extracting coefficients: {e}")
            import traceback
            traceback.print_exc()
            return

    non_zero_mask = np.array(coefficients) != 0
    non_zero_features_names = [str(name) for name, mask_val in zip(feature_names, non_zero_mask) if mask_val]
    non_zero_weights = np.array(coefficients)[non_zero_mask]

    print(f"Total features in model/file: {len(coefficients)}")
    print(f"Non-zero features: {len(non_zero_weights)}")

    if len(non_zero_weights) == 0:
        print("No non-zero weights found!")
        return

    parsed_features = []
    unparsed_feature_examples = [] # New list to store examples
    for name, weight in zip(non_zero_features_names, non_zero_weights):
        parsed = parse_genomic_coordinates(name, weight)
        if parsed:
            parsed_features.append(parsed)
        elif len(unparsed_feature_examples) < 10: # Store up to 10 examples
            unparsed_feature_examples.append(name)
    
    print(f"Successfully parsed {len(parsed_features)} features with coordinates.")
    if not parsed_features and non_zero_features_names: # Check if there were features to parse
        print("No features could be parsed into genomic coordinates. Cannot generate plot.")
        print("Sample feature names that were not parsed (up to 10 shown):")
        for i, feature_name_example in enumerate(unparsed_feature_examples):
            print(f"  {i+1}: {feature_name_example}")
        return
    elif not parsed_features:
        print("No non-zero features to parse or plot.")
        return
    elif unparsed_feature_examples: # If some were parsed but others failed
        print(f"{len(unparsed_feature_examples)} other feature(s) could not be parsed. Examples (up to 10 shown): ")
        for i, feature_name_example in enumerate(unparsed_feature_examples):
            print(f"  {i+1}: {feature_name_example}")

    # Calculate cumulative lengths for x-axis
    chromosome_cumulative_lengths = {}
    current_offset = 0
    for chrom_name in ORDERED_CHROMS:
        chromosome_cumulative_lengths[chrom_name] = current_offset
        current_offset += CHR_LENGTHS_HG19[chrom_name]
    total_genome_length = current_offset

    # Assign absolute genomic start/end and midpoints for sorting and plotting
    for feature in parsed_features:
        chrom_offset = chromosome_cumulative_lengths[feature['chr']]
        feature['abs_start_pos'] = chrom_offset + feature['start_genomic']
        feature['abs_end_pos'] = chrom_offset + feature['end_genomic']
        feature['plot_pos_on_genome'] = chrom_offset + (feature['start_genomic'] + feature['end_genomic']) / 2


    # Sort features by their absolute genomic position
    parsed_features.sort(key=lambda x: (x['plot_pos_on_genome']))
    
    fig, ax = plt.subplots(figsize=(30, 12)) # Increased figure size, added height for legend

    # Separate features for plotting
    whole_arm_features = [f for f in parsed_features if f.get('is_whole_arm', False)]
    segment_features = [f for f in parsed_features if not f.get('is_whole_arm', False)]

    # Plot whole arms as rectangles
    for feature in whole_arm_features:
        x_rect = feature['abs_start_pos']
        width_rect = feature['abs_end_pos'] - feature['abs_start_pos']
        if width_rect <= 0: # Safety check
            # print(f"Warning: Skipping whole arm feature with non-positive width: {feature['feature_name']}")
            continue
        y_rect = 0
        height_rect = feature['weight']
        color_rect = 'lightcoral' if height_rect > 0 else 'lightskyblue'
        rect = Rectangle((x_rect, y_rect), width_rect, height_rect, color=color_rect, alpha=0.5, ec='none')
        ax.add_patch(rect)

    # Plot segments as thin bars
    if segment_features:
        # New: Plot segments as Rectangles with width proportional to genomic length
        for feature in segment_features:
            x_rect_seg = feature['abs_start_pos']
            width_rect_seg = feature['abs_end_pos'] - feature['abs_start_pos']
            
            # Ensure minimum width for very small segments to be visible, can be a small fraction of total or a fixed small value.
            # This is tricky: if too small, they vanish. If too large, they misrepresent.
            # Let's try a very small minimum genomic unit, e.g. 1/10000th of genome width, or ensure it's at least 1 pixel on a 3000px wide image.
            # min_pixel_width_on_genome = total_genome_length / 3000 # Roughly 1 pixel on typical output
            # if width_rect_seg < min_pixel_width_on_genome:
            #     width_rect_seg = min_pixel_width_on_genome
            # For now, let actual genomic length dictate width. User wants length proportional.
            # If width is zero or negative (should not happen with valid start < end), skip or give minimal default.
            if width_rect_seg <= 0:
                # print(f"Warning: Segment feature with non-positive width: {feature['feature_name']}. Assigning minimal width.")
                # This might happen if start_genomic == end_genomic. Assign a tiny width.
                width_rect_seg = total_genome_length / 10000 # A very small default width for point-like features

            y_rect_seg = 0
            height_rect_seg = feature['weight']
            color_rect_seg = 'maroon' if height_rect_seg > 0 else 'navy'
            
            segment_rect = Rectangle((x_rect_seg, y_rect_seg), width_rect_seg, height_rect_seg, 
                                     color=color_rect_seg, alpha=0.9, ec='none')
            ax.add_patch(segment_rect)
    
    ax.set_ylabel('Weight', fontsize=14)
    ax.set_title(plot_title, fontsize=16, fontweight='bold')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, axis='y')
    ax.tick_params(axis='y', labelsize=14) # Set Y-tick font size to 14

    # Set x-axis limits and ticks
    ax.set_xlim(0, total_genome_length)
    # Set fixed y-limits
    ax.set_ylim(-0.2, 0.2)
    
    chr_tick_positions = []
    chr_tick_labels = []
    last_chrom_end = 0
    for chrom_name in ORDERED_CHROMS:
        chrom_start_abs = chromosome_cumulative_lengths[chrom_name]
        chrom_end_abs = chrom_start_abs + CHR_LENGTHS_HG19[chrom_name]
        
        # Vertical line at the start of each chromosome (except the first one)
        if chrom_start_abs > 0:
            ax.axvline(x=chrom_start_abs, color='grey', linestyle='--', linewidth=0.8, alpha=0.7)
        
        # Chromosome label centered in its region
        label_pos = chrom_start_abs + (CHR_LENGTHS_HG19[chrom_name] / 2)
        chr_tick_positions.append(label_pos)
        chr_tick_labels.append(chrom_name) # Changed from f'Chr{chrom_name}'
        
        # Draw p and q arm boundary if not acrocentric
        if chrom_name not in ['13', '14', '15', '21', '22', 'Y']: # Y is also effectively acrocentric for this
            p_len = CYTOBAND_HG19.get(chrom_name, 0)
            centromere_pos_abs = chrom_start_abs + p_len
            ax.axvline(x=centromere_pos_abs, color='lightgrey', linestyle=':', linewidth=0.6, alpha=0.5)
        last_chrom_end = chrom_end_abs

    ax.set_xticks(chr_tick_positions)
    ax.set_xticklabels(chr_tick_labels, fontsize=14, rotation=0, ha="center") # Changed X-tick fontsize to 14
    ax.set_xlabel('Chromosome', fontsize=14) # Changed label text


    legend_elements = [
        Patch(facecolor='maroon', alpha=0.9, label='Positive segment'),
        Patch(facecolor='navy', alpha=0.9, label='Negative segment'),
        Patch(facecolor='lightcoral', alpha=0.5, label='Positive whole arm'),
        Patch(facecolor='lightskyblue', alpha=0.5, label='Negative whole arm')
    ]
    # Position legend at the bottom center
    ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.18), ncol=4, fontsize=14) # Changed legend fontsize to 14

    plt.subplots_adjust(bottom=0.2) # Adjust bottom margin to make space for legend
    # plt.tight_layout() # tight_layout might conflict with subplots_adjust, call before or manage manually

    # Define output directory
    output_dir = "adjusted_output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # Generate output filename
    base_name = model_file.split('/')[-1]
    if excel_data is not None:
         # model_file is a label like "Excel_Data_Full"
        plot_file_name = base_name.replace(' ', '_') 
    else:
        plot_file_name = base_name.replace('.pkl', '').replace('.joblib', '')

    # Add a suffix for plots from the L1 search
    output_file = os.path.join(output_dir, f'{plot_file_name}_L1search_genomic_bars.png') 
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")

    # Print statistics
    final_plotted_weights = [f['weight'] for f in parsed_features]
    final_plotted_names = [f['feature_name'] for f in parsed_features]
    
    print(f"\nStatistics for plotted features:")
    print(f"Total features plotted: {len(final_plotted_weights)}")
    if final_plotted_weights:
        print(f"Positive weights: {sum(1 for w in final_plotted_weights if w > 0)}")
        print(f"Negative weights: {sum(1 for w in final_plotted_weights if w < 0)}")
        print(f"Weight range: {min(final_plotted_weights):.6f} to {max(final_plotted_weights):.6f}")

        # Print top 10 features by absolute weight
        abs_weights_tuples = sorted([(abs(w), name, w) for w, name in zip(final_plotted_weights, final_plotted_names)], reverse=True)
        print(f"\nTop 10 plotted features by absolute weight:")
        for i, (abs_w, name, w) in enumerate(abs_weights_tuples[:10]):
            print(f"{i+1:2d}. {name[:70]}... : {w:.6f}")
    else:
        print("No features were plotted.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot feature landscape from a trained model file or Excel data.")
    parser.add_argument("model_file", help="Path to the model file (.pkl or .joblib) or Excel file (.xlsx)")
    parser.add_argument("--output_csv", help="Optional path to save non-zero features and weights to a CSV file.") # New argument
    args = parser.parse_args()

    if args.model_file.endswith('.xlsx'):
        try:
            excel_df = pd.read_excel(args.model_file)
            # Assuming columns are 'Genomic Region' and 'Weight'
            if 'Genomic Region' not in excel_df.columns or 'Weight' not in excel_df.columns:
                print(f"Error: Excel file {args.model_file} must contain 'Genomic Region' and 'Weight' columns.")
            else:
                plot_trained_model_bars(args.model_file, excel_data=excel_df) # Pass filename as label
        except Exception as e:
            print(f"Error processing Excel file {args.model_file}: {e}")
            import traceback
            traceback.print_exc()
    else:
        # For model files, we can extract and save features if --output_csv is given
        if args.output_csv:
            print(f"Attempting to extract features to CSV: {args.output_csv}")
            # Minimal model loading logic to get features and weights for CSV
            # This duplicates some loading logic from plot_trained_model_bars but is cleaner for this specific task
            try:
                try:
                    model_data_for_csv = joblib.load(args.model_file)
                except:
                    with open(args.model_file, 'rb') as f_csv:
                        try:
                            model_data_for_csv = pickle.load(f_csv)
                        except:
                            f_csv.seek(0)
                            import pickle5
                            model_data_for_csv = pickle5.load(f_csv)
                
                csv_model, csv_feature_names, csv_coefficients = None, [], None

                if isinstance(model_data_for_csv, dict):
                    # For CSV extraction, also use the corrected keys
                    csv_coefficients = model_data_for_csv.get('coefficients_on_original_scale')
                    csv_feature_names = model_data_for_csv.get('selected_feature_names', [])
                    # csv_model = model_data_for_csv.get('model_object') # if other model props needed
                else: # Legacy raw model handling for CSV
                    csv_model = model_data_for_csv
                    if hasattr(csv_model, 'coef_'):
                        csv_coefficients = csv_model.coef_
                        if len(csv_coefficients.shape) > 1:
                            if csv_coefficients.shape[0] == 1:
                                csv_coefficients = csv_coefficients[0]
                            else:
                                print(f"Warning for CSV: Model has {csv_coefficients.shape[0]} sets of coefficients. Using the first set.")
                                csv_coefficients = csv_coefficients[0]
                    
                    if hasattr(csv_model, 'feature_names_in_'):
                        csv_feature_names = csv_model.feature_names_in_
                    elif hasattr(csv_model, 'features_'):
                        csv_feature_names = csv_model.features_
                    elif csv_coefficients is not None:
                         csv_feature_names = [f"feature_{i}" for i in range(len(csv_coefficients))]
                    else:
                        csv_feature_names = []

                if csv_coefficients is not None and list(csv_feature_names) and len(csv_feature_names) == len(csv_coefficients):
                    non_zero_mask_csv = np.array(csv_coefficients) != 0
                    nz_feat_names_csv = [str(name) for name, mask_val in zip(csv_feature_names, non_zero_mask_csv) if mask_val]
                    nz_weights_csv = np.array(csv_coefficients)[non_zero_mask_csv]

                    if nz_feat_names_csv:
                        features_df = pd.DataFrame({
                            'FeatureName': nz_feat_names_csv,
                            'Weight': nz_weights_csv
                        })
                        # Sort by absolute weight, descending, for consistency with plot top features
                        features_df['AbsWeight'] = features_df['Weight'].abs()
                        features_df = features_df.sort_values(by='AbsWeight', ascending=False).drop(columns=['AbsWeight'])
                        
                        output_csv_path = args.output_csv
                        # Ensure the output directory for the CSV exists if it's in a subdirectory
                        csv_output_dir = os.path.dirname(output_csv_path)
                        if csv_output_dir and not os.path.exists(csv_output_dir):
                            os.makedirs(csv_output_dir)
                            print(f"Created directory for CSV: {csv_output_dir}")
                        features_df.to_csv(output_csv_path, index=False)
                        print(f"Non-zero features saved to {output_csv_path}")
                    else:
                        print("No non-zero features found to save to CSV.")
                else:
                    print("Could not extract features or coefficients for CSV output, or lengths mismatch.")
            except Exception as e_csv:
                print(f"Error during CSV feature extraction: {e_csv}")
                import traceback
                traceback.print_exc()
        
        # Always attempt to plot after (optional) CSV generation
        plot_trained_model_bars(args.model_file)