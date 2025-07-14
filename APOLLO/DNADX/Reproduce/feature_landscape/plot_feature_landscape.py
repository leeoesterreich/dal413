import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import re
import joblib # Added joblib import

# Load the model
model_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/feature_landscape/GSEA_Median_GP17_Basal_signaling.r_0.958_SMID_BREAST_CANCER_BASAL_UP_scaleX_scaleY_l1_0.05_elasticnet_model_components.pkl"
print("Loading model from: {}".format(model_path))
try:
    model_data = joblib.load(model_path)
except Exception as e:
    print("Error loading PKL model: {}".format(e))
    exit()

# Extract feature names and coefficients
if isinstance(model_data, dict) and 'model' in model_data and 'feature_names' in model_data:
    print("Model loaded as a dictionary. Extracting 'model' and 'feature_names'.")
    model = model_data['model']
    feature_names = model_data['feature_names']
    if hasattr(model, 'coef_'):
        coefficients = model.coef_
    elif hasattr(model, 'feature_importances_'): # For models like RandomForest
        coefficients = model.feature_importances_
    else:
        print("Error: Couldn't extract coefficients from the model inside the dictionary.")
        exit()
elif hasattr(model_data, 'coef_') and hasattr(model_data, 'feature_names_in_'): # Direct scikit-learn model
    print("Model loaded as a scikit-learn model object.")
    model = model_data
    feature_names = model.feature_names_in_
    coefficients = model.coef_
else:
    print("Error: Model structure not recognized. Cannot extract features or coefficients.")
    print("Model type: {}".format(type(model_data)))
    if isinstance(model_data, dict):
        print("Dictionary keys: {}".format(model_data.keys()))
    exit()

if coefficients.ndim > 1: # Handle cases where coef_ might be 2D (e.g. multi-class)
    coefficients = coefficients[0]


print("Found {} feature names".format(len(feature_names)))
print("Found {} coefficients".format(len(coefficients)))
print("Found {} non-zero coefficients".format(np.sum(coefficients != 0)))

# Print first few features and their coefficients for inspection
print("\nFirst few features and their coefficients:")
for i in range(min(10, len(feature_names))):
    print("Feature: {:<50} Coefficient: {}".format(str(feature_names[i]), coefficients[i]))

# Read the segment annotation file
segment_file = "CNA_segments.utf8.gmt"
print("Reading segment annotations from: {}".format(segment_file))

segment_annotations = {}
with open(segment_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            segment_name = parts[0]
            # Parse genomic coordinates from segment name
            coord_match = re.search(r'chr(\d+|X|Y):(\d+)-(\d+)', segment_name)
            if coord_match:
                chrom = coord_match.group(1)
                start = int(coord_match.group(2))
                end = int(coord_match.group(3))
                segment_annotations[segment_name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'is_wholearm': 'wholearm' in segment_name.lower()
                }
            elif 'wholearm' in segment_name.lower():
                # Handle whole-arm features
                arm_match = re.match(r'(\d+|X|Y)\.(p|q)\s+wholearm', segment_name)
                if arm_match:
                    chrom = arm_match.group(1)
                    arm = arm_match.group(2)
                    segment_annotations[segment_name] = {
                        'chrom': chrom,
                        'arm': arm,
                        'is_wholearm': True
                    }

# Create mapping of feature names to coefficients and segment info
features_data = []
for i, (feature, coef) in enumerate(zip(feature_names, coefficients)):
    if pd.notnull(coef) and coef != 0 and isinstance(feature, str):
        # Try to find matching segment annotation
        segment_info = None
        for segment_name, info in segment_annotations.items():
            if segment_name == feature:  # Exact match
                segment_info = info
                break
        
        if segment_info:
            features_data.append({
                'name': feature,
                'coef': float(coef),
                'chrom': segment_info['chrom'],
                'start': segment_info.get('start', 0),
                'end': segment_info.get('end', 0),
                'arm': segment_info.get('arm', ''),
                'is_wholearm': segment_info['is_wholearm']
            })

print("\nProcessed {} non-zero features".format(len(features_data)))

# Print feature distribution by chromosome
chrom_counts = {}
for feature in features_data:
    chrom = feature['chrom']
    if chrom not in chrom_counts:
        chrom_counts[chrom] = {'total': 0, 'p_arm': 0, 'q_arm': 0}
    chrom_counts[chrom]['total'] += 1
    if feature['arm'] == 'p':
        chrom_counts[chrom]['p_arm'] += 1
    elif feature['arm'] == 'q':
        chrom_counts[chrom]['q_arm'] += 1

print("\nFeature distribution by chromosome:")
for chrom in sorted(chrom_counts.keys()):
    counts = chrom_counts[chrom]
    print("Chromosome {}: Total={}, p-arm={}, q-arm={}".format(
        chrom, counts['total'], counts['p_arm'], counts['q_arm']))

# Print some example parsed features
print("\nExample parsed features:")
for i, feature in enumerate(features_data[:5]):
    print("Feature {}: {}".format(i+1, feature))

# Define chromosome order and lengths
chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y']

# Group features by chromosome
chrom_features = {}
for chrom in chrom_order:
    chrom_features[chrom] = {
        'p_arm': {'wholearm': [], 'segments': []},
        'q_arm': {'wholearm': [], 'segments': []}
    }

# Organize features by chromosome and arm
for feature in features_data:
    chrom = feature['chrom']
    
    if chrom not in chrom_features:
        continue
    
    # Determine arm bucket ('p_arm' or 'q_arm')
    if feature['is_wholearm']:
        arm_bucket = 'p_arm' if feature.get('arm') == 'p' else 'q_arm'
    else:
        # For non-wholearm features, use genomic coordinates to determine arm
        # Assuming centromere position is roughly at 1/3 of chromosome length
        centromere_pos = feature['start'] + (feature['end'] - feature['start']) / 3
        arm_bucket = 'p_arm' if feature['start'] < centromere_pos else 'q_arm'
    
    # Determine feature type bucket ('wholearm' or 'segments')
    type_bucket = 'wholearm' if feature['is_wholearm'] else 'segments'
    
    # Add feature to the appropriate bucket
    chrom_features[chrom][arm_bucket][type_bucket].append(feature)

# Create the plot
plt.figure(figsize=(15, 5))

# Set up positions for each chromosome on the x-axis
chrom_positions = {chrom: i for i, chrom in enumerate(chrom_order)}

# Set up color scheme
pos_color = 'red'
neg_color = 'blue'
pos_bg_color = 'pink'
neg_bg_color = 'lightblue'

# Calculate y-axis limits
all_coefs = [f['coef'] for f in features_data]
max_coef = max(all_coefs) if all_coefs else 0.4
min_coef = min(all_coefs) if all_coefs else -0.4
y_range = max(abs(max_coef), abs(min_coef))
y_max = 0.15  # Fixed y-axis maximum
y_min = -0.15  # Fixed y-axis minimum

# Plot features chromosome by chromosome
for chrom in chrom_order:
    x_pos = chrom_positions[chrom]
    
    # Check for whole-arm features in p and q arms
    p_wholearm = chrom_features[chrom]['p_arm']['wholearm']
    q_wholearm = chrom_features[chrom]['q_arm']['wholearm']
    p_segments = chrom_features[chrom]['p_arm']['segments']
    q_segments = chrom_features[chrom]['q_arm']['segments']
    
    # Highlight chromosome background for whole-arm features
    p_has_pos = any(f['coef'] > 0 for f in p_wholearm)
    p_has_neg = any(f['coef'] < 0 for f in p_wholearm)
    q_has_pos = any(f['coef'] > 0 for f in q_wholearm)
    q_has_neg = any(f['coef'] < 0 for f in q_wholearm)
    
    # Background for p-arm (left half of chromosome column)
    if p_has_pos:
        plt.axvspan(x_pos - 0.45, x_pos, ymin=0, ymax=1, color=pos_bg_color, alpha=0.3, zorder=-1)
    elif p_has_neg:
        plt.axvspan(x_pos - 0.45, x_pos, ymin=0, ymax=1, color=neg_bg_color, alpha=0.3, zorder=-1)
    
    # Background for q-arm (right half of chromosome column)
    if q_has_pos:
        plt.axvspan(x_pos, x_pos + 0.45, ymin=0, ymax=1, color=pos_bg_color, alpha=0.3, zorder=-1)
    elif q_has_neg:
        plt.axvspan(x_pos, x_pos + 0.45, ymin=0, ymax=1, color=neg_bg_color, alpha=0.3, zorder=-1)
    
    # Plot p-arm wholearm features
    for feature in p_wholearm:
        color = pos_color if feature['coef'] > 0 else neg_color
        plt.vlines(x_pos - 0.25, 0, feature['coef'], color=color, linewidth=2)
    
    # Plot q-arm wholearm features
    for feature in q_wholearm:
        color = pos_color if feature['coef'] > 0 else neg_color
        plt.vlines(x_pos + 0.25, 0, feature['coef'], color=color, linewidth=2)
    
    # Plot p-arm segment features
    if p_segments:
        # Calculate positions for p-arm segments based on genomic coordinates
        for feature in p_segments:
            rel_pos = (feature['start'] + feature['end']) / 2 / 24776  # Scale to [0,1]
            x_offset = -0.4 + rel_pos * 0.3  # Map to [-0.4, -0.1]
            color = pos_color if feature['coef'] > 0 else neg_color
            plt.vlines(x_pos + x_offset, 0, feature['coef'], color=color, linewidth=1.5)
    
    # Plot q-arm segment features
    if q_segments:
        # Calculate positions for q-arm segments based on genomic coordinates
        for feature in q_segments:
            rel_pos = (feature['start'] + feature['end']) / 2 / 24776  # Scale to [0,1]
            x_offset = 0.1 + rel_pos * 0.3  # Map to [0.1, 0.4]
            color = pos_color if feature['coef'] > 0 else neg_color
            plt.vlines(x_pos + x_offset, 0, feature['coef'], color=color, linewidth=1.5)

# Add vertical grid lines between chromosomes
for pos in range(len(chrom_order) + 1):
    plt.axvline(x=pos - 0.5, color='gray', linestyle='-', linewidth=0.5, alpha=0.7)

# Add horizontal line at y=0
plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)

# Set x-ticks and labels
plt.xticks(range(len(chrom_order)), chrom_order)
plt.xlim(-0.5, len(chrom_order) - 0.5)

# Set y-limits
plt.ylim(y_min, y_max)

# Add labels and title
plt.ylabel('Weight', fontsize=12)
plt.xlabel('Chromosome', fontsize=12)
plt.title('Basal signaling (New Model)', fontsize=14)

# Add grid and finalize plot
plt.grid(True, axis='x', alpha=0.3)
plt.grid(False, axis='y')

# Save the finalized plot
plt.tight_layout()
plt.savefig('Basal_signaling_feature_landscape_new_model.png', dpi=300)
print("\nFeature landscape plotted and saved as Basal_signaling_feature_landscape_new_model.png") 