import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set the backend to Agg for non-interactive environments
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Wedge

# Read the CSV file
df = pd.read_csv('IO_summary.csv')

# Reorder columns to match desired sequence
column_order = ['Trial', 'Biomarker', 'Molecular subtype', 'Stage', 'IO Drug', 'Primary Endpoint', 'Stratification']
df = df[column_order]

# Calculate initial rotation to make PD-L1 horizontal
pd_l1_data = df[df['Biomarker'] == 'PD-L1']
if not pd_l1_data.empty:
    total_count = len(df)
    # Find the starting position of PD-L1
    pd_l1_start = len(df[df['Biomarker'].isin(df['Biomarker'].unique()[df['Biomarker'].unique() < 'PD-L1'])])
    pd_l1_angle = (pd_l1_start / total_count) * 360
    # Calculate rotation needed to make PD-L1 start at 0 degrees
    initial_rotation = -pd_l1_angle
else:
    initial_rotation = 0

# Set up the figure
fig, ax = plt.subplots(figsize=(15, 12))  # Reduced from (20, 15) to make plot relatively larger
ax.set_aspect('equal')

# Base colors for each level - using a more harmonious color palette
base_colors = {
    'Biomarker': '#2E4057',      # Deep blue-gray
    'Molecular subtype': '#16DB93',  # Mint green
    'Stage': '#EFEA5A',          # Soft yellow
    'IO Drug': '#048BA8',        # Teal
    'Primary Endpoint': '#F29E4C',  # Warm orange
    'Stratification': '#8B8B8B'   # Neutral gray
}

# Function to create color variations
def create_color_variations(base_color, n):
    # Convert hex to RGB
    rgb = matplotlib.colors.hex2color(base_color)
    # Create n variations by adjusting lightness and saturation
    variations = []
    for i in range(n):
        # Create variations that maintain color identity while providing contrast
        factor = 0.4 + (0.6 * i / (n-1 if n > 1 else 1))  # Range from 0.4 to 1.0
        # Adjust both lightness and saturation
        variation = tuple(min(1.0, c * factor + (1 - factor) * 0.8) for c in rgb)
        variations.append(variation)
    return variations

# Function to handle values without splitting
def process_values(series):
    # Create a dictionary to store values and their counts
    value_dict = {}
    
    for value in series:
        if pd.isna(value):
            value = 'NA'
        else:
            value = str(value).strip()
        value_dict[value] = value_dict.get(value, 0) + 1
    
    return pd.Series(value_dict)

# Function to process molecular subtypes without splitting
def process_molecular_subtypes(df):
    # Create a dictionary to store subtypes
    subtype_dict = {}
    
    for value in df['Molecular subtype']:
        if pd.isna(value):
            value = 'NA'
        else:
            value = str(value).strip()
        subtype_dict[value] = subtype_dict.get(value, 0) + 1
    
    return pd.Series(subtype_dict)

# Create legend elements
legend_elements = []

# First, get the biomarker groups and their counts
biomarker_groups = df.groupby('Biomarker')
biomarker_counts = df['Biomarker'].value_counts()
total_biomarkers = len(df)

# Plot each level
n_levels = len(df.columns) - 1  # Subtract 1 to exclude Trial
max_radius = 0  # Track the maximum radius for dividing lines

# Calculate biomarker circle dimensions first - enlarged
trial_inner = 1.4  # New innermost circle for trials
trial_outer = trial_inner + 0.3  # Thin ring for trials
biomarker_inner = trial_outer + 0.3  # Adjusted to make room for trial circle
biomarker_outer = biomarker_inner + 1.0  # Keep same width as before

# Function to calculate text rotation for consistent downward facing text
def get_text_rotation(angle):
    # Normalize angle to 0-360 range
    angle = (angle + initial_rotation) % 360  # Add initial rotation
    # Base rotation is 90 degrees counter-clockwise from radius
    rotation = angle + 90
    # Adjust text to always face downward
    if rotation > 90 and rotation < 270:
        rotation = rotation - 180
    return rotation

# Store biomarker angles for alignment
biomarker_angle_map = {}
current_angle = 0

# First pass to store biomarker angles
for biomarker, group in biomarker_groups:
    biomarker_count = len(group)
    # Add spacing between biomarker segments by reducing the angle
    spacing_angle = 3  # degrees of spacing between biomarker segments
    biomarker_angle = (360 * (biomarker_count / total_biomarkers)) - spacing_angle
    biomarker_angle_map[biomarker] = {
        'start': current_angle + spacing_angle/2,  # Add half spacing at start
        'end': current_angle + biomarker_angle + spacing_angle/2,  # Add half spacing at end
        'trials': group['Trial'].unique()
    }
    current_angle += biomarker_angle + spacing_angle  # Add full spacing

# Create color maps for each column to ensure consistent colors for same values
color_maps = {}
for col in df.columns:
    if col not in ['Trial', 'Biomarker', 'Molecular subtype']:
        # Convert all values to strings first
        unique_values = [str(val) if not pd.isna(val) else 'NA' for val in df[col].unique()]
        unique_values = sorted(unique_values)
        colors = create_color_variations(base_colors[col], len(unique_values))
        color_maps[col] = dict(zip(unique_values, colors))

for i, col in enumerate(df.columns):
    if col == 'Trial':  # Skip the Trial column
        continue
        
    # Calculate the radius for this level
    if col == 'Trial':  # Add trial circle as innermost circle
        inner_radius = trial_inner
        outer_radius = trial_outer
    elif col == 'Biomarker':
        inner_radius = biomarker_inner
        outer_radius = biomarker_outer
    else:
        # More compact spacing for outer circles, starting after biomarker circle
        inner_radius = biomarker_outer + 0.2 + (i-2) * 0.3  # Adjusted i-2 since we added Trial
        outer_radius = inner_radius + 0.2  # Thinner rings for outer circles
    
    max_radius = max(max_radius, outer_radius)
    
    # Add to legend
    legend_elements.append(Patch(facecolor='white', edgecolor='none', label=f'\n{col}:'))
    
    # Plot segments
    current_angle = 0
    category_counts = {}  # Track total counts for each category
    
    if col == 'Trial':  # Special handling for trial circle
        # Draw trial segments aligned with biomarker segments
        for biomarker, info in biomarker_angle_map.items():
            biomarker_start = info['start']
            biomarker_angle = info['end'] - info['start']
            trials = info['trials']
            
            # Calculate angle per trial
            trial_angle = biomarker_angle / len(trials)
            
            # Draw each trial segment
            for idx, trial in enumerate(trials):
                start_angle = biomarker_start + (idx * trial_angle)
                end_angle = start_angle + trial_angle
                
                wedge = Wedge(
                    center=(0, 0),
                    r=outer_radius,
                    width=outer_radius-inner_radius,
                    theta1=start_angle + initial_rotation,
                    theta2=end_angle + initial_rotation,
                    facecolor='#E6E6E6',  # Light gray for all trial segments
                    edgecolor='white',
                    linewidth=1
                )
                ax.add_patch(wedge)
                
                # Add trial name if segment is large enough
                if trial_angle > 10:
                    text_angle = start_angle + trial_angle/2
                    text_radius = (inner_radius + outer_radius) / 2
                    text_x = text_radius * np.cos(np.radians(text_angle + initial_rotation))
                    text_y = text_radius * np.sin(np.radians(text_angle + initial_rotation))
                    
                    rotation = get_text_rotation(text_angle)
                    plt.text(text_x, text_y, str(trial),
                            rotation=rotation,
                            ha='center', va='center',
                            fontsize=6,  # Smaller font for trial names
                            weight='bold',
                            color='black')
                
                category_counts[trial] = 1
    
    elif col in ['Biomarker', 'Molecular subtype']:
        # Keep original layout for biomarker and molecular subtype
        for biomarker, group in biomarker_groups:
            biomarker_count = len(group)
            spacing_angle = 3  # degrees of spacing between biomarker segments
            biomarker_angle = (360 * (biomarker_count / total_biomarkers)) - spacing_angle
            
            value_counts = process_values(group[col])
            colors = create_color_variations(base_colors[col], len(value_counts))
            
            group_start_angle = current_angle + spacing_angle/2  # Add half spacing at start
            for (value, count), color in zip(value_counts.items(), colors):
                category_counts[str(value)] = category_counts.get(str(value), 0) + count
                angle = biomarker_angle * (count / biomarker_count)
                
                wedge = Wedge(
                    center=(0, 0),
                    r=outer_radius,
                    width=outer_radius-inner_radius,
                    theta1=group_start_angle + initial_rotation,
                    theta2=group_start_angle + angle + initial_rotation,
                    facecolor=color,
                    edgecolor='white',
                    linewidth=1
                )
                ax.add_patch(wedge)
                
                if angle > 10:
                    text_angle = group_start_angle + angle/2
                    text_radius = (inner_radius + outer_radius) / 2
                    text_x = text_radius * np.cos(np.radians(text_angle + initial_rotation))
                    text_y = text_radius * np.sin(np.radians(text_angle + initial_rotation))
                    
                    rotation = get_text_rotation(text_angle)
                    text_color = 'white' if col == 'Biomarker' else 'black'
                    plt.text(text_x, text_y, str(value),
                            rotation=rotation,
                            ha='center', va='center',
                            fontsize=10 if col == 'Biomarker' else 7,  # Larger font for biomarker
                            weight='bold',  # Bold font for all text
                            color=text_color)
                
                group_start_angle += angle
            
            if col == 'Biomarker':
                rad_angle = np.radians(current_angle + biomarker_angle + spacing_angle + initial_rotation)
                x = max_radius * np.cos(rad_angle)
                y = max_radius * np.sin(rad_angle)
                # Remove the line drawing
                # plt.plot([0, x], [0, y], 'k-', linewidth=0.5, alpha=0.2)
            
            current_angle += biomarker_angle + spacing_angle  # Add full spacing
    else:
        # For outer circles, align with trials within each biomarker segment
        for biomarker, info in biomarker_angle_map.items():
            biomarker_start = info['start']
            biomarker_angle = info['end'] - info['start']
            trials = info['trials']
            
            # Get all values for this biomarker's trials in the correct order
            trial_values = []
            trial_data = []
            for trial in trials:
                trial_row = df[df['Trial'] == trial].iloc[0]
                trial_data.append(trial_row)
                value = trial_row[col]
                if pd.isna(value):
                    value = 'NA'
                else:
                    value = str(value).strip()
                trial_values.append(value)
            
            # Calculate angle per trial
            trial_angle = biomarker_angle / len(trials)
            
            # Draw segments
            trial_idx = 0
            while trial_idx < len(trials):
                value = trial_values[trial_idx]
                
                # Find how many consecutive segments have the same value
                same_value_count = 1
                for next_idx in range(trial_idx + 1, len(trial_values)):
                    if trial_values[next_idx] == value:
                        same_value_count += 1
                    else:
                        break
                
                # Calculate angles for the merged segment
                start_angle = biomarker_start + (trial_idx * trial_angle)
                end_angle = start_angle + (trial_angle * same_value_count)
                
                # Update category counts
                category_counts[value] = category_counts.get(value, 0) + 1
                
                # Get color from color map
                color = color_maps[col][value]
                
                # Draw a single wedge for all consecutive same values
                wedge = Wedge(
                    center=(0, 0),
                    r=outer_radius,
                    width=outer_radius-inner_radius,
                    theta1=start_angle + initial_rotation,
                    theta2=end_angle + initial_rotation,
                    facecolor=color,
                    edgecolor='white',
                    linewidth=1
                )
                ax.add_patch(wedge)
                
                # Add text if the segment is large enough
                total_angle = trial_angle * same_value_count
                if total_angle > 10:
                    text_angle = start_angle + total_angle/2
                    text_radius = (inner_radius + outer_radius) / 2
                    text_x = text_radius * np.cos(np.radians(text_angle + initial_rotation))
                    text_y = text_radius * np.sin(np.radians(text_angle + initial_rotation))
                    
                    rotation = get_text_rotation(text_angle)
                    plt.text(text_x, text_y, str(value),
                            rotation=rotation,
                            ha='center', va='center',
                            fontsize=7,
                            weight='bold',
                            color='black')
                
                # Skip to the next different value
                trial_idx += same_value_count
    
    # Add aggregated category counts to legend
    if col in ['Biomarker', 'Molecular subtype']:
        colors = create_color_variations(base_colors[col], len(category_counts))
        for (value, count), color in zip(sorted(category_counts.items()), colors):
            legend_elements.append(Patch(facecolor=color, edgecolor='white',
                                       label=f'  {value} ({count})'))
    else:
        for value, count in sorted(category_counts.items()):
            color = color_maps[col][value]
            legend_elements.append(Patch(facecolor=color, edgecolor='white',
                                       label=f'  {value} ({count})'))

# Set limits and remove axes
ax.set_xlim(-5, 5)  # Reduced from (-7, 7) to make plot larger
ax.set_ylim(-5, 5)  # Reduced from (-7, 7) to make plot larger
ax.axis('off')

# Add legend
plt.legend(handles=legend_elements,
          loc='center left',
          bbox_to_anchor=(1.05, 0.5),  # Moved legend closer to plot (from 1.1 to 1.05)
          fontsize=8,
          frameon=False)

# Add title
plt.title('IO Trials Circular Visualization', pad=20, fontsize=16)

# Save as PNG with extra width for legend
plt.savefig('io_sunburst.png', bbox_inches='tight', dpi=600)
print("Visualization has been created and saved as 'io_sunburst.png'") 