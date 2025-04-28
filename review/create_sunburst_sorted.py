import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Set the backend to Agg for non-interactive environments
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Wedge
from itertools import groupby

# Read the sorted CSV file
df = pd.read_csv('IO_summary_sorted.csv')

# Set up the figure
fig, ax = plt.subplots(figsize=(18, 13))  # Reduced from (20, 15) to make plot fill more space
ax.set_aspect('equal')

# Base colors for each level - using a soft rainbow palette
base_colors = {
    'Biomarker': '#FF69B4',      # Hot pink (keeping as is)
    'Molecular subtype': '#FF6B6B',  # Soft red
    'Stage': '#FFA07A',          # Light salmon (soft orange)
    'Immunotherapy Drug': '#FFD700',  # Gold (soft yellow)
    'Primary Endpoint': '#98FB98',  # Pale green
    'Stratification': '#87CEEB'   # Sky blue
}

# Create drug name mapping
drug_full_names = {
    'pembro': 'pembrolizumab',
    'atezo': 'atezolizumab',
    'nivo': 'nivolumab',
    'durva': 'durvalumab',
    'ipi': 'ipilimumab',
    'nivo+ipi': 'nivolumab, ipilimumab',
    'camre': 'camrelizumab'  # Added camrelizumab
}

# Create endpoint name mapping
endpoint_full_names = {
    'ORR': 'Overall response rate',
    'PFS': 'Progression free survival',
    'pCR': 'pathological\ncomplete response',
    'pCR, EFS': 'Event free survival',
    'PFS, OS': 'Overall survival'
}

def get_color(base_color, index, total):
    """Generate a color with more contrast between shades"""
    # Convert hex to RGB
    r = int(base_color[1:3], 16)
    g = int(base_color[3:5], 16)
    b = int(base_color[5:7], 16)
    
    # Calculate shade with more contrast
    # Reverse the order so earlier segments are more saturated
    shade = 1.0 - (0.8 * (index / total))  # Range from 1.0 to 0.2
    
    # Apply shade
    r = int(r * shade)
    g = int(g * shade)
    b = int(b * shade)
    
    # Convert back to hex
    return f'#{r:02x}{g:02x}{b:02x}'

def create_color_variations(base_color, n):
    rgb = matplotlib.colors.hex2color(base_color)
    variations = []
    for i in range(n):
        # Reverse the order so earlier segments are more saturated
        factor = 1.0 - (0.7 * i / (n-1 if n > 1 else 1))  # Range from 1.0 to 0.3
        variation = tuple(min(1.0, c * factor + (1 - factor) * 0.8) for c in rgb)
        variations.append(variation)
    return variations

def get_text_rotation(angle):
    rotation = angle + 90
    if rotation > 90 and rotation < 270:
        rotation = rotation - 180
    return rotation

# Get total number of trials
total_trials = len(df)
angle_per_trial = 360 / total_trials
spacing_angle = 4  # Reduced from 6 to 4 for smaller gaps between segments

# Calculate the rotation needed to make PD-L1 horizontal
# Find the index of PD-L1 in the data
pdl1_index = df[df['Biomarker'] == 'PD-L1'].index[0]
# Calculate the angle where PD-L1 appears
pdl1_angle = pdl1_index * angle_per_trial
# Calculate the rotation needed to make it horizontal (0 degrees)
initial_rotation = -pdl1_angle

# Calculate dimensions
inner_radius = 1.2  # Reduced from 1.4
outer_radius = inner_radius + 0.6  # Increased from 0.4 to 0.6 for wider biomarker circle
circle_spacing = 0.2  # Reduced from 0.3 to 0.2 for smaller spacing between circles
circle_width = 0.3    # Consistent width for all outer circles

# Remove the sorting step to maintain original CSV order
# df = df.sort_values(['Biomarker', 'Molecular subtype', 'Stage', 'IO Drug', 'Primary Endpoint', 'Stratification'])

# Create color maps for each column
color_maps = {}
for col in df.columns:
    if col != 'Trial':
        # Get unique values in the order they appear in the CSV
        unique_values = []
        seen = set()
        for val in df[col]:
            val_str = str(val).strip() if not pd.isna(val) else 'NA'
            if val_str not in seen:
                unique_values.append(val_str)
                seen.add(val_str)
        colors = create_color_variations(base_colors[col], len(unique_values))
        color_maps[col] = dict(zip(unique_values, colors))

# Plot each level
for i, col in enumerate(['Biomarker', 'Molecular subtype', 'Stage', 'Immunotherapy Drug', 'Primary Endpoint', 'Stratification']):
    # Calculate radius for this level
    if col == 'Biomarker':
        current_inner = inner_radius
        current_outer = outer_radius
    else:
        # Use consistent spacing and width for all outer circles
        current_inner = outer_radius + circle_spacing + (i-1) * (circle_spacing + circle_width)
        current_outer = current_inner + circle_width

    # Group consecutive segments with the same value
    current_angle = initial_rotation
    values = [str(val).strip() if not pd.isna(val) else 'NA' for val in df[col]]
    
    # Group consecutive values
    for value, group in groupby(enumerate(values), lambda x: x[1]):
        # Convert group iterator to list to get count
        group_indices = list(group)
        count = len(group_indices)
        segment_angle = angle_per_trial * count
        
        color = color_maps[col].get(value, base_colors[col])
        
        # Draw the aggregated wedge
        wedge = Wedge(
            center=(0, 0),
            r=current_outer,
            width=current_outer-current_inner,
            theta1=current_angle + spacing_angle/2,
            theta2=current_angle + segment_angle - spacing_angle/2,
            facecolor=color,
            edgecolor='white',
            linewidth=1
        )
        ax.add_patch(wedge)
        
        # Add text label if segment is large enough
        if segment_angle > 10:
            text_angle = current_angle + segment_angle/2
            text_radius = (current_inner + current_outer) / 2
            text_x = text_radius * np.cos(np.radians(text_angle))
            text_y = text_radius * np.sin(np.radians(text_angle))
            
            rotation = get_text_rotation(text_angle)
            text_color = 'black'
            # Use larger font size for biomarker circle (inner circle)
            font_size = 14 if col == 'Biomarker' else 9  # Decreased from 10 to 9 for outer circles
            plt.text(text_x, text_y, value,
                    rotation=rotation,
                    ha='center', va='center',
                    fontsize=font_size,
                    weight='bold',
                    color=text_color)
        
        current_angle += segment_angle

# Add trial names as the outermost circle
current_inner = current_outer + circle_spacing
current_outer = current_inner + circle_width

# Group consecutive segments with the same value
current_angle = initial_rotation
values = [str(val).strip() if not pd.isna(val) else 'NA' for val in df['Trial']]

# Group consecutive values
for value, group in groupby(enumerate(values), lambda x: x[1]):
    # Convert group iterator to list to get count
    group_indices = list(group)
    count = len(group_indices)
    segment_angle = angle_per_trial * count
    
    # Use light grey for trial segments
    wedge = Wedge(
        center=(0, 0),
        r=current_outer,
        width=current_outer-current_inner,
        theta1=current_angle + spacing_angle/2,
        theta2=current_angle + segment_angle - spacing_angle/2,
        facecolor='#F0F0F0',  # Light grey
        edgecolor='white',
        linewidth=1
    )
    ax.add_patch(wedge)
    
    # Add text label if segment is large enough
    if segment_angle > 10:
        text_angle = current_angle + segment_angle/2
        text_radius = (current_inner + current_outer) / 2
        text_x = text_radius * np.cos(np.radians(text_angle))
        text_y = text_radius * np.sin(np.radians(text_angle))
        
        rotation = get_text_rotation(text_angle)
        
        # Handle special cases for trial names
        if value == "DART, SWOG S1609":
            value = "DART, SWOG\nS1609"
        elif value == "SAFIR-02-Breast IMMUNO":
            value = "SAFIR-02-\nBreast IMMUNO"
        
        plt.text(text_x, text_y, value,
                rotation=rotation,
                ha='center', va='center',
                fontsize=8,  # Decreased from 10 to 8 to match other outer circles
                color='black',
                weight='bold')
    
    current_angle += segment_angle

# Create legend elements
legend_elements = []
for col in ['Biomarker', 'Molecular subtype', 'Stage', 'Immunotherapy Drug', 'Primary Endpoint', 'Stratification']:
    # Add category label with definition for Stratification
    if col == 'Stratification':
        legend_elements.append(Patch(facecolor='white', edgecolor='none', 
            label=f'\n{col}:\n(primary endpoint differences\nof biomarker level high vs. low)'))
    else:
        legend_elements.append(Patch(facecolor='white', edgecolor='none', label=f'\n{col}:'))
    
    # Add color variations for this category
    unique_values = []
    seen = set()
    
    # Add to unique values in order of appearance
    for val in df[col]:
        val_str = str(val).strip() if not pd.isna(val) else 'NA'
        if val_str not in seen:
            unique_values.append(val_str)
            seen.add(val_str)
    
    # Add each unique value with its color
    for value in unique_values:
        color = color_maps[col].get(value, base_colors[col])
        # Add full drug names for IO Drug category
        if col == 'Immunotherapy Drug' and value in drug_full_names:
            if value == 'nivo+ipi':
                label = f"{value}\n({drug_full_names[value]})"
            else:
                label = f"{value} ({drug_full_names[value]})"
        # Add endpoint definitions for Primary Endpoint category
        elif col == 'Primary Endpoint' and value in endpoint_full_names:
            label = f"{value} ({endpoint_full_names[value]})"
        else:
            label = value
        legend_elements.append(Patch(facecolor=color, edgecolor='white', label=label))

# Add legend
ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

# Remove axes
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# Set limits to ensure the plot is centered and fills the space
max_radius = current_outer + 0.1  # Reduced from 0.2 to 0.1
ax.set_xlim(-max_radius, max_radius)
ax.set_ylim(-max_radius, max_radius)

# Add a light pink center circle
center_circle = plt.Circle((0, 0), inner_radius - 0.2, color='#FFE4E1', alpha=0.5)  # Misty Rose color
ax.add_patch(center_circle)

# Save the plot in both PNG and PDF formats
plt.savefig('sunburst_sorted.png', bbox_inches='tight', dpi=300, pad_inches=0.05)  # Reduced padding from 0.1 to 0.05
plt.savefig('sunburst_sorted.pdf', bbox_inches='tight', pad_inches=0.05)  # Reduced padding from 0.1 to 0.05
plt.close()

print("Sunburst plot has been created and saved as 'sunburst_sorted.png' and 'sunburst_sorted.pdf'") 