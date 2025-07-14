import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime

# Define the input file and output directory
INPUT_CSV_FILE = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/results_xgb_retrained_60_20_20_split/summaries/all_signatures_retrained_results.csv'
OUTPUT_DIR = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/parameter_tuning/normalization_and_model_training/results_xgb_retrained_60_20_20_split/summaries/'
OUTPUT_PLOT_FILENAME = 'retrained_tercile_auc_boxplot_final.png'

def load_metrics(file_path):
    """Load model metrics from the CSV file."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Metrics file not found: {file_path}")
    
    metrics_df = pd.read_csv(file_path)
    
    # Columns for this specific CSV
    required_cols = ['Signature', 'Train_AUC_Tercile', 'Test_AUC_Tercile', 'Val_AUC_Tercile'] 
    for col in required_cols:
        if col not in metrics_df.columns:
            raise ValueError(f"Required column '{col}' not found in the CSV file.")
            
    print(f"Loaded metrics for {len(metrics_df)} signatures from {os.path.basename(file_path)}")
    return metrics_df

def create_auc_boxplot(metrics_df, output_dir, output_filename):
    """Create box plots for Train, Test, and Validation Tercile AUCs."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    output_path = os.path.join(output_dir, output_filename)

    plt.figure(figsize=(8, 8))
    
    plot_data = []
    # AUC columns for this specific CSV and desired plot order
    auc_columns = {
        'Train_AUC_Tercile': 'Training', # Label for plot
        'Test_AUC_Tercile': 'Testing',   # Label for plot
        'Val_AUC_Tercile': 'Validation'  # Label for plot
    }
    
    # Ensure correct order for processing and plotting
    # The order of items in auc_columns will determine processing if iterated directly.
    # For explicit plot order, use xtick_order defined later.
    
    for col_name, type_label in auc_columns.items():
        for _, row in metrics_df.iterrows():
            # Use 'Signature' as the column name from this CSV
            if pd.notna(row[col_name]) and 'Signature' in row and pd.notna(row['Signature']):
                plot_data.append({
                    'Type': type_label,
                    'AUC': row[col_name],
                    'Signature': row['Signature'] 
                })
            else:
                # Print signature name if available, otherwise indicate missing signature info
                sig_name = row['Signature'] if 'Signature' in row and pd.notna(row['Signature']) else "Unknown Signature"
                print(f"Warning: NaN value or missing signature found for {sig_name} in column {col_name}. Skipping this entry for plot.")

    if not plot_data:
        print("No data available for plotting after filtering NaNs. Exiting plot generation.")
        plt.close()
        return

    plot_df = pd.DataFrame(plot_data)
    
    # Define the explicit order for x-axis categories
    xtick_order = ['Training', 'Testing', 'Validation']
    
    sns.boxplot(data=plot_df, x='Type', y='AUC', 
                order=xtick_order,
                color='white',
                boxprops=dict(edgecolor='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                width=0.4,
                showfliers=False)
    
    sns.stripplot(data=plot_df, x='Type', y='AUC',
                  order=xtick_order,
                  color='darkgrey',
                  size=4,
                  alpha=0.5,
                  jitter=True,
                  dodge=True) # dodge=True is important if you have hue, not strictly needed here but good practice

    plt.axhline(y=0.75, color='gray', linestyle='--', alpha=0.5)
    
    highlight_signatures = {
        'UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450': 'red',
        'GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP': 'blue',
        'GSEA_Median_GP7_Estrogen_signaling.r=0.97_SMID_BREAST_CANCER_BASAL_DN': 'darkgreen' 
    }
    
    legend_handles = []
    legend_labels = []
    
    for signature_name_part, color in highlight_signatures.items():
        # Ensure we are looking for 'Signature' column in plot_df
        mask = plot_df['Signature'].str.contains(signature_name_part, na=False)
        if any(mask):
            sig_data = plot_df[mask]
            point_plotted_for_legend = False
            # Iterate based on the defined xtick_order for correct dot placement
            for i, type_label_in_order in enumerate(xtick_order):
                # Filter sig_data for the current type in the defined order
                auc_value_series = sig_data[sig_data['Type'] == type_label_in_order]['AUC']
                if not auc_value_series.empty:
                    for value in auc_value_series.values: # Should be a single value if data is unique per sig/type
                        plt.plot(i, value, 'o', color=color, markersize=8, 
                               alpha=0.9, markeredgecolor='black', label=f'_{signature_name_part}_point_{type_label_in_order}', zorder=5) 
                    
                    if not point_plotted_for_legend: 
                        if 'UNC_RB_LOH' in signature_name_part:
                            label = 'RB_LOH'
                        elif 'Basal_signaling' in signature_name_part:
                            label = 'Basal signaling'
                        elif 'Estrogen_signaling' in signature_name_part:
                            label = 'Estrogen signaling'
                        else:
                            label = signature_name_part 
                        
                        if label not in legend_labels:
                             legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=label,
                                markerfacecolor=color, markersize=8, markeredgecolor='black'))
                             legend_labels.append(label)
                        point_plotted_for_legend = True

    plt.title('XGBoost ROC AUC', fontsize=16, pad=20)
    plt.xlabel('', fontsize=14) 
    plt.ylabel('Tercile AUC Score', fontsize=14)
    
    # Set custom x-tick labels based on xtick_order
    # plt.xticks(ticks=range(len(xtick_order)), labels=xtick_order, fontsize=12) # This was how it was before
    # Ensure current_ticks aligns with the number of categories if using get_xticks()
    # Simpler to just set them directly if order is fixed
    new_xtick_labels = ['Training(60%)', 'Testing(20%)', 'Validation(20%)']
    plt.xticks(ticks=np.arange(len(xtick_order)), labels=new_xtick_labels, fontsize=12)


    plt.yticks(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7, axis='y') 
    
    if legend_handles:
        plt.legend(handles=legend_handles, labels=legend_labels, loc='lower right',
                   fontsize=10, title="Highlighted Signatures", title_fontsize=11)
    
    plt.ylim(0.4, 1.0) # Updated Y-axis limits
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.95]) 
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nBoxplot saved to: {output_path}")

def main():
    try:
        metrics_data = load_metrics(INPUT_CSV_FILE)
        # Check for the 'Error' column and filter out rows where it's not NaN (i.e., an error occurred)
        if 'Error' in metrics_data.columns:
            original_count = len(metrics_data)
            metrics_data = metrics_data[metrics_data['Error'].isna()]
            errors_filtered_count = original_count - len(metrics_data)
            if errors_filtered_count > 0:
                print(f"Filtered out {errors_filtered_count} signatures with errors during training.")

        create_auc_boxplot(metrics_data, OUTPUT_DIR, OUTPUT_PLOT_FILENAME)
        print("\nVisualization complete!")
        
    except FileNotFoundError as e:
        print(f"Error: {str(e)}")
        print("Please ensure the input CSV file exists at the specified path.")
    except ValueError as e:
        print(f"Error: {str(e)}")
        print("Please check the CSV file format and column names.")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    main() 