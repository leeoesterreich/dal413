import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def create_grouped_plot(df, sample_id_col, output_dir):
    """Create a grouped bar plot for Coverage and GB Dragen."""
    plt.figure(figsize=(20, 10))
    
    # Set the width of each bar and positions of the bars
    width = 0.35
    x = np.arange(len(df)) * 1.5  # Multiply by 1.5 to add more space between sample groups
    
    # Create the bars
    coverage_data = df["Avg Alignment Coverage (QC1)"]
    gb_dragen_data = df["GB DRAGEN"]
    
    # Create bars
    plt.bar(x - width/2, coverage_data, width, label='Coverage', 
            color=['red' if sample == 'TP22-M527' else '#1f77b4' for sample in df[sample_id_col]])
    plt.bar(x + width/2, gb_dragen_data, width, label='GB DRAGEN',
            color=['darkred' if sample == 'TP22-M527' else '#2ca02c' for sample in df[sample_id_col]])
    
    # Customize the plot
    plt.xlabel("Collab. Sample ID", fontsize=12)
    plt.ylabel("Value", fontsize=14)
    plt.title("Coverage and GB DRAGEN by Sample ID", fontsize=14)
    plt.xticks(x, df[sample_id_col], rotation=90, fontsize=10)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=12)
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = os.path.join(output_dir, "Coverage_and_GB_DRAGEN_plot.png")
    plt.savefig(plot_filename, dpi=300)
    plt.close()
    
    return "Coverage_and_GB_DRAGEN"

def create_plot(df, sample_id_col, metric_col, output_dir):
    """Create a bar plot for a specific metric."""
    plt.figure(figsize=(18, 10))
    
    # Create color array, set all to default color except TP22-M527
    colors = ['red' if x == 'TP22-M527' else '#1f77b4' for x in df[sample_id_col]]
    
    # Convert to millions if it's the unique mapped reads metric
    y_values = df[metric_col]
    ylabel = metric_col
    if metric_col == "Num Unq Mapped Reads":
        y_values = y_values / 1_000_000  # Convert to millions
        ylabel = "Number of Uniquely Mapped Reads (Millions)"
    
    plt.bar(df[sample_id_col], y_values, color=colors)
    plt.xlabel("Collab. Sample ID", fontsize=12)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(f"{metric_col} by Sample ID", fontsize=14)
    
    # Increase font size for tick labels
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=14)
    
    plt.tight_layout()

    # Save the plot
    safe_metric_name = metric_col.replace(" ", "_").replace("(", "").replace(")", "").replace("%", "percent").replace("/", "_per_")
    plot_filename = os.path.join(output_dir, f"{safe_metric_name}_plot.png")
    plt.savefig(plot_filename, dpi=300)
    plt.close()
    
    return safe_metric_name

def convert_to_numeric(series):
    """Convert a series to numeric, handling various formats."""
    # First try direct conversion
    numeric_series = pd.to_numeric(series, errors='coerce')
    
    # If we have any non-null values, return as is
    if not numeric_series.isna().all():
        return numeric_series
        
    # Try removing commas and converting
    try:
        return pd.to_numeric(series.astype(str).str.replace(",", ""), errors='coerce')
    except:
        return numeric_series

def main():
    # Define the path to the CSV file and output directory
    file_path = "/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/Extraction/ctDNA_shWGS.csv"
    output_dir = "plots"
    
    print(f"Reading file from: {file_path}")

    # Read the CSV file
    try:
        df = pd.read_csv(file_path)
        print(f"Successfully read CSV file with {len(df)} rows")
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    # Define the sample ID column and special metrics
    sample_id_col = "Collab. Sample ID"
    coverage_metric = "Avg Alignment Coverage (QC1)"
    gb_dragen_metric = "GB DRAGEN"
    efficiency_metric = "X Coverage / GB DRAGEN"
    
    # Initialize dictionary to store means
    means_dict = {}
    
    # Process each numeric column except the sample ID
    numeric_columns = []
    for col in df.columns:
        if col != sample_id_col:
            # Try to convert to numeric
            df[col] = convert_to_numeric(df[col])
            if not df[col].isna().all():  # Only include if at least some values converted successfully
                numeric_columns.append(col)
    
    print(f"\nProcessing metrics...")
    
    # First create the grouped plot for Coverage and GB DRAGEN
    print("\nCreating grouped plot for Coverage and GB DRAGEN")
    metric_df = df.dropna(subset=[coverage_metric, gb_dragen_metric])
    if not metric_df.empty:
        safe_metric_name = create_grouped_plot(metric_df, sample_id_col, output_dir)
        print(f"Plot saved as {safe_metric_name}_plot.png")
        
        # Calculate means for these metrics
        means_dict[coverage_metric] = metric_df[coverage_metric].mean()
        means_dict[gb_dragen_metric] = metric_df[gb_dragen_metric].mean()
    
    # Create separate plot for X Coverage / GB DRAGEN
    print(f"\nProcessing X Coverage / GB DRAGEN")
    metric_df = df.dropna(subset=[efficiency_metric])
    if not metric_df.empty:
        safe_metric_name = create_plot(metric_df, sample_id_col, efficiency_metric, output_dir)
        print(f"Plot saved as {safe_metric_name}_plot.png")
        means_dict[efficiency_metric] = metric_df[efficiency_metric].mean()
    
    # Process remaining metrics (excluding the special ones)
    remaining_metrics = [col for col in numeric_columns if col not in 
                        [coverage_metric, gb_dragen_metric, efficiency_metric]]
    
    for metric_col in remaining_metrics:
        print(f"\nProcessing metric: {metric_col}")
        metric_df = df.dropna(subset=[metric_col])
        
        if not metric_df.empty:
            safe_metric_name = create_plot(metric_df, sample_id_col, metric_col, output_dir)
            print(f"Plot saved as {safe_metric_name}_plot.png")
            
            mean_value = metric_df[metric_col].mean()
            if metric_col == "Num Unq Mapped Reads":
                mean_value = mean_value / 1_000_000
                means_dict[metric_col + " (Millions)"] = mean_value
            else:
                means_dict[metric_col] = mean_value
            print(f"Mean value: {mean_value:.4f}")
    
    # Create summary CSV
    summary_df = pd.DataFrame(list(means_dict.items()), columns=['Metric', 'Mean Value'])
    summary_csv_path = os.path.join(output_dir, 'metric_means_summary.csv')
    summary_df.to_csv(summary_csv_path, index=False)
    print(f"\nSummary of means saved to {summary_csv_path}")

if __name__ == "__main__":
    main()