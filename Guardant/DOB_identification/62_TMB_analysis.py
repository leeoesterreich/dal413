import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import os

def analyze_tmb_scores(ilc_file, idc_file, output_dir):
    """
    Analyzes and compares TMB scores from ILC and IDC genomic data,
    using only one representative test per unique patient.

    Args:
        ilc_file (str): Path to the ILC genomic infinity CSV file.
        idc_file (str): Path to the IDC genomic infinity CSV file.
        output_dir (str): Directory to save the report and plot.
    """
    # Load data
    try:
        ilc_df = pd.read_csv(ilc_file)
        idc_df = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error loading files: {e}")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    # --- 1. Filter for viable TMB Score and de-duplicate patients ---
    tmb_col = 'TMB Score'
    patient_col = 'Effective Patient ID'

    # ILC
    ilc_df[tmb_col] = pd.to_numeric(ilc_df[tmb_col], errors='coerce')
    ilc_viable_tmb = ilc_df.dropna(subset=[tmb_col])
    # De-duplicate to keep only the first test for each patient
    ilc_unique_df = ilc_viable_tmb.drop_duplicates(subset=[patient_col], keep='first')
    ilc_unique_patients = ilc_unique_df[patient_col].nunique()
    ilc_tmb_scores = ilc_unique_df[tmb_col]

    # IDC
    idc_df[tmb_col] = pd.to_numeric(idc_df[tmb_col], errors='coerce')
    idc_viable_tmb = idc_df.dropna(subset=[tmb_col])
    # De-duplicate to keep only the first test for each patient
    idc_unique_df = idc_viable_tmb.drop_duplicates(subset=[patient_col], keep='first')
    idc_unique_patients = idc_unique_df[patient_col].nunique()
    idc_tmb_scores = idc_unique_df[tmb_col]

    # --- 2. Perform statistical test ---
    # Mann-Whitney U test is appropriate as we don't assume a normal distribution.
    if not ilc_tmb_scores.empty and not idc_tmb_scores.empty:
        stat, p_value = mannwhitneyu(ilc_tmb_scores, idc_tmb_scores, alternative='two-sided')
        ilc_median = ilc_tmb_scores.median()
        idc_median = idc_tmb_scores.median()
    else:
        stat, p_value, ilc_median, idc_median = (None, None, None, None)

    # --- Generate Report ---
    report_path = os.path.join(output_dir, 'tmb_analysis_report_unique_patients.txt')
    with open(report_path, 'w') as f:
        f.write("TMB Score Analysis Report (Unique Patients Only)\n")
        f.write("="*50 + "\n\n")
        f.write("Method: Each patient is represented by their first test with a valid TMB score.\n\n")
        f.write(f"ILC Data File: {os.path.basename(ilc_file)}\n")
        f.write(f"Number of unique patients with viable TMB score: {ilc_unique_patients}\n\n")
        f.write(f"IDC Data File: {os.path.basename(idc_file)}\n")
        f.write(f"Number of unique patients with viable TMB score: {idc_unique_patients}\n\n")
        f.write("Statistical Comparison\n")
        f.write("-" * 30 + "\n")
        f.write("Test: Mann-Whitney U test\n")
        if stat is not None:
            f.write(f"U-statistic: {stat:.4f}\n")
            f.write(f"P-value: {p_value:.4f}\n")
            f.write(f"ILC Median TMB Score: {ilc_median:.4f}\n")
            f.write(f"IDC Median TMB Score: {idc_median:.4f}\n")
            if p_value < 0.05:
                higher_group = "IDC" if idc_median > ilc_median else "ILC"
                f.write(f"Conclusion: There is a statistically significant difference between the TMB scores, with the {higher_group} group having a higher median TMB score.\n")
            else:
                f.write("Conclusion: There is no statistically significant difference between the TMB scores of the ILC and IDC groups.\n")
        else:
            f.write("Could not perform statistical test due to lack of viable TMB data in one or both groups.\n")

    print(f"Analysis report saved to {report_path}")

    # --- 3. Plot boxplot ---
    if ilc_tmb_scores.empty and idc_tmb_scores.empty:
        print("Cannot generate plot: No viable TMB data available.")
        return
        
    plt.figure(figsize=(8, 6))

    # Combine data for plotting
    plot_data = pd.concat([
        pd.DataFrame({'TMB Score': ilc_tmb_scores, 'Group': 'ILC'}),
        pd.DataFrame({'TMB Score': idc_tmb_scores, 'Group': 'IDC'})
    ])
    
    # Apply log10(x+1) transformation
    plot_data['log10_TMB'] = np.log10(plot_data['TMB Score'] + 1)


    # Boxplot
    sns.boxplot(x='Group', y='log10_TMB', data=plot_data, palette=['#D3D3D3', '#A9A9A9'], showfliers=False)
    # Swarmplot/stripplot for individual data points
    sns.stripplot(x='Group', y='log10_TMB', data=plot_data, color='lightgrey', jitter=True, alpha=0.6, edgecolor='gray', linewidth=0.5)

    plt.title('Comparison of log10(TMB Scores) between ILC and IDC Groups (Unique Patients)')
    plt.ylabel('log10(TMB Score + 1)')
    plt.xlabel('Group')
    plt.tight_layout()

    plot_path = os.path.join(output_dir, 'tmb_comparison_boxplot_logscale_unique_patients.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Boxplot saved to {plot_path}")


if __name__ == '__main__':
    # Define file paths
    base_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts'
    ilc_infinity_file = os.path.join(base_path, 'ILC_genomic_infinity.csv')
    idc_infinity_file = os.path.join(base_path, 'IDC_genomic_infinity.csv')

    # Create output directory if it doesn't exist
    if not os.path.exists(base_path):
        os.makedirs(base_path)

    # Run the analysis
    analyze_tmb_scores(ilc_infinity_file, idc_infinity_file, base_path) 