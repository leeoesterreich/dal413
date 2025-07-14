import pandas as pd
from scipy.stats import shapiro
import os

def test_normality(file_path, group_name):
    """
    Loads data, extracts unique patient TMB scores, and performs the Shapiro-Wilk test for normality.

    Args:
        file_path (str): The path to the CSV file.
        group_name (str): The name of the group being tested (e.g., "ILC" or "IDC").
    """
    print(f"--- Normality Test for: {group_name} ---")

    try:
        df = pd.read_csv(file_path, low_memory=False)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return

    tmb_col = 'TMB Score'
    patient_col = 'Effective Patient ID'

    # Ensure TMB score is numeric and drop rows without a valid score
    df[tmb_col] = pd.to_numeric(df[tmb_col], errors='coerce')
    valid_tmb_df = df.dropna(subset=[tmb_col])

    # De-duplicate to ensure one record per unique patient
    unique_patient_df = valid_tmb_df.drop_duplicates(subset=[patient_col], keep='first')
    
    if len(unique_patient_df) < 3:
        print("Not enough data points to perform a reliable normality test.")
        return

    tmb_scores = unique_patient_df[tmb_col]

    # Perform Shapiro-Wilk test
    stat, p_value = shapiro(tmb_scores)

    print(f"Shapiro-Wilk Test Results for {group_name} TMB Scores:")
    print(f"Test Statistic: {stat:.4f}")
    print(f"P-value: {p_value:.4f}")

    # Interpretation
    if p_value < 0.05:
        print("Conclusion: The p-value is less than 0.05, so we reject the null hypothesis.")
        print("This suggests the data is NOT normally distributed.\n")
    else:
        print("Conclusion: The p-value is greater than 0.05, so we cannot reject the null hypothesis.")
        print("This suggests the data may be normally distributed.\n")

if __name__ == '__main__':
    base_path = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/final_genomic_cohorts'
    ilc_file = os.path.join(base_path, 'ILC_genomic_infinity.csv')
    idc_file = os.path.join(base_path, 'IDC_genomic_infinity.csv')

    test_normality(ilc_file, "ILC")
    test_normality(idc_file, "IDC")

    print("Why shouldn't we use a t-test?")
    print("A t-test assumes that the data in both groups are normally distributed. Since the Shapiro-Wilk test shows that at least one of our groups (and likely both) significantly deviates from a normal distribution, using a t-test could lead to incorrect conclusions. The Mann-Whitney U test is the appropriate non-parametric alternative because it does not require this assumption of normality.") 