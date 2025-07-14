import pandas as pd
import os

def analyze_ilc_tmb_distribution(file_path):
    """
    Performs a detailed analysis of the TMB score distribution for a given file.

    Args:
        file_path (str): The path to the CSV file.
    """
    print(f"--- Detailed TMB Analysis for: {file_path} ---")

    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return

    try:
        df = pd.read_csv(file_path, low_memory=False)

        # Ensure required columns exist
        if 'TMB Score' not in df.columns or 'Effective Patient ID' not in df.columns:
            print("Error: 'TMB Score' or 'Effective Patient ID' column not found.")
            return

        # --- Data Cleaning and Preparation ---
        # Convert 'TMB Score' to numeric, coercing errors
        df['TMB Score'] = pd.to_numeric(df['TMB Score'], errors='coerce')
        # Drop rows where TMB Score is not a number
        df.dropna(subset=['TMB Score'], inplace=True)
        # Get data for unique patients only (taking the first entry for each)
        unique_patients_df = df.drop_duplicates(subset=['Effective Patient ID'])

        # --- Calculations ---
        total_unique_patients = len(unique_patients_df)
        patients_tmb_gt_10 = unique_patients_df[unique_patients_df['TMB Score'] > 10]
        count_tmb_gt_10 = len(patients_tmb_gt_10)
        
        if total_unique_patients > 0:
            percentage_gt_10 = (count_tmb_gt_10 / total_unique_patients) * 100
        else:
            percentage_gt_10 = 0

        # --- Statistical Analysis ---
        stats = unique_patients_df['TMB Score'].describe()
        
        # --- Results ---
        print(f"\nTotal unique patients with a valid TMB score: {total_unique_patients}")
        print(f"Unique patients with TMB score > 10: {count_tmb_gt_10}")
        print(f"Percentage of unique patients with TMB score > 10: {percentage_gt_10:.2f}%\n")
        
        print("--- TMB Score Distribution ---")
        print(f"Mean (Average): {stats['mean']:.2f}")
        print(f"25th Percentile (Q1): {stats['25%']:.2f}")
        print(f"Median (50th Percentile): {stats['50%']:.2f}")
        print(f"75th Percentile (Q3): {stats['75%']:.2f}")
        print(f"Minimum: {stats['min']:.2f}")
        print(f"Maximum: {stats['max']:.2f}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    ilc_file_path = 'final_genomic_cohorts/ILC_genomic_infinity.csv'
    analyze_ilc_tmb_distribution(ilc_file_path) 