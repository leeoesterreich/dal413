import pandas as pd
import numpy as np
import os
import argparse

def load_and_compare_coefficients(paper_coeffs_path, your_model_coeffs_path, rna_signature_name_to_compare):
    """
    Loads coefficients from the main paper and your single specified model,
    then compares them for the target RNA signature.
    """
    print(f"Loading paper's master coefficient file from: {paper_coeffs_path}")
    try:
        paper_df = pd.read_csv(paper_coeffs_path)
        if paper_df.columns[0] == 'Unnamed: 0' or paper_df.columns[0] == '':
            paper_df = paper_df.rename(columns={paper_df.columns[0]: 'Signature_Name'})
        paper_df = paper_df.set_index('Signature_Name')
        print(f"Paper's coefficient data loaded successfully. Shape: {paper_df.shape}")
        # print(f"Paper's columns (features, first 5): {paper_df.columns.tolist()[:5]}...")
    except FileNotFoundError:
        print(f"ERROR: Paper's coefficient file not found at {paper_coeffs_path}")
        return
    except Exception as e:
        print(f"ERROR: Could not load or process paper's CSV: {e}")
        return

    print(f"\n--- Comparing coefficients for RNA Signature: {rna_signature_name_to_compare} ---")
    
    print(f"  Loading your model's coefficients from: {your_model_coeffs_path}")
    try:
        your_df = pd.read_csv(your_model_coeffs_path)
        
        feature_col_from_your_file = 'FeatureName'
        coeff_col_from_your_file = 'Weight'

        if feature_col_from_your_file not in your_df.columns:
            print(f"  ERROR: Expected feature column '{feature_col_from_your_file}' not found in {your_model_coeffs_path}. Available: {your_df.columns.tolist()}")
            return
        if coeff_col_from_your_file not in your_df.columns:
            print(f"  ERROR: Expected coefficient column '{coeff_col_from_your_file}' not found in {your_model_coeffs_path}. Available: {your_df.columns.tolist()}")
            return
        
        print(f"    Using column '{feature_col_from_your_file}' for features and '{coeff_col_from_your_file}' for coefficients (assumed original scale).")

        your_df = your_df[[feature_col_from_your_file, coeff_col_from_your_file]].rename(
            columns={feature_col_from_your_file: 'Feature', coeff_col_from_your_file: 'Your_Coefficient'}
        )
        your_df = your_df.set_index('Feature')
        your_df_non_zero = your_df[your_df['Your_Coefficient'] != 0].copy()
        print(f"    Your model: Found {len(your_df)} features, {len(your_df_non_zero)} non-zero coefficients.")

    except FileNotFoundError:
        print(f"  ERROR: Your coefficient file not found: {your_model_coeffs_path}")
        return
    except Exception as e:
        print(f"  ERROR: Could not load or process your coefficient file {your_model_coeffs_path}: {e}")
        return

    if rna_signature_name_to_compare not in paper_df.index:
        print(f"  ERROR: RNA Signature '{rna_signature_name_to_compare}' not found in the paper's coefficient file.")
        print(f"  Available signatures in paper file start with: {paper_df.index.tolist()[:5]}...")
        return
    
    paper_coeffs_series = paper_df.loc[rna_signature_name_to_compare].copy()
    paper_coeffs_series.name = 'Paper_Coefficient'
    paper_coeffs_series.index = paper_coeffs_series.index.str.strip().str.replace(r'\s+', ' ', regex=True).str.replace(u'\xa0', u' ')
    paper_coeffs_non_zero = paper_coeffs_series[paper_coeffs_series != 0].copy()
    print(f"    Paper model: Found {len(paper_coeffs_series)} features, {len(paper_coeffs_non_zero)} non-zero coefficients.")

    your_df_non_zero.index = your_df_non_zero.index.str.strip().str.replace(r'\s+', ' ', regex=True).str.replace(u'\xa0', u' ')

    paper_coeffs_df_for_merge = paper_coeffs_non_zero.reset_index()
    paper_coeffs_df_for_merge.columns = ['Feature', 'Paper_Coefficient']
    
    your_coeffs_df_for_merge = your_df_non_zero.reset_index()
    
    comparison_df = pd.merge(paper_coeffs_df_for_merge, your_coeffs_df_for_merge, on='Feature', how='outer')
    comparison_df = comparison_df.set_index('Feature')
    comparison_df = comparison_df.fillna(0)

    comparison_df['Ratio (Paper/Your)'] = np.where(
        (comparison_df['Paper_Coefficient'] != 0) & (comparison_df['Your_Coefficient'] != 0),
        comparison_df['Paper_Coefficient'] / comparison_df['Your_Coefficient'],
        np.nan
    )
    comparison_df['Abs_Difference'] = (comparison_df['Your_Coefficient'] - comparison_df['Paper_Coefficient']).abs()

    print("\n    --- Comparison of Non-Zero Coefficients (Top 40 by Abs Difference, then by Paper's Coeff Mag): ---")
    sorted_comparison = comparison_df.copy()
    sorted_comparison = sorted_comparison[(sorted_comparison['Paper_Coefficient'] != 0) | (sorted_comparison['Your_Coefficient'] != 0)]
    display_columns = ['Paper_Coefficient', 'Your_Coefficient', 'Ratio (Paper/Your)', 'Abs_Difference']
    sorted_comparison['_abs_paper_coeff'] = sorted_comparison['Paper_Coefficient'].abs()
    sorted_comparison = sorted_comparison.sort_values(by=['Abs_Difference', '_abs_paper_coeff'], ascending=[False, False])
    del sorted_comparison['_abs_paper_coeff']
    
    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
        print(sorted_comparison[display_columns].head(40).to_string())

    common_features_both_non_zero = comparison_df[
        (comparison_df['Paper_Coefficient'] != 0) & (comparison_df['Your_Coefficient'] != 0)
    ].copy()

    if not common_features_both_non_zero.empty:
        sum_abs_paper = common_features_both_non_zero['Paper_Coefficient'].abs().sum()
        sum_abs_your = common_features_both_non_zero['Your_Coefficient'].abs().sum()
        
        print("\n    --- Summary for Common Non-Zero Features ---")
        print(f"    Number of common non-zero features: {len(common_features_both_non_zero)}")
        if sum_abs_your != 0:
            overall_magnitude_ratio_paper_over_your = sum_abs_paper / sum_abs_your
            print(f"    Ratio of Sum(|Paper Coefficients|) / Sum(|Your Coefficients|): {overall_magnitude_ratio_paper_over_your:.4f}")
        else:
            print("    Your model's sum of absolute coefficients is zero, cannot compute overall magnitude ratio (Paper/Your).")
        
        median_individual_ratio = common_features_both_non_zero['Ratio (Paper/Your)'].median()
        non_nan_ratios = common_features_both_non_zero['Ratio (Paper/Your)'].dropna()
        num_ratios_for_median = len(non_nan_ratios)
        print(f"    Median of (Paper Coeff / Your Coeff) for common non-zeros: {median_individual_ratio:.4f} (calculated over {num_ratios_for_median} features)")

        sign_differs = common_features_both_non_zero[
            np.sign(common_features_both_non_zero['Paper_Coefficient']) != np.sign(common_features_both_non_zero['Your_Coefficient'])
        ]
        print(f"    Number of common non-zero features where sign differs: {len(sign_differs)}")
        if 0 < len(sign_differs) <= 10:
             print("      Features with differing signs:")
             with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
                print(sign_differs[['Paper_Coefficient', 'Your_Coefficient']].to_string())
    else:
        print("\n    No common non-zero features found to compare magnitudes directly.")

    print(f"--- End of comparison for {rna_signature_name_to_compare} ---\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare your model's coefficients against a paper's reference coefficients.")
    parser.add_argument("your_coeffs_path", type=str, help="Path to your model's generated coefficient CSV file (e.g., containing FeatureName and Weight).")
    parser.add_argument("--paper_coeffs_path", type=str, 
                        default="/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/adjusted_output/41467_2019_13588_MOESM4_ESM_Elastic_Net_gene_signatures.csv",
                        help="Path to the paper's master coefficient CSV file.")
    parser.add_argument("--rna_signature_name", type=str, 
                        default="GSEA_Median_GP17_Basal_signaling.r=0.958_SMID_BREAST_CANCER_BASAL_UP", 
                        help="The specific RNA signature name to extract from the paper's CSV and compare.")
    
    args = parser.parse_args()

    if not os.path.isfile(args.your_coeffs_path):
        print(f"ERROR: Your coefficient file not found at: {args.your_coeffs_path}")
    elif not os.path.isfile(args.paper_coeffs_path):
        print(f"ERROR: Paper's coefficient file not found at: {args.paper_coeffs_path}")
    else:
        print(f"Comparing your model ({args.your_coeffs_path}) with paper ({args.paper_coeffs_path}) for signature: '{args.rna_signature_name}'")
        load_and_compare_coefficients(args.paper_coeffs_path, args.your_coeffs_path, args.rna_signature_name) 