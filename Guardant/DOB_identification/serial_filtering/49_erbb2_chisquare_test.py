import numpy as np
from scipy.stats import chi2_contingency
import sys

def perform_chi_square_test_erbb2():
    """
    Performs a chi-square test on the provided ERBB2 mutation counts
    for ILC and IDC cohorts to test for significant differences.
    """
    try:
        # --- 1. Define the Contingency Table ---
        # Data format: [[Group1_Positive, Group1_Negative],
        #               [Group2_Positive, Group2_Negative]]
        
        # IDC: 82 positive, 715 - 82 = 633 negative
        # ILC: 21 positive, 123 - 21 = 102 negative
        
        observed_counts = np.array([
            [82, 633],  # IDC cohort
            [21, 102]     # ILC cohort
        ])
        
        # --- 2. Perform the Chi-Square Test ---
        chi2, p_value, dof, expected = chi2_contingency(observed_counts)
        
        # --- 3. Report the Findings ---
        print("--- Chi-Square Test for ERBB2 Mutation Rates (ILC vs. IDC) ---")
        
        print("\\nObjective: To determine if the ERBB2 mutation rates between cohorts are significantly different.")
        
        print("\\nContingency Table:")
        print("          ERBB2+    ERBB2-    Total")
        print(f"IDC       {observed_counts[0,0]:<9} {observed_counts[0,1]:<9} {np.sum(observed_counts[0,:]):<5}")
        print(f"ILC       {observed_counts[1,0]:<9} {observed_counts[1,1]:<9} {np.sum(observed_counts[1,:]):<5}")
        
        print("\\nTest Results:")
        print(f"Chi-Square Statistic: {chi2:.4f}")
        print(f"P-value: {p_value:.4f}")
        print(f"Degrees of Freedom: {dof}")
        
        print("\\n--- Conclusion ---")
        alpha = 0.05
        if p_value < alpha:
            print(f"The p-value ({p_value:.4f}) is less than the significance level ({alpha}).")
            print("We reject the null hypothesis.")
            print("This indicates that there is a statistically significant association between cohort type and ERBB2 mutation status.")
            
            # Compare proportions to confirm direction
            idc_proportion = observed_counts[0,0] / np.sum(observed_counts[0,:])
            ilc_proportion = observed_counts[1,0] / np.sum(observed_counts[1,:])
            direction = "higher" if ilc_proportion > idc_proportion else "lower"
            print(f"The mutation rate in the ILC cohort ({ilc_proportion:.2%}) is significantly {direction} than in the IDC cohort ({idc_proportion:.2%}).")
        else:
            print(f"The p-value ({p_value:.4f}) is greater than the significance level ({alpha}).")
            print("We fail to reject the null hypothesis.")
            print("This indicates that the observed difference in ERBB2 mutation rates between the ILC and IDC cohorts is not statistically significant.")
            
    except ImportError:
        print("Error: The 'scipy' library is required. Please install it using 'pip install scipy'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    perform_chi_square_test_erbb2() 