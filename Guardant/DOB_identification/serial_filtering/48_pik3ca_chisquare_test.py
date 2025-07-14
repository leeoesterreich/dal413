import numpy as np
from scipy.stats import chi2_contingency
import sys

def perform_chi_square_test():
    """
    Performs a chi-square test on the provided PIK3CA mutation counts
    for ILC and IDC cohorts to test for significant differences.
    """
    try:
        # --- 1. Define the Contingency Table ---
        # Data format: [[Group1_Positive, Group1_Negative],
        #               [Group2_Positive, Group2_Negative]]
        
        # IDC: 268 positive, 715 - 268 = 447 negative
        # ILC: 61 positive,  123 - 61  = 62 negative
        
        observed_counts = np.array([
            [268, 447],  # IDC cohort
            [61, 62]     # ILC cohort
        ])
        
        # --- 2. Perform the Chi-Square Test ---
        chi2, p_value, dof, expected = chi2_contingency(observed_counts)
        
        # --- 3. Report the Findings ---
        print("--- Chi-Square Test for PIK3CA Mutation Rates (ILC vs. IDC) ---")
        
        print("\\nObjective: To determine if the higher PIK3CA mutation rate in ILC patients is statistically significant.")
        
        print("\\nContingency Table:")
        print("          PIK3CA+   PIK3CA-   Total")
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
            print("This indicates that there is a statistically significant association between cohort type and PIK3CA mutation status.")
            
            # Compare proportions to confirm direction
            idc_proportion = observed_counts[0,0] / np.sum(observed_counts[0,:])
            ilc_proportion = observed_counts[1,0] / np.sum(observed_counts[1,:])
            print(f"The mutation rate in the ILC cohort ({ilc_proportion:.2%}) is significantly higher than in the IDC cohort ({idc_proportion:.2%}).")
        else:
            print(f"The p-value ({p_value:.4f}) is greater than the significance level ({alpha}).")
            print("We fail to reject the null hypothesis.")
            print("This indicates that the observed difference in PIK3CA mutation rates between the ILC and IDC cohorts is not statistically significant.")
            
    except ImportError:
        print("Error: The 'scipy' library is required. Please install it using 'pip install scipy'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    perform_chi_square_test() 