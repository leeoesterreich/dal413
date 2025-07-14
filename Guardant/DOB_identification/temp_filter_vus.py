import pandas as pd
import os

# Define paths
vus_dir = 'Vus_filtering'
idc_input = 'final_genomic_cohorts/IDC_genomic.csv'
ilc_input = 'final_genomic_cohorts/ILC_genomic.csv'
idc_output = os.path.join(vus_dir, 'IDC_genomic.csv')
ilc_output = os.path.join(vus_dir, 'ILC_genomic.csv')

# Filter IDC
try:
    idc_df = pd.read_csv(idc_input, low_memory=False)
    idc_df['Vus'] = idc_df['Vus'].astype(str).str.upper()
    idc_filtered = idc_df[idc_df['Vus'] == 'FALSE'].copy()
    idc_filtered.to_csv(idc_output, index=False)
    print(f"Filtered IDC data saved to {idc_output} ({len(idc_filtered)} rows)")
except FileNotFoundError:
    print(f"Could not find {idc_input}")

# Filter ILC
try:
    ilc_df = pd.read_csv(ilc_input, low_memory=False)
    ilc_df['Vus'] = ilc_df['Vus'].astype(str).str.upper()
    ilc_filtered = ilc_df[ilc_df['Vus'] == 'FALSE'].copy()
    ilc_filtered.to_csv(ilc_output, index=False)
    print(f"Filtered ILC data saved to {ilc_output} ({len(ilc_filtered)} rows)")
except FileNotFoundError:
    print(f"Could not find {ilc_input}") 