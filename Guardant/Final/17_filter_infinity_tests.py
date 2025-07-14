import pandas as pd
import numpy as np

# Read the input files
idc_df = pd.read_csv('IDC_genomic_corrected.csv')
ilc_df = pd.read_csv('ILC_genomic_corrected.csv')

# Filter for Guardant360 Infinity tests
idc_infinity = idc_df[idc_df['Test Reported'] == 'Guardant360 Infinity']
ilc_infinity = ilc_df[ilc_df['Test Reported'] == 'Guardant360 Infinity']

# Save filtered data
idc_infinity.to_csv('IDC_genomic_infinity_corrected.csv', index=False)
ilc_infinity.to_csv('ILC_genomic_infinity_corrected.csv', index=False)

# Generate summaries
def generate_summary(df, cohort):
    summary = {
        f'Total {cohort} Infinity Tests': len(df),
        f'Unique {cohort} Patients': len(df['Effective Patient ID'].unique()),
        'Gene Mutations': {
            'PIK3CA': len(df[df['Gene'] == 'PIK3CA']),
            'ESR1': len(df[df['Gene'] == 'ESR1']),
            'ERBB2': len(df[df['Gene'] == 'ERBB2']),
            'BRCA1': len(df[df['Gene'] == 'BRCA1']),
            'BRCA2': len(df[df['Gene'] == 'BRCA2']),
            'PALB2': len(df[df['Gene'] == 'PALB2'])
        }
    }
    
    # Convert to DataFrame
    summary_df = pd.DataFrame([
        ['Total Tests', summary[f'Total {cohort} Infinity Tests']],
        ['Unique Patients', summary[f'Unique {cohort} Patients']],
        ['PIK3CA Mutations', summary['Gene Mutations']['PIK3CA']],
        ['ESR1 Mutations', summary['Gene Mutations']['ESR1']],
        ['ERBB2 Mutations', summary['Gene Mutations']['ERBB2']],
        ['BRCA1 Mutations', summary['Gene Mutations']['BRCA1']],
        ['BRCA2 Mutations', summary['Gene Mutations']['BRCA2']],
        ['PALB2 Mutations', summary['Gene Mutations']['PALB2']]
    ], columns=['Metric', 'Count'])
    
    return summary_df

# Generate and save summaries
idc_summary = generate_summary(idc_infinity, 'IDC')
ilc_summary = generate_summary(ilc_infinity, 'ILC')

idc_summary.to_csv('IDC_infinity_summary_corrected.csv', index=False)
ilc_summary.to_csv('ILC_infinity_summary_corrected.csv', index=False)

# Print summaries to console
print("\nIDC Infinity Summary:")
print(idc_summary.to_string(index=False))
print("\nILC Infinity Summary:")
print(ilc_summary.to_string(index=False)) 