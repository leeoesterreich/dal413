import pandas as pd
import re

def clean_sample_name(x):
    return re.sub(r'\s*\([^)]*\)', '', str(x)).strip()

df_conc = pd.read_csv('Concentration_summary_IDC.csv')
df_pat = pd.read_csv('Patient_Information_IDC.csv')
df_conc['Sample_clean'] = df_conc['Sample_name'].apply(clean_sample_name)

count_yes = 0
for idx, row in df_pat.iterrows():
    sample_names = [clean_sample_name(row[col]) for col in df_pat.columns if col.startswith('Draw') and pd.notnull(row[col]) and str(row[col]).strip()]
    overlap = [s for s in sample_names if s in set(df_conc['Sample_clean'])]
    if row['RA status'].strip() == 'Yes' and overlap:
        count_yes += 1
print(f'Rows with RA status Yes and at least one overlapping sample: {count_yes}') 