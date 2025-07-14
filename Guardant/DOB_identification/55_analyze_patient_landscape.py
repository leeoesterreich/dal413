import pandas as pd
import numpy as np
from datetime import datetime

# Read the genomic data files
print("Reading files...")
idc_data = pd.read_csv('IDC_genomic.csv')
ilc_data = pd.read_csv('ILC_genomic.csv')

def analyze_dataset(df, name):
    print("\nAnalyzing {} dataset:".format(name))
    print("-" * 50)
    
    # Basic counts
    total_records = len(df)
    unique_patients = df['Effective Patient ID'].nunique()
    unique_gh_ids = df['GH_ID'].nunique()
    
    print("Total records: {}".format(total_records))
    print("Unique patients: {}".format(unique_patients))
    print("Unique GH IDs: {}".format(unique_gh_ids))
    
    # Test characteristics
    print("\nTest characteristics:")
    print("Test types present:")
    print(df['Test Reported'].value_counts())
    
    # Date range
    df['Sample Received Date'] = pd.to_datetime(df['Sample Received Date'], format='%m/%d/%Y', errors='coerce')
    min_date = df['Sample Received Date'].min()
    max_date = df['Sample Received Date'].max()
    print("\nDate range:")
    print("First sample: {}".format(min_date.strftime('%Y-%m-%d')))
    print("Last sample: {}".format(max_date.strftime('%Y-%m-%d')))
    
    # Age distribution
    print("\nAge distribution:")
    print("Mean age: {:.1f}".format(df['Patient Age'].mean()))
    print("Age range: {}-{}".format(df['Patient Age'].min(), df['Patient Age'].max()))
    
    return unique_patients, unique_gh_ids

# Analyze overlap between datasets
overlapping_ids = set(idc_data['Effective Patient ID']) & set(ilc_data['Effective Patient ID'])
print("\nOverlap analysis:")
print("-" * 50)
print("Patients present in both IDC and ILC: {}".format(len(overlapping_ids)))

# Analyze each dataset
idc_patients, idc_gh = analyze_dataset(idc_data, "IDC")
ilc_patients, ilc_gh = analyze_dataset(ilc_data, "ILC")

# Check GH_ID uniqueness across datasets
all_gh_ids_idc = set(idc_data['GH_ID'])
all_gh_ids_ilc = set(ilc_data['GH_ID'])
overlapping_gh = all_gh_ids_idc & all_gh_ids_ilc

print("\nGH_ID Analysis:")
print("-" * 50)
print("GH_IDs present in both IDC and ILC: {}".format(len(overlapping_gh)))
if len(overlapping_gh) > 0:
    print("\nSample of overlapping GH_IDs:")
    overlap_sample = pd.concat([
        idc_data[idc_data['GH_ID'].isin(list(overlapping_gh)[:5])][['GH_ID', 'Effective Patient ID', 'Patient DOB']],
        ilc_data[ilc_data['GH_ID'].isin(list(overlapping_gh)[:5])][['GH_ID', 'Effective Patient ID', 'Patient DOB']]
    ])
    print(overlap_sample)

# Analyze ICD-O-3 descriptions
print("\nICD-O-3 Description Analysis:")
print("-" * 50)
print("\nIDC descriptions:")
print(idc_data['Histo/Behavior ICD-O-3-Desc'].value_counts().head())
print("\nILC descriptions:")
print(ilc_data['Histo/Behavior ICD-O-3-Desc'].value_counts().head())

# Analyze alterations
def get_alteration_stats(df):
    total_patients = df['Effective Patient ID'].nunique()
    patients_with_alterations = df[df['Alteration'].notna()]['Effective Patient ID'].nunique()
    return total_patients, patients_with_alterations

idc_total, idc_with_alt = get_alteration_stats(idc_data)
ilc_total, ilc_with_alt = get_alteration_stats(ilc_data)

print("\nAlteration Statistics:")
print("-" * 50)
print("IDC patients with alterations: {}/{} ({:.1f}%)".format(
    idc_with_alt, idc_total, (idc_with_alt/float(idc_total))*100))
print("ILC patients with alterations: {}/{} ({:.1f}%)".format(
    ilc_with_alt, ilc_total, (ilc_with_alt/float(ilc_total))*100))

# Save summary to file
with open('patient_landscape_summary.txt', 'w') as f:
    f.write("Patient Cohort Landscape Summary\n")
    f.write("=" * 30 + "\n\n")
    f.write("Total unique patients: {}\n".format(idc_patients + ilc_patients - len(overlapping_ids)))
    f.write("- IDC patients: {}\n".format(idc_patients))
    f.write("- ILC patients: {}\n".format(ilc_patients))
    f.write("- Overlapping patients: {}\n".format(len(overlapping_ids)))
    f.write("\nTotal unique GH_IDs: {}\n".format(len(all_gh_ids_idc | all_gh_ids_ilc)))
    f.write("- IDC GH_IDs: {}\n".format(len(all_gh_ids_idc)))
    f.write("- ILC GH_IDs: {}\n".format(len(all_gh_ids_ilc)))
    f.write("- Overlapping GH_IDs: {}\n".format(len(overlapping_gh))) 