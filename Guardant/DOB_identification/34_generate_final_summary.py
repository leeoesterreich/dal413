import pandas as pd
import numpy as np

def generate_final_summary(ilc_path, idc_path, genomic_path, output_path):
    """
    Generates a final summary report of the patient cohorts and their overlap
    with genomic data.

    Args:
        ilc_path (str): Path to the ILC cohort CSV with Effective Patient IDs.
        idc_path (str): Path to the IDC cohort CSV with Effective Patient IDs.
        genomic_path (str): Path to the first positive tests CSV.
        output_path (str): Path to save the summary text file.
    """
    try:
        ilc_df = pd.read_csv(ilc_path)
        idc_df = pd.read_csv(idc_path)
        genomic_df = pd.read_csv(genomic_path)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found: {e.filename}")
        return

    # Get unique patient IDs from each cohort
    ilc_patients = set(ilc_df['Effective Patient ID'].unique())
    idc_patients = set(idc_df['Effective Patient ID'].unique())
    
    # After the correction, there should be no overlap, but we calculate it just in case.
    histology_patients = ilc_patients.union(idc_patients)
    overlap_patients = ilc_patients.intersection(idc_patients)

    # Get unique patient IDs from the genomic data (positive tests)
    genomic_patients = set(genomic_df['Effective Patient ID'].unique())

    # --- Calculations ---
    # Intersection: Patients in histology AND genomic data
    intersection_set = histology_patients.intersection(genomic_patients)

    # Histology Only: Patients in histology but NOT in genomic data
    histology_only_set = histology_patients.difference(genomic_patients)

    # Genomic Only: Patients in genomic data but NOT in histology data
    genomic_only_set = genomic_patients.difference(histology_patients)
    
    # Breakdown of intersection
    intersection_ilc = intersection_set.intersection(ilc_patients)
    intersection_idc = intersection_set.intersection(idc_patients)

    # --- Generate Summary Text ---
    summary = []
    summary.append("Final Cohort Counts:")
    summary.append(f"ILC Cohort: {len(ilc_patients)} unique patients")
    summary.append(f"IDC Cohort: {len(idc_patients)} unique patients")
    if overlap_patients:
        summary.append(f"Overlap: {len(overlap_patients)} patient(s) present in both histology cohorts: {', '.join(list(overlap_patients))}")
    else:
        summary.append("Overlap: 0 patients are present in both histology cohorts.")
    summary.append(f"Total Histology: {len(histology_patients)} unique patients")
    summary.append("\nOverlap with Genomic Data:")
    summary.append(f"Intersection: We found {len(intersection_set)} patients who are in the histology cohorts and have at least one test with a detected alteration.")
    summary.append(f"Histology Only: There are {len(histology_only_set)} patients in the histology cohorts who have no record of a positive test in our genomic data.")
    summary.append(f"Genomic Only: {len(genomic_only_set)} patients. This is a key finding: every single patient with a detected alteration should be accounted for in our histology cohorts.")
    summary.append("\nIntersection Breakdown:")
    summary.append(f"Of the {len(intersection_set)} patients in the intersection:")
    summary.append(f"{len(intersection_ilc)} are from the ILC cohort.")
    summary.append(f"{len(intersection_idc)} are from the IDC cohort.")
    
    note = f"(Note: {len(intersection_ilc)} + {len(intersection_idc)} = {len(intersection_ilc) + len(intersection_idc)}. This should equal the intersection total. If not, it indicates an overlap in the intersection groups.)"
    summary.append(note)

    # --- Write to file ---
    with open(output_path, 'w') as f:
        f.write('\n'.join(summary))
    
    print(f"Final analysis summary saved to '{output_path}'")
    for line in summary:
        print(line)


if __name__ == "__main__":
    ILC_FILE = "ilc_cohort_with_ids.csv"
    IDC_FILE = "idc_cohort_with_ids.csv"
    GENOMIC_FILE = "first_positive_tests_with_all_results.csv"
    OUTPUT_FILE = "analysis_summary.txt"
    
    generate_final_summary(ILC_FILE, IDC_FILE, GENOMIC_FILE, OUTPUT_FILE) 