import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sys

def analyze_cohort_overlap(genomic_file, ilc_file, idc_file):
    """
    Analyzes and visualizes the overlap between the genomic cohort and a
    filtered histology cohort (only including patients with true alterations).

    Args:
        genomic_file (str): Path to the first-true-alteration genomic data.
        ilc_file (str): Path to the ILC cohort data.
        idc_file (str): Path to the IDC cohort data.
    """
    try:
        df_genomic = pd.read_csv(genomic_file)
        df_ilc = pd.read_csv(ilc_file, low_memory=False)
        df_idc = pd.read_csv(idc_file, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    # --- 1. Prepare Patient ID Sets ---
    genomic_ids = set(df_genomic['Effective Patient ID'].unique())
    ilc_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_ids = set(df_idc['Effective Patient ID'].unique())
    
    # Define the histology cohort as patients who are in ILC/IDC AND have a true alteration
    all_histology_ids = ilc_ids.union(idc_ids)
    histology_with_true_ids = all_histology_ids.intersection(genomic_ids)

    # --- 2. Calculate Overlaps for Venn Diagram ---
    intersection_ids = genomic_ids.intersection(histology_with_true_ids)
    genomic_only_ids = genomic_ids - histology_with_true_ids
    histology_only_ids = histology_with_true_ids - genomic_ids # This will be 0 by definition

    # --- 3. Calculate ILC/IDC breakdown of the intersection ---
    ilc_in_intersection = intersection_ids.intersection(ilc_ids)
    idc_in_intersection = intersection_ids.intersection(idc_ids)

    print("--- Cohort Overlap Analysis (Histology Filtered) ---")
    print(f"Total patients in genomic data (first true alteration): {len(genomic_ids)}")
    print(f"Total patients in histology data with a true alteration: {len(histology_with_true_ids)}")
    print(f"  - ILC patients in intersection: {len(ilc_in_intersection)}")
    print(f"  - IDC patients in intersection: {len(idc_in_intersection)}")
    print(f"Patients in genomic data only (no histology match): {len(genomic_only_ids)}")
    
    # --- 4. Plot Venn Diagram and Pie Chart ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    # Venn Diagram
    venn2(subsets=(len(genomic_only_ids), len(histology_only_ids), len(intersection_ids)),
          set_labels=('Genomic Cohort', 'Histology Cohort\n(with True Alteration)'), ax=ax1)
    ax1.set_title('Overlap of Genomic and Filtered Histology Cohorts', fontsize=14)
    
    # Pie Chart for the intersection
    pie_labels = [f'ILC ({len(ilc_in_intersection)})', f'IDC ({len(idc_in_intersection)})']
    pie_sizes = [len(ilc_in_intersection), len(idc_in_intersection)]
    ax2.pie(pie_sizes, labels=pie_labels, autopct='%1.1f%%', startangle=140, textprops={'fontsize': 12})
    ax2.set_title(f'Breakdown of {len(intersection_ids)} Intersection Patients', fontsize=14)
    ax2.axis('equal')

    fig.suptitle('Cohort Analysis: Overlap with Filtered Histology and Intersection Breakdown', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plot_file = 'cohort_analysis_filtered_venn_pie.png'
    plt.savefig(plot_file)
    print(f"\nCombined plot saved to {plot_file}")

    # --- 5. Filter and Save CSVs ---
    # Genomic only
    df_genomic_only = df_genomic[df_genomic['Effective Patient ID'].isin(genomic_only_ids)]
    df_genomic_only.to_csv('genomic_only_patients.csv', index=False)
    print(f"Saved {len(df_genomic_only)} records to genomic_only_patients.csv")

    # Intersection
    df_intersection = df_genomic[df_genomic['Effective Patient ID'].isin(intersection_ids)]
    df_intersection.to_csv('intersection_patients.csv', index=False)
    print(f"Saved {len(df_intersection)} records to intersection_patients.csv")

if __name__ == "__main__":
    GENOMIC_FILE = 'first_true_test_per_patient.csv'
    ILC_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_FILE = 'idc_cohort_with_patient_id.csv'
    
    analyze_cohort_overlap(GENOMIC_FILE, ILC_FILE, IDC_FILE) 