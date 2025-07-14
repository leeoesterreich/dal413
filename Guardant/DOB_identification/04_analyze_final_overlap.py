import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sys

def analyze_final_overlap(genomic_file, ilc_file, idc_file):
    """
    Performs a final analysis using definitive IDs, counting cohorts and
    visualizing the overlap between genomic and histology data.

    Args:
        genomic_file (str): Path to the corrected first-true-alteration data.
        ilc_file (str): Path to the corrected ILC cohort data.
        idc_file (str): Path to the corrected IDC cohort data.
    """
    try:
        df_genomic = pd.read_csv(genomic_file)
        df_ilc = pd.read_csv(ilc_file)
        df_idc = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}", file=sys.stderr)
        sys.exit(1)

    # --- 1. Prepare Patient ID Sets ---
    genomic_ids = set(df_genomic['Effective Patient ID'].unique())
    ilc_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_ids = set(df_idc['Effective Patient ID'].unique())
    
    # --- 2. Report on Histology Cohorts ---
    histology_overlap_ids = ilc_ids.intersection(idc_ids)
    histology_union_ids = ilc_ids.union(idc_ids)
    
    print("--- Histology Cohort Counts (Definitive IDs) ---")
    print(f"Unique patients in ILC cohort: {len(ilc_ids)}")
    print(f"Unique patients in IDC cohort: {len(idc_ids)}")
    print(f"Patients in BOTH ILC and IDC cohorts: {len(histology_overlap_ids)}")
    if len(histology_overlap_ids) > 0:
        print(f"  -> Overlapping ID(s): {sorted(list(histology_overlap_ids))}")
    print(f"Total unique histology patients: {len(histology_union_ids)}")
    
    # --- 3. Calculate Overlaps for Venn Diagram ---
    intersection_ids = genomic_ids.intersection(histology_union_ids)
    genomic_only_ids = genomic_ids - histology_union_ids
    histology_only_ids = histology_union_ids - genomic_ids

    print("\n--- Venn Diagram Data ---")
    print(f"Genomic Data Only (no histology match): {len(genomic_only_ids)}")
    print(f"Histology Data Only (no true alteration): {len(histology_only_ids)}")
    print(f"Intersection (in both data sets): {len(intersection_ids)}")
    
    # --- 4. Calculate ILC/IDC breakdown of the intersection ---
    ilc_in_intersection = intersection_ids.intersection(ilc_ids)
    idc_in_intersection = intersection_ids.intersection(idc_ids)

    print("\n--- Intersection Breakdown ---")
    print(f"ILC patients in intersection: {len(ilc_in_intersection)}")
    print(f"IDC patients in intersection: {len(idc_in_intersection)}")

    # --- 5. Plot Venn Diagram and Pie Chart ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    venn2(subsets=(len(genomic_only_ids), len(histology_only_ids), len(intersection_ids)),
          set_labels=('Genomic Cohort\n(At least one true alteration)', 'Histology Cohort\n(ILC or IDC)'), ax=ax1)
    ax1.set_title('Overlap of Genomic and Histology Cohorts', fontsize=14)
    
    pie_labels = [f'ILC ({len(ilc_in_intersection)})', f'IDC ({len(idc_in_intersection)})']
    pie_sizes = [len(ilc_in_intersection), len(idc_in_intersection)]
    ax2.pie(pie_sizes, labels=pie_labels, autopct='%1.1f%%', startangle=140, textprops={'fontsize': 12})
    ax2.set_title(f'Breakdown of {len(intersection_ids)} Intersection Patients', fontsize=14)
    ax2.axis('equal')

    fig.suptitle('Final Cohort Analysis with Definitive Patient IDs', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plot_file = 'final_cohort_analysis.png'
    plt.savefig(plot_file)
    print(f"\nAnalysis complete. Plot saved to {plot_file}")

if __name__ == "__main__":
    GENOMIC_FILE = 'first_true_test_per_patient.csv'
    ILC_FILE = 'ilc_cohort_with_ids.csv'
    IDC_FILE = 'idc_cohort_with_ids.csv'
    
    analyze_final_overlap(GENOMIC_FILE, ILC_FILE, IDC_FILE) 