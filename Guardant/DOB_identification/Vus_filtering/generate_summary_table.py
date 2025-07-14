import pandas as pd
import numpy as np

def get_patient_count(df, genes):
    """Counts unique patients with mutations in a given list of genes."""
    return df[df['Gene'].isin(genes)]['Effective Patient ID'].nunique()

def main():
    # --- Configuration ---
    IDC_FILE = 'IDC_genomic_filtered.csv'
    ILC_FILE = 'ILC_genomic_filtered.csv'
    OUTPUT_FILE = 'mutation_summary.csv'

    IDC_TOTAL_PATIENTS = 715
    ILC_TOTAL_PATIENTS = 123
    TOTAL_PATIENTS = IDC_TOTAL_PATIENTS + ILC_TOTAL_PATIENTS

    GENE_GROUPS = {
        'PIK3CA': ['PIK3CA'],
        'ESR1': ['ESR1'],
        'ERBB2': ['ERBB2'],
        'BRCA 1/2': ['BRCA1', 'BRCA2'],
        'PALB2': ['PALB2']
    }
    ALL_TARGET_GENES = [gene for sublist in GENE_GROUPS.values() for gene in sublist]

    # --- Load Data ---
    try:
        idc_df = pd.read_csv(IDC_FILE, low_memory=False)
        ilc_df = pd.read_csv(ILC_FILE, low_memory=False)
        combined_df = pd.concat([idc_df, ilc_df], ignore_index=True)
    except FileNotFoundError as e:
        print(f"Error: Could not find file {e.filename}. Please ensure the filtered files exist.")
        return

    # --- Analysis ---
    results = []

    # 1. Clinically Relevant Mutations (any of the target genes)
    idc_relevant_count = get_patient_count(idc_df, ALL_TARGET_GENES)
    ilc_relevant_count = get_patient_count(ilc_df, ALL_TARGET_GENES)
    total_relevant_count = get_patient_count(combined_df, ALL_TARGET_GENES)
    results.append({
        'Mutation': 'Clinically Relevant Mutations',
        'IDC_Count': idc_relevant_count,
        'ILC_Count': ilc_relevant_count,
        'Total_Count': total_relevant_count
    })

    # 2. Individual Gene/Group Mutations
    for group_name, genes in GENE_GROUPS.items():
        idc_count = get_patient_count(idc_df, genes)
        ilc_count = get_patient_count(ilc_df, genes)
        total_count = get_patient_count(combined_df, genes)
        results.append({
            'Mutation': group_name,
            'IDC_Count': idc_count,
            'ILC_Count': ilc_count,
            'Total_Count': total_count
        })
        
    # --- Formatting ---
    summary_df = pd.DataFrame(results)

    # Calculate and format percentages
    summary_df[f'IDC (n = {IDC_TOTAL_PATIENTS})'] = summary_df.apply(
        lambda row: f"{row['IDC_Count']} ({((row['IDC_Count']/IDC_TOTAL_PATIENTS)*100):.1f})", axis=1
    )
    summary_df[f'ILC (n = {ILC_TOTAL_PATIENTS})'] = summary_df.apply(
        lambda row: f"{row['ILC_Count']} ({((row['ILC_Count']/ILC_TOTAL_PATIENTS)*100):.1f})", axis=1
    )
    summary_df[f'Total (n = {TOTAL_PATIENTS})'] = summary_df.apply(
        lambda row: f"{row['Total_Count']} ({((row['Total_Count']/TOTAL_PATIENTS)*100):.1f})", axis=1
    )

    # Final table formatting
    final_table = summary_df[['Mutation', f'IDC (n = {IDC_TOTAL_PATIENTS})', f'ILC (n = {ILC_TOTAL_PATIENTS})', f'Total (n = {TOTAL_PATIENTS})']]
    final_table = final_table.rename(columns={'Mutation': 'Number of patients (%)'})


    # --- Output ---
    final_table.to_csv(OUTPUT_FILE, index=False)
    print(f"Summary table saved to '{OUTPUT_FILE}'")
    print("\n--- Summary Table ---")
    print(final_table.to_string(index=False))


if __name__ == "__main__":
    main() 