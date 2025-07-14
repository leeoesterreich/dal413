import pandas as pd

def get_mutation_counts(df, genes):
    """Counts the total number of mutations for a given list of genes."""
    return df[df['Gene'].isin(genes)].shape[0]

def main():
    # --- Configuration ---
    IDC_FILE = 'IDC_genomic.csv'
    ILC_FILE = 'ILC_genomic.csv'
    OUTPUT_FILE = 'mutation_counts_summary.csv'

    # Genes to analyze, matching the user's table
    GENES_TO_COUNT = ['PIK3CA', 'ESR1', 'ERBB2', 'BRCA1', 'BRCA2', 'PALB2']

    # --- Load Data ---
    try:
        idc_df = pd.read_csv(IDC_FILE, low_memory=False)
        ilc_df = pd.read_csv(ILC_FILE, low_memory=False)
    except FileNotFoundError as e:
        print(f"Error: Could not find file {e.filename}. Please ensure the original files exist.")
        return

    # --- Analysis ---
    results = []
    for gene in GENES_TO_COUNT:
        idc_count = get_mutation_counts(idc_df, [gene])
        ilc_count = get_mutation_counts(ilc_df, [gene])
        total_count = idc_count + ilc_count
        results.append({
            'Gene': gene,
            'IDC': idc_count,
            'ILC': ilc_count,
            'Total': total_count
        })

    # --- Formatting and Output ---
    summary_df = pd.DataFrame(results)

    summary_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Mutation counts summary saved to '{OUTPUT_FILE}'")
    print("\n--- Mutation Counts ---")
    # Set pandas to display all rows for clean printing
    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
        print(summary_df)


if __name__ == "__main__":
    main() 