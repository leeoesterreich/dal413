import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os
import traceback
import re
import sys
import argparse

# Parameters
GENE_SEGMENTS_FILE = 'gene_segments.bed'
QVAL_CUTOFF = 0.01

# Define signatures
SIGNATURES = [
    {
        'short': 'Estrogen_signaling',
        'drivers': ['MAPK3', 'CREBBP', 'ZFHX3']
    },
    {
        'short': 'HER1_C2',
        'drivers': ['HER1', 'EGFR', 'HER2', 'HER4', 'TGFA', 'AREG', 'EGF', 'KRAS', 'PIK3ca', 'AKT1', 'MEK1', 'ERK1']
    },
    {
        'short': 'Basal_signaling',
        'drivers': ['RAD17', 'RAD50', 'PALB2', 'TP53', 'BRCA1']
    }
]

def sanitize_filename(name):
    return re.sub(r'[^\w\-_\. ]', '_', name)

def load_gene_segments():
    """Load gene positions from BED file"""
    print(f"Loading gene segments from {GENE_SEGMENTS_FILE}...")
    if not os.path.exists(GENE_SEGMENTS_FILE):
        print(f"Warning: Gene segments file {GENE_SEGMENTS_FILE} not found. Will use gene indices instead of genomic coordinates.")
        return None
    
    try:
        gene_df = pd.read_csv(GENE_SEGMENTS_FILE, sep='\t', header=None, 
                              names=['chrom', 'start', 'end', 'gene', 'score', 'strand'])
        print(f"Loaded {len(gene_df)} gene segments")
        return gene_df
    except Exception as e:
        print(f"Error loading gene segments: {str(e)}")
        return None

def plot_association_results(signature_short, driver_genes, input_csv=None, output_png=None, 
                             max_genes=None, dpi=300, adjust_labels=True, q_cutoff=None, 
                             plot_type='bar'):
    """
    Plot association results for a signature from a CSV file.
    
    Parameters:
    - signature_short: Short name of the signature
    - driver_genes: List of driver genes for the signature
    - input_csv: Path to input CSV file (default: signature_short + "_association_results.csv")
    - output_png: Path to output PNG file (default: signature_short + "_association_landscape_genomewide.png")
    - max_genes: Maximum number of genes to plot (for memory optimization)
    - dpi: Resolution of the output image
    - adjust_labels: Whether to adjust driver gene labels
    - q_cutoff: Custom q-value cutoff (default: use QVAL_CUTOFF)
    - plot_type: Type of plot to use ('bar', 'vlines', or 'scatter')
    """
    if q_cutoff is None:
        q_cutoff = QVAL_CUTOFF
        
    if input_csv is None:
        input_csv = f"{sanitize_filename(signature_short)}_association_results.csv"
    
    if output_png is None:
        output_png = f"{sanitize_filename(signature_short)}_association_landscape_genomewide.png"
    
    print(f"\n=== Plotting results for {signature_short} ===")
    print(f"Loading data from {input_csv}...")
    
    # Check if input file exists
    if not os.path.exists(input_csv):
        print(f"Error: Input file {input_csv} does not exist!")
        return False
    
    try:
        # Load association results
        res_df = pd.read_csv(input_csv)
        
        # Only plot genes where any q-value is < q_cutoff
        sig_mask = (
            (res_df['spearman_pos_p_adj'] < q_cutoff) |
            (res_df['spearman_neg_p_adj'] < q_cutoff) |
            (res_df['gain_p_adj'] < q_cutoff) |
            (res_df['loss_p_adj'] < q_cutoff)
        )
        sig_df = res_df[sig_mask].copy()
        
        # Sort genes to ensure consistent order
        sig_df = sig_df.sort_values('gene')
        sig_df = sig_df.reset_index(drop=True)
        
        # Limit the number of genes if max_genes is specified
        if max_genes is not None and len(sig_df) > max_genes:
            print(f"Limiting to {max_genes} genes for memory optimization (from {len(sig_df)} significant genes)")
            # Sort by the most significant q-values
            min_q_values = sig_df[['spearman_pos_p_adj', 'spearman_neg_p_adj', 'gain_p_adj', 'loss_p_adj']].min(axis=1)
            sig_df['min_q'] = min_q_values
            sig_df = sig_df.sort_values('min_q').head(max_genes)
            sig_df = sig_df.drop('min_q', axis=1)
            sig_df = sig_df.reset_index(drop=True)
        
        if len(sig_df) == 0:
            print(f"No significant genes found for {signature_short}, skipping plot.")
            return False
            
        print(f"Plotting results for {signature_short} ({len(sig_df)} significant genes)...")
        
        # Create a new figure
        plt.figure(figsize=(18, 5))
        
        # Create masks for different categories
        pos_mask = (sig_df['spearman_pos_p_adj'] < q_cutoff) & (sig_df['gain_p_adj'] < q_cutoff)
        neg_mask = (sig_df['spearman_neg_p_adj'] < q_cutoff) & (sig_df['loss_p_adj'] < q_cutoff)
        gain_mask = (sig_df['gain_p_adj'] < q_cutoff)
        loss_mask = (sig_df['loss_p_adj'] < q_cutoff)
        
        print(f"  Number of positive correlations: {sum(pos_mask)}")
        print(f"  Number of negative correlations: {sum(neg_mask)}")
        print(f"  Number of gains: {sum(gain_mask)}")
        print(f"  Number of losses: {sum(loss_mask)}")
        
        # Create x positions for plotting
        x_pos = np.arange(len(sig_df))
        
        # Plot based on selected type
        if plot_type == 'bar':
            # Use bar plots (original approach but optimized)
            print("  Using bar plot")
            bar_width = 0.8  # Width of the bars
            
            # Only plot genes with significant values
            if sum(pos_mask) > 0:
                print("  Plotting positive correlations...")
                pos_indices = np.where(pos_mask)[0]
                pos_values = -np.log10(sig_df.loc[pos_mask, 'spearman_pos_p_adj'].values)
                plt.bar(x_pos[pos_indices], pos_values, width=bar_width, color='red', 
                        alpha=0.7, label='Positive correlation')
            
            if sum(neg_mask) > 0:
                print("  Plotting negative correlations...")
                neg_indices = np.where(neg_mask)[0]
                neg_values = np.log10(sig_df.loc[neg_mask, 'spearman_neg_p_adj'].values)
                plt.bar(x_pos[neg_indices], neg_values, width=bar_width, color='darkblue', 
                        alpha=0.7, label='Negative correlation')
            
            if sum(gain_mask) > 0:
                print("  Plotting gains...")
                gain_indices = np.where(gain_mask)[0]
                gain_values = -np.log10(sig_df.loc[gain_mask, 'gain_p_adj'].values)
                plt.bar(x_pos[gain_indices], gain_values, width=bar_width, color='orange', 
                        alpha=0.5, label='Gain')
            
            if sum(loss_mask) > 0:
                print("  Plotting losses...")
                loss_indices = np.where(loss_mask)[0]
                loss_values = np.log10(sig_df.loc[loss_mask, 'loss_p_adj'].values)
                plt.bar(x_pos[loss_indices], loss_values, width=bar_width, color='lightblue', 
                        alpha=0.5, label='Loss')
                
        elif plot_type == 'vlines':
            # Use vertical lines (faster alternative)
            print("  Using vertical lines plot")
            if sum(pos_mask) > 0:
                print("  Plotting positive correlations...")
                pos_indices = np.where(pos_mask)[0]
                pos_x = x_pos[pos_indices]
                pos_y = -np.log10(sig_df.loc[pos_mask, 'spearman_pos_p_adj'].values)
                plt.vlines(pos_x, np.zeros_like(pos_y), pos_y, color='red', alpha=0.7, label='Positive correlation')
            
            if sum(neg_mask) > 0:
                print("  Plotting negative correlations...")
                neg_indices = np.where(neg_mask)[0]
                neg_x = x_pos[neg_indices]
                neg_y = np.log10(sig_df.loc[neg_mask, 'spearman_neg_p_adj'].values)
                plt.vlines(neg_x, np.zeros_like(neg_y), neg_y, color='darkblue', alpha=0.7, label='Negative correlation')
            
            if sum(gain_mask) > 0:
                print("  Plotting gains...")
                gain_indices = np.where(gain_mask)[0]
                gain_x = x_pos[gain_indices]
                gain_y = -np.log10(sig_df.loc[gain_mask, 'gain_p_adj'].values)
                plt.vlines(gain_x, np.zeros_like(gain_y), gain_y, color='orange', alpha=0.5, label='Gain')
            
            if sum(loss_mask) > 0:
                print("  Plotting losses...")
                loss_indices = np.where(loss_mask)[0]
                loss_x = x_pos[loss_indices]
                loss_y = np.log10(sig_df.loc[loss_mask, 'loss_p_adj'].values)
                plt.vlines(loss_x, np.zeros_like(loss_y), loss_y, color='lightblue', alpha=0.5, label='Loss')
                
        else:  # scatter
            # Use scatter plot
            print("  Using scatter plot")
            if sum(pos_mask) > 0:
                print("  Plotting positive correlations...")
                pos_indices = np.where(pos_mask)[0]
                plt.scatter(x_pos[pos_indices], -np.log10(sig_df.loc[pos_mask, 'spearman_pos_p_adj']), 
                           color='red', s=30, alpha=0.7, label='Positive correlation')
            
            if sum(neg_mask) > 0:
                print("  Plotting negative correlations...")
                neg_indices = np.where(neg_mask)[0]
                plt.scatter(x_pos[neg_indices], np.log10(sig_df.loc[neg_mask, 'spearman_neg_p_adj']), 
                           color='darkblue', s=30, alpha=0.7, label='Negative correlation')
            
            if sum(gain_mask) > 0:
                print("  Plotting gains...")
                gain_indices = np.where(gain_mask)[0]
                plt.scatter(x_pos[gain_indices], -np.log10(sig_df.loc[gain_mask, 'gain_p_adj']), 
                           color='orange', s=20, alpha=0.5, label='Gain')
            
            if sum(loss_mask) > 0:
                print("  Plotting losses...")
                loss_indices = np.where(loss_mask)[0]
                plt.scatter(x_pos[loss_indices], np.log10(sig_df.loc[loss_mask, 'loss_p_adj']), 
                           color='lightblue', s=20, alpha=0.5, label='Loss')
        
        # Annotate driver genes with repelled labels
        print("  Annotating driver genes...")
        texts = []
        for gene in driver_genes:
            if gene in sig_df['gene'].values:
                gene_idx = sig_df[sig_df['gene'] == gene].index[0]
                texts.append(plt.text(gene_idx, 0, gene, fontsize=16, fontweight='bold', ha='center', va='bottom'))
        
        if adjust_labels and texts:
            try:
                print("  Adjusting text positions...")
                adjust_text(texts, arrowprops=dict(arrowstyle='-|>', color='black', shrinkA=5), 
                            force_text=1.5, force_points=1.5)
            except Exception as e:
                print(f"  Warning: Error adjusting text positions: {str(e)}")
                print("  Continuing with unadjusted text...")
        
        # Add significance threshold lines
        plt.axhline(-np.log10(q_cutoff), color='k', linestyle='dashed')
        plt.axhline(np.log10(q_cutoff), color='k', linestyle='dashed')
        
        # Labels and formatting
        plt.xlabel('Genes (ordered)', fontsize=18)
        plt.ylabel('-log10(q)', fontsize=18)
        plt.title(f'{signature_short} association landscape (significant genes)', fontsize=20)
        plt.xticks([])  # Hide x-ticks since there are too many genes
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14, loc='upper right')
        plt.tight_layout()
        
        # Save figure
        print(f"  Saving plot to {output_png} with DPI={dpi}...")
        plt.savefig(output_png, dpi=dpi, bbox_inches='tight')
        plt.close()
        print(f"Plot saved as {output_png}")
        return True
        
    except Exception as e:
        print(f"Error during plotting for {signature_short} signature:")
        print(traceback.format_exc())
        return False

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Plot association results from CSV files')
    parser.add_argument('signatures', nargs='*', help='Names of signatures to plot (default: all signatures)')
    parser.add_argument('--max-genes', type=int, help='Maximum number of genes to plot (for memory optimization)')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for output images (default: 300)')
    parser.add_argument('--no-adjust-labels', action='store_true', help='Disable label adjustment for driver genes')
    parser.add_argument('--q-cutoff', type=float, default=QVAL_CUTOFF, help=f'Q-value cutoff for significance (default: {QVAL_CUTOFF})')
    parser.add_argument('--plot-type', choices=['bar', 'vlines', 'scatter'], default='bar', 
                        help='Type of plot (bar, vlines, or scatter, default: bar)')
    return parser.parse_args()

def main():
    """Plot association results for all signatures or specified ones."""
    args = parse_arguments()
    
    # Check if signature names were provided as command-line arguments
    if args.signatures:
        sig_names = args.signatures
        print(f"Plotting results for specified signatures: {', '.join(sig_names)}")
        # Filter signatures based on command-line arguments
        signatures_to_plot = [sig for sig in SIGNATURES if sig['short'] in sig_names]
        if not signatures_to_plot:
            print("Warning: None of the specified signatures were found!")
            print(f"Available signatures: {', '.join(sig['short'] for sig in SIGNATURES)}")
            return
    else:
        print("Plotting results for all signatures")
        signatures_to_plot = SIGNATURES
    
    # Plot results for each signature
    for sig in signatures_to_plot:
        plot_association_results(
            sig['short'], 
            sig['drivers'],
            max_genes=args.max_genes,
            dpi=args.dpi,
            adjust_labels=not args.no_adjust_labels,
            q_cutoff=args.q_cutoff,
            plot_type=args.plot_type
        )
    
    print("\nPlotting complete!")

if __name__ == "__main__":
    main() 