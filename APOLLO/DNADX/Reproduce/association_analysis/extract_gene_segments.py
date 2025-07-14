import pandas as pd

# Read refFlat.txt (tab-separated, no header)
refFlat = pd.read_csv('refFlat.txt', sep='\t', header=None)
refFlat.columns = [
    'geneName', 'transcriptName', 'chrom', 'strand',
    'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
    'exonCount', 'exonStarts', 'exonEnds'
]

# Group by geneName and chrom, get min(txStart) and max(txEnd) for each gene
# (Some genes may be on multiple chromosomes, rare but possible)
gene_coords = refFlat.groupby(['geneName', 'chrom']).agg({
    'txStart': 'min',
    'txEnd': 'max'
}).reset_index()

# Prepare BED-like output: chrom, txStart, txEnd, geneName
gene_coords = gene_coords[['chrom', 'txStart', 'txEnd', 'geneName']]

# Save to BED file
gene_coords.to_csv('gene_segments.bed', sep='\t', header=False, index=False)

print('BED-like gene segment file saved as gene_segments.bed') 