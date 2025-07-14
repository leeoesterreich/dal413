#!/usr/bin/env python3
"""Simple script to annotate CNV segments with gene names from a BED file."""
import argparse
import sys
from collections import defaultdict

def load_bed_genes(bed_file):
    """Load gene locations from BED file into a chromosome-based dictionary."""
    genes = defaultdict(list)
    try:
        with open(bed_file) as f:
            for line in f:
                if line.startswith('#'):  # Skip comment lines
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:  # Ensure we have chr, start, end, and gene name
                    chrom, start, end, gene = fields[0:4]
                    # Store chromosome without 'chr' prefix for consistency
                    chrom = chrom.replace('chr', '')
                    genes[chrom].append((int(start), int(end), gene))
    except Exception as e:
        print(f"Error reading BED file: {e}", file=sys.stderr)
        sys.exit(1)
    return genes

def find_overlapping_genes(chrom, start, end, genes):
    """Find all genes that overlap with the given segment."""
    # Remove 'chr' prefix if present for consistency
    chrom = chrom.replace('chr', '')
    if chrom not in genes:
        return []
    
    overlapping = []
    for gene_start, gene_end, gene_name in genes[chrom]:
        if start <= gene_end and end >= gene_start:
            overlapping.append(gene_name)
    return overlapping

def process_cnv_file(cnv_file, genes, output_file):
    """Process CNV file and add gene annotations."""
    try:
        with open(cnv_file) as f_in, open(output_file, 'w') as f_out:
            # Read header and add genes column
            header = f_in.readline().strip()
            f_out.write(f"{header}\tgenes\n")
            
            # Process each line
            for line in f_in:
                fields = line.strip().split('\t')
                if len(fields) < 3:  # Skip malformed lines
                    continue
                
                try:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    # Find overlapping genes
                    overlapping = find_overlapping_genes(chrom, start, end, genes)
                    gene_annotation = ';'.join(overlapping) if overlapping else '-'
                    
                    # Write the original line plus gene annotation
                    f_out.write(f"{line.strip()}\t{gene_annotation}\n")
                except ValueError as e:
                    print(f"Warning: Skipping malformed line: {line.strip()}", file=sys.stderr)
                    continue
                
    except Exception as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bed_file', help="BED file with gene annotations")
    parser.add_argument('cnv_file', help="CNV segment file (.cna.seg)")
    parser.add_argument('-o', '--output', help="Output file name", required=True)
    args = parser.parse_args()

    print(f"Loading gene annotations from {args.bed_file}...", file=sys.stderr)
    genes = load_bed_genes(args.bed_file)
    
    print(f"Processing {args.cnv_file}...", file=sys.stderr)
    process_cnv_file(args.cnv_file, genes, args.output)
    print(f"Done. Output written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
