#!/usr/bin/env python

from __future__ import print_function
import sys
import requests
import json

def get_ncbi_ids(gene_symbols):
    """Convert gene symbols to NCBI gene IDs using MyGene.info API"""
    gene_map = {}
    
    # Query each gene individually for exact matches
    for symbol in gene_symbols:
        # Query MyGene.info API with exact symbol match
        url = 'https://mygene.info/v3/query?q=symbol:"{}"&species=human&fields=entrezgene'.format(symbol)
        response = requests.get(url)
        data = response.json()
        
        # Check for exact matches
        for hit in data.get('hits', []):
            if 'entrezgene' in hit:
                gene_map[symbol] = str(hit['entrezgene'])
                break
    
    return gene_map

def main():
    # Read input file
    input_file = sys.argv[1]
    output_file = input_file.rsplit('.', 1)[0] + '_ncbi.txt'
    
    # Read gene symbols
    with open(input_file, 'r') as f:
        gene_symbols = [line.strip() for line in f if line.strip()]
    
    print("Converting {} gene symbols to NCBI IDs...".format(len(gene_symbols)))
    
    # Get NCBI IDs
    gene_map = get_ncbi_ids(gene_symbols)
    
    # Write output file
    with open(output_file, 'w') as f:
        for symbol in gene_symbols:
            if symbol in gene_map:
                f.write("{}\t{}\n".format(symbol, gene_map[symbol]))
            else:
                print("Warning: No NCBI ID found for {}".format(symbol))
    
    print("\nConverted IDs written to: {}".format(output_file))
    print("Found NCBI IDs for {} out of {} genes".format(len(gene_map), len(gene_symbols)))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python convert_gene_ids.py <input_file>")
        sys.exit(1)
    main() 