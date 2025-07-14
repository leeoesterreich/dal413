import pandas as pd

# Read the genomic data
genomic_file = 'Guardant Project - Genomic Data(UPMC data - Initial Test Only).csv'
genomic_df = pd.read_csv(genomic_file, encoding='latin-1')

# Print column names
print("\nColumns in genomic data:")
print(genomic_df.columns.tolist())

# Print first few rows
print("\nFirst few rows of genomic data:")
print(genomic_df.head().to_string()) 