import pandas as pd
import os

# Define the base directory
BASE_DIR = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/DOB_identification/Vus_filtering'

# Read the data file
try:
    idc_genomic = pd.read_csv(os.path.join(BASE_DIR, 'IDC_genomic.csv'), encoding='latin-1', low_memory=False)
    print("File loaded successfully.")
    print("Columns are:")
    print(idc_genomic.columns)
except FileNotFoundError as e:
    print(f"Error: An input file was not found. {e}")
except Exception as e:
    print(f"An error occurred: {e}") 