import pandas as pd
import numpy as np
from scipy import stats

def print_file_info(file_path, file_type):
    print(f"Inspecting {file_type} file: {file_path}")
    try:
        data = pd.read_pickle(file_path)
        print(f"  Successfully loaded. Shape: {data.shape}")
        print(f"  Data type: {type(data)}")
        if isinstance(data, pd.DataFrame):
            print(f"  DataFrame head:\n{data.head()}")
            print(f"  DataFrame describe (overall stats for numeric columns):\n{data.describe(include=np.number)}")
            if data.shape[0] > 0 and data.shape[1] > 0:
                 print(f"  Value of first element: {data.iloc[0,0]}")
        elif isinstance(data, pd.Series):
            print(f"  Series head:\n{data.head()}")
            print(f"  Series describe:\n{data.describe()}")
            if data.shape[0] > 0:
                print(f"  Value of first element: {data.iloc[0]}")
        else:
            print(f"  Data is of an unexpected type: {type(data)}")
        
    except FileNotFoundError:
        print(f"  Error: File not found at {file_path}")
    except Exception as e:
        print(f"  Error loading or inspecting file: {e}")
    print("-"*50)

cna_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/training_data/cna_segment_score_mean_no_norm.pkl'
rna_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/scaling/training_data/rna_signature_score_median_no_norm.pkl'

print_file_info(cna_file, 'CNA score')
print_file_info(rna_file, 'RNA score') 