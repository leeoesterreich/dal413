import pandas as pd
import numpy as np

def filter_first_positive_tests(data_path, master_list_path, output_path):
    """
    For each patient, finds the date of their first positive test and extracts all
    records for that patient on that specific date.

    Args:
        data_path (str): Path to the input genomic data CSV.
        master_list_path (str): Path to the patient master list CSV.
        output_path (str): Path to save the filtered output CSV.
    """
    # Load the genomic data with 'utf-8-sig' encoding to handle BOM
    try:
        df = pd.read_csv(data_path, encoding='utf-8-sig')
        master_list = pd.read_csv(master_list_path)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found.")
        print(e)
        return

    # Merge with master list to get Effective Patient ID
    df = pd.merge(df, master_list[['GH_ID', 'Effective Patient ID']], on='GH_ID', how='left')

    # Ensure date columns are in datetime format
    df['Sample Received Date'] = pd.to_datetime(df['Sample Received Date'])
    
    # Correct potential typo in column name
    if 'Aleration Detected?' in df.columns:
        df.rename(columns={'Aleration Detected?': 'Alteration Detected?'}, inplace=True)

    # Filter for positive tests
    positive_tests = df[df['Alteration Detected?'] == True].copy()

    if positive_tests.empty:
        print("No positive tests found in the dataset.")
        return

    # Find the first positive test date for each patient
    first_positive_dates = positive_tests.loc[positive_tests.groupby('Effective Patient ID')['Sample Received Date'].idxmin()]
    first_positive_dates = first_positive_dates[['Effective Patient ID', 'Sample Received Date']].rename(columns={'Sample Received Date': 'First Positive Date'})
    
    # Merge this back to the original dataframe to get all records on that date
    merged_df = pd.merge(df, first_positive_dates, on='Effective Patient ID')
    
    # Filter for rows where the test date matches the first positive date
    result_df = merged_df[merged_df['Sample Received Date'] == merged_df['First Positive Date']].copy()

    # Clean up the dataframe
    result_df.drop(columns=['First Positive Date'], inplace=True)
    
    # Save the result
    result_df.to_csv(output_path, index=False)
    print(f"Filtered data saved to '{output_path}'")
    print(f"Found first positive test records for {result_df['Effective Patient ID'].nunique()} unique patients.")
    print(f"Total records saved: {len(result_df)}")


if __name__ == "__main__":
    INPUT_FILE = "Guardant Project - Genomic Data_DOB_identification.csv"
    MASTER_LIST_FILE = "patient_master_list.csv"
    OUTPUT_FILE = "first_positive_tests_with_all_results.csv"
    filter_first_positive_tests(INPUT_FILE, MASTER_LIST_FILE, OUTPUT_FILE) 