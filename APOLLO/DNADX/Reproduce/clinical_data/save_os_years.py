import pandas as pd

def load_clinical_data():
    """Load clinical data and convert OS to years"""
    clinical_data = pd.read_csv('data_clinical_patient.txt', sep='\t', skiprows=4)
    # Convert overall survival in months to years
    clinical_data['OS_YEARS'] = clinical_data['OS_MONTHS'] / 12
    # Set the index to PATIENT_ID to ensure alignment with scores
    clinical_data.set_index('PATIENT_ID', inplace=True)
    return clinical_data

def main():
    # Load clinical data
    clinical_data = load_clinical_data()
    
    # Save OS_YEARS and PATIENT_ID to a CSV file
    clinical_data[['OS_YEARS']].to_csv('os_years_and_patient_id.csv')
    print("OS_YEARS and PATIENT_ID have been saved to 'os_years_and_patient_id.csv'.")

if __name__ == "__main__":
    main() 