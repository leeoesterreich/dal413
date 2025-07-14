import pandas as pd
import numpy as np

def analyze_serial_tests(input_file, output_file, report_file):
    """
    Analyzes serial genomic data to verify test completeness for each patient.

    This script reads a CSV file with serial genomic data, assigns a unique
    Patient ID to each patient based on their Date of Birth (DOB), and checks
    if the 'Serial Test' numbers for each patient form a complete sequence
    from 1 to their 'Total Test' count.

    Args:
        input_file (str): The path to the input CSV file.
        output_file (str): The path to save the output CSV file with Patient IDs.
        report_file (str): The path to save the analysis report.
    """
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file {input_file} was not found.")
        return

    # Create a unique Patient ID for each unique Patient DOB
    df['Patient ID'] = pd.factorize(df['Patient DOB'])[0] + 1
    df['Patient ID'] = 'P' + df['Patient ID'].astype(str)

    # Move 'Patient ID' to the first column
    cols = ['Patient ID'] + [col for col in df if col != 'Patient ID']
    df = df[cols]

    # Group by the new Patient ID
    grouped = df.groupby('Patient ID')

    report_content = []
    report_content.append("Analyzing patient test series...")
    print("Analyzing patient test series...")
    
    discrepancy_patients = []

    for patient_id, group in grouped:
        total_tests = group['Total Test'].iloc[0]
        serial_tests = set(group['Serial Test'])
        expected_tests = set(range(1, total_tests + 1))

        missing_tests = sorted(list(expected_tests - serial_tests))
        extra_tests = sorted(list(serial_tests - expected_tests))

        if missing_tests or extra_tests:
            discrepancy_patients.append({
                'Patient ID': patient_id,
                'Patient DOB': group['Patient DOB'].iloc[0],
                'Total Tests': total_tests,
                'Actual Tests': sorted(list(serial_tests)),
                'Missing Tests': missing_tests,
                'Overloaded Tests': extra_tests
            })

    if discrepancy_patients:
        msg = "\nPatients with incomplete or overloaded test series found:"
        report_content.append(msg)
        print(msg)
        for patient in discrepancy_patients:
            msg = (f"  Patient ID: {patient['Patient ID']} (DOB: {patient['Patient DOB']})\n"
                   f"    - Expected {patient['Total Tests']} tests.\n"
                   f"    - Found tests: {patient['Actual Tests']}")
            if patient['Missing Tests']:
                msg += f"\n    - Missing tests: {patient['Missing Tests']}"
            if patient['Overloaded Tests']:
                msg += f"\n    - Overloaded tests: {patient['Overloaded Tests']}"
            report_content.append(msg)
            print(msg)
        
        summary_msg = f"\nTotal patients with discrepancies: {len(discrepancy_patients)}"
        report_content.append(summary_msg)
        print(summary_msg)

    else:
        msg = "\nAll patients have a complete series of tests."
        report_content.append(msg)
        print(msg)

    # Save the dataframe with the new Patient ID column
    df.to_csv(output_file, index=False)
    msg = f"\nData with Patient IDs saved to {output_file}"
    report_content.append(msg)
    print(msg)

    # Write report to file
    with open(report_file, 'w') as f:
        f.write('\n'.join(report_content))
    print(f"Analysis report saved to {report_file}")


if __name__ == "__main__":
    INPUT_CSV = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/Guardant/Guardant Project - Serial Genomic Data.csv'
    OUTPUT_CSV = 'serial_genomic_with_patient_id.csv'
    REPORT_TXT = 'serial_test_analysis_report.txt'
    analyze_serial_tests(INPUT_CSV, OUTPUT_CSV, REPORT_TXT) 