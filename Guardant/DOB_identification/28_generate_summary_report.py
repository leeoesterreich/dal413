import pandas as pd
import sys

def generate_summary_report(serial_file, ilc_file, idc_file, report_output_file):
    """
    Generates a summary report verifying cohort presence in genomic data and
    counting patients with no detected alterations.

    Args:
        serial_file (str): Path to the main serial genomic data CSV.
        ilc_file (str): Path to the ILC cohort CSV.
        idc_file (str): Path to the IDC cohort CSV.
        report_output_file (str): Path to save the summary report TXT file.
    """
    try:
        df_serial = pd.read_csv(serial_file)
        df_ilc = pd.read_csv(ilc_file)
        df_idc = pd.read_csv(idc_file)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found - {e}")
        sys.exit(1)

    report_content = []
    report_content.append("Histology and Genomic Data Cross-Verification Report")
    report_content.append("="*55)
    report_content.append("\n")

    # --- Part 1: Verification of Patient Cohorts ---
    report_content.append("Part 1: Cohort Presence Verification")
    report_content.append("---------------------------------------")
    
    serial_patient_ids = set(df_serial['Effective Patient ID'].unique())
    ilc_patient_ids = set(df_ilc['Effective Patient ID'].unique())
    idc_patient_ids = set(df_idc['Effective Patient ID'].unique())

    ilc_missing = ilc_patient_ids - serial_patient_ids
    if not ilc_missing:
        report_content.append("✅ All patients from the ILC cohort are present in the serial test data.")
    else:
        report_content.append(f"❌ Found {len(ilc_missing)} patients from the ILC cohort NOT present in the serial test data.")

    idc_missing = idc_patient_ids - serial_patient_ids
    if not idc_missing:
        report_content.append("✅ All patients from the IDC cohort are present in the serial test data.")
    else:
        report_content.append(f"❌ Found {len(idc_missing)} patients from the IDC cohort NOT present in the serial test data.")
    
    report_content.append("\n")

    # --- Part 2: Analysis of Patients with No Alterations ---
    report_content.append("Part 2: Patients Without Detected Alterations")
    report_content.append("-------------------------------------------------")
    
    patients_with_true = set(df_serial[df_serial['Aleration Detected?'] == True]['Effective Patient ID'].unique())
    
    ilc_all_negative = ilc_patient_ids - patients_with_true
    idc_all_negative = idc_patient_ids - patients_with_true

    report_content.append(f"Number of ILC patients with no alterations detected: {len(ilc_all_negative)}")
    report_content.append(f"Number of IDC patients with no alterations detected: {len(idc_all_negative)}")

    # --- Write Report to File ---
    try:
        with open(report_output_file, 'w') as f:
            f.write('\n'.join(report_content))
        print(f"Successfully generated summary report.")
        print(f"Report saved to {report_output_file}")
    except IOError as e:
        print(f"Error writing report file: {e}")


if __name__ == "__main__":
    SERIAL_GENOMIC_FILE = 'serial_genomic_with_patient_id.csv'
    ILC_COHORT_FILE = 'ilc_cohort_with_patient_id.csv'
    IDC_COHORT_FILE = 'idc_cohort_with_patient_id.csv'
    REPORT_FILE = 'histology_genomic_summary_report.txt'
    
    generate_summary_report(SERIAL_GENOMIC_FILE, ILC_COHORT_FILE, IDC_COHORT_FILE, REPORT_FILE) 