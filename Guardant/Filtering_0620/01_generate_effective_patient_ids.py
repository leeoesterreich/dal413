import pandas as pd
import sys

def generate_definitive_patient_ids(ilc_file, idc_file, master_output, ilc_output, idc_output):
    """
    Generates a definitive 'Effective Patient ID' that accounts for multiple
    patients sharing the same DOB, by appending a letter suffix based on GH_ID.
    """
    try:
        print("Step 1: Reading cohort files...")
        df_ilc = pd.read_csv(ilc_file, low_memory=False)
        df_idc = pd.read_csv(idc_file, low_memory=False)
        print("...done.")

        print("Step 2: Combining data and creating Base PIDs...")
        df_ilc.columns = df_ilc.columns.str.strip()
        df_idc.columns = df_idc.columns.str.strip()
        df_combined = pd.concat([df_ilc, df_idc], ignore_index=True)

        df_combined['Patient DOB'] = pd.to_datetime(df_combined['Date of Birth'], format='%Y%m%d', errors='coerce').dt.strftime('%-m/%-d/%Y')
        df_combined.dropna(subset=['Patient DOB'], inplace=True)

        unique_dobs = df_combined['Patient DOB'].unique()
        dob_to_int_id = {dob: i + 1 for i, dob in enumerate(unique_dobs)}
        df_combined['Base PID'] = 'P' + df_combined['Patient DOB'].map(dob_to_int_id).astype(str)
        print("...done.")

        print("Step 3: Generating suffix-based PID map...")
        grouped = df_combined.groupby('Patient DOB')
        pid_map = {}
        for dob, group in grouped:
            unique_gh_ids = sorted(group['GH_ID'].unique())
            base_pid = group['Base PID'].iloc[0]
            if len(unique_gh_ids) > 1:
                for i, gh_id in enumerate(unique_gh_ids):
                    suffix = chr(ord('A') + i)
                    pid_map[gh_id] = f"{base_pid}{suffix}"
            else:
                pid_map[unique_gh_ids[0]] = base_pid
        print("...done.")

        print("Step 4: Creating and saving master patient list...")
        master_df = df_combined[['GH_ID', 'Patient DOB']].drop_duplicates(subset=['GH_ID'])
        master_df['Effective Patient ID'] = master_df['GH_ID'].map(pid_map)
        master_df = master_df[['Effective Patient ID', 'GH_ID', 'Patient DOB']]
        master_df.to_csv(master_output, index=False)
        print(f"-> Master patient list saved to {master_output}")

        print("Step 5: Updating and saving cohort files...")
        
        print("  - Processing ILC cohort...")
        df_ilc['Effective Patient ID'] = df_ilc['GH_ID'].map(pid_map)
        cols_ilc = ['Effective Patient ID'] + [col for col in df_ilc if col != 'Effective Patient ID']
        df_ilc_final = df_ilc[cols_ilc]
        df_ilc_final.to_csv(ilc_output, index=False)
        print(f"  -> Updated ILC cohort saved to {ilc_output}")

        print("  - Processing IDC cohort...")
        df_idc['Effective Patient ID'] = df_idc['GH_ID'].map(pid_map)
        cols_idc = ['Effective Patient ID'] + [col for col in df_idc if col != 'Effective Patient ID']
        df_idc_final = df_idc[cols_idc]
        df_idc_final.to_csv(idc_output, index=False)
        print(f"  -> Updated IDC cohort saved to {idc_output}")

        print("\nScript finished successfully!")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    ILC_FILE = 'ilc_cohort.csv'
    IDC_FILE = 'idc_cohort.csv'
    MASTER_OUTPUT = 'patient_master_list.csv'
    ILC_OUTPUT = 'ilc_cohort_with_ids.csv'
    IDC_OUTPUT = 'idc_cohort_with_ids.csv'
    
    generate_definitive_patient_ids(ILC_FILE, IDC_FILE, MASTER_OUTPUT, ILC_OUTPUT, IDC_OUTPUT) 