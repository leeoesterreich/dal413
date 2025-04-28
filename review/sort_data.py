import pandas as pd

# Read the CSV file
df = pd.read_csv('IO_summary.csv')

# Sort the dataframe by multiple columns in the specified order
df_sorted = df.sort_values(
    by=['Biomarker', 'Molecular subtype', 'Stage', 'IO Drug', 'Primary Endpoint', 'Stratification'],
    ascending=[True, True, True, True, True, True]
)

# Save the sorted data back to CSV
df_sorted.to_csv('IO_summary_sorted.csv', index=False)

print("Data has been sorted and saved to IO_summary_sorted.csv") 