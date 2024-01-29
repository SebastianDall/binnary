import pandas as pd
import os

# List of Excel files to convert
excel_files = [
    'bin-motifs.xlsx',
    'bins.xlsx',
    'motifs-scored.xlsx',
    'test_data.xlsx'
]

# Directory where the Excel files are located
input_dir = 'data'  # Adjust this to the path where your Excel files are

# Directory where you want to save the TSV files
output_dir = 'data'  # Adjust this to where you want to save the TSV files

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate over each file and convert it
for file_name in excel_files:
    # Construct full file path
    input_file_path = os.path.join(input_dir, file_name)
    # Generate output file name by replacing '.xlsx' with '.tsv'
    output_file_name = file_name.replace('.xlsx', '.tsv')
    output_file_path = os.path.join(output_dir, output_file_name)

    # Read the Excel file
    df = pd.read_excel(input_file_path)

    # Write to a TSV file, using tab as the separator
    df.to_csv(output_file_path, sep='\t', index=False)

    print(f"Converted '{input_file_path}' to '{output_file_path}'")
