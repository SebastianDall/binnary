import os
import pytest
import pandas as pd
from src.data_processing import generate_output

def test_generate_output_new_directory(tmp_path):
    """
    Test generate_output when the output directory does not exist.
    It should create the directory and the output file.
    """
    # Create a DataFrame to write to a file
    df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
    
    # Define the output path using a directory that does not exist
    output_path = tmp_path / "new_directory" / "output.tsv"
    
    # Call the function to test
    generate_output(df, str(output_path))
    
    # Check if the new directory was created
    assert output_path.parent.is_dir(), "Output directory was not created."
    
    # Check if the file was created and contains the expected content
    assert output_path.is_file(), "Output file was not created."
    # Read the content to verify
    result_df = pd.read_csv(output_path, sep="\t")
    pd.testing.assert_frame_equal(result_df, df, check_dtype=False)

def test_generate_output_existing_directory(tmp_path):
    """
    Test generate_output when the output directory already exists.
    """
    # Create a DataFrame to write to a file
    df = pd.DataFrame({'col1': [5, 6], 'col2': [7, 8]})
    
    # Use an existing directory
    output_path = tmp_path / "output.tsv"
    
    # Call the function to test
    generate_output(df, str(output_path))
    
    # Check if the file was created and contains the expected content
    assert output_path.is_file(), "Output file was not created in the existing directory."
    # Read the content to verify
    result_df = pd.read_csv(output_path, sep="\t")
    pd.testing.assert_frame_equal(result_df, df, check_dtype=False)

