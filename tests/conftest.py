import pytest
from src import data_processing

class MockArgs:
    def __init__(self):
        self.motifs_scored = "data/motifs-scored.tsv"
        self.bin_motifs = "data/bin-motifs.tsv"
        self.contig_bins = "data/bins.tsv"
        self.assembly_stats = "data/assembly_info.txt"
        self.assembly_file = "data/assembly_file.fasta"

@pytest.fixture(scope="session")
def loaded_data():
    # Instantiate MockArgs
    mock_args = MockArgs()

    # Load the data using the mock_args object
    data = data_processing.load_data(mock_args)

    # Unpack the data tuple to individual variables if needed
    motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file = data

    # Return the data as a dictionary or as individual variables, depending on your preference
    return {
        "motifs_scored": motifs_scored,
        "bin_motifs": bin_motifs,
        "contig_bins": contig_bins,
        "assembly_stats": assembly_stats,
        "assembly_file": assembly_file
    }