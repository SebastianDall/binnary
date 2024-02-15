import pytest
from src import data_processing

class MockArgs:
    def __init__(self):
        self.motifs_scored = "data/test_data/motifs-scored.tsv"
        self.bin_motifs = "data/test_data/bin-motifs.tsv"
        self.contig_bins = "data/test_data/bins.tsv"
        self.assembly_stats = "data/test_data/assembly_info.txt"
        self.mean_methylation_cutoff = 0.25
        self.n_motif_bin_cutoff = 500
        self.n_motif_contig_cutoff = 10
        self.ambiguous_motif_percentage_cutoff = 0.40

@pytest.fixture(scope="session")
def loaded_data():
    # Instantiate MockArgs
    mock_args = MockArgs()

    # Load the data using the mock_args object
    data = data_processing.load_data(mock_args)

    # Unpack the data tuple to individual variables if needed
    motifs_scored, bin_motifs, contig_bins, assembly_stats = data #assembly_file

    # Return the data as a dictionary or as individual variables, depending on your preference
    return {
        "motifs_scored": motifs_scored,
        "bin_motifs": bin_motifs,
        "contig_bins": contig_bins,
        "assembly_stats": assembly_stats
    }
    
    
@pytest.fixture(scope="session")
def motifs_scored_in_bins_and_bin_motifs(loaded_data):
    args = MockArgs()
    
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]
    assembly_stats = loaded_data["assembly_stats"]
    
    bin_motif_binary = data_processing.prepare_bin_motifs_binary(bin_motifs, args)
    
    motifs_in_bins = bin_motif_binary["motif_mod"].unique()
    
    # Step 2: create motifs_scored_in_bins
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored, motifs_in_bins, contig_bins, assembly_stats, 
    )
    
    return {
        "motifs_scored_in_bins": motifs_scored_in_bins,
        "bin_motif_binary": bin_motif_binary
    }