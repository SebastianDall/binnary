import pytest
from src import data_processing

class MockArgs:
    def __init__(self):
        self.motifs_scored = "data/motifs-scored.tsv"
        self.bin_motifs = "data/bin-motifs.tsv"
        self.contig_bins = "data/bins.tsv"
        self.assembly_stats = "data/assembly_info.txt"
        self.assembly_file = "data/assembly_file.fasta"
        self.mean_methylation_cutoff = 0.5
        self.n_motif_cutoff = 6

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
    
    
@pytest.fixture(scope="session")
def motifs_scored_in_bins_and_bin_motifs(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]
    assembly_stats = loaded_data["assembly_stats"]
    
    # Alter bin_motifs to include motif_mod and mean
    bin_motifs["motif_mod"] = bin_motifs["motif"] + "_" + bin_motifs["mod_type"]
    # Calculate n_motifs and mean methylation
    bin_motifs["n_motifs"] = bin_motifs["n_mod_bin"] + bin_motifs["n_nomod_bin"]
    bin_motifs["mean"] = bin_motifs["n_mod_bin"] / bin_motifs["n_motifs"]
    
    
    
    # Step 2: create motifs_scored_in_bins
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored, bin_motifs, contig_bins, assembly_stats, 
    )
    
    return {
        "motifs_scored_in_bins": motifs_scored_in_bins,
        "bin_motifs": bin_motifs
    }