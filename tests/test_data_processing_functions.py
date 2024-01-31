import pytest
from src import data_processing



def test_feature_with_loaded_data(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]
    assembly_stats = loaded_data["assembly_stats"]
    assembly_file = loaded_data["assembly_file"]
    
    
    
    # Now you can use the data in your test assertions
    assert motifs_scored is not None
    
    # bin motifs
    assert bin_motifs is not None
    assert set(bin_motifs["bin"].unique()) == {"b1", "b2", "b3"}
    
    # contig bins
    assert contig_bins is not None
    assert set(contig_bins["bin"].unique()) == {"b1", "b2", "b3", "b4"}
    
    contig_set = {f"contig_{i}" for i in range(1, 14)} | {"contig_16"}
    assert set(contig_bins["contig"].unique()) == contig_set

    
    # assembly stats
    assert assembly_stats is not None
    contig_set = {f"contig_{i}" for i in range(1, 17)}
    assert set(assembly_stats["contig"].unique()) == contig_set
    
    # assembly file
    assert assembly_file is not None
    assert len(assembly_file) == 16