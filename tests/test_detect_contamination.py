import pytest
from src import detect_contamination
from src import data_processing as dp
from .conftest import MockArgs

def test_detect_contamination(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data
    WHEN detect_contamination is called
    THEN assert that the output contains only contig 3 and the expected columns
    """
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    
    bin_motifs_from_motifs_scored_in_bins = dp.construct_bin_consensus_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        args
    )
    
    contaminated_contigs = detect_contamination.detect_contamination(
        motifs_scored_in_bins,
        bin_motifs_from_motifs_scored_in_bins,
        args
    )
    
    assert contaminated_contigs is not None
    assert contaminated_contigs["bin"].unique().tolist() == ["b3"]
    assert sorted(contaminated_contigs["contig"].unique().tolist()) == ["contig_12", "contig_13", "contig_6"]
    assert contaminated_contigs[contaminated_contigs["contig"] == "contig_6"]["binary_methylation_missmatch_score"].values[0] == 3.0
    assert contaminated_contigs[contaminated_contigs["contig"] == "contig_6"]["alternative_bin"].unique() == ["b2"]
    