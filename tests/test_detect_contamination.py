import pytest
from src import detect_contamination
from .conftest import MockArgs

def test_detect_contamination(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data
    WHEN detect_contamination is called
    THEN assert that the output contains only contig 3 and the expected columns
    """
    args = MockArgs()
    
    contaminated_contigs = detect_contamination.detect_contamination(
        motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"],
        motifs_scored_in_bins_and_bin_motifs["bin_motifs"],
        args
    )
    
    assert contaminated_contigs is not None
    assert contaminated_contigs["bin"].unique() == ["b3"]
    assert contaminated_contigs["contig"].unique() == ["contig_6"]
    assert contaminated_contigs["binary_methylation_missmatch_score"].unique() == [2.0]
    assert contaminated_contigs["alternative_bin"].unique() == ["b2"]
    