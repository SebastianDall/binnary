import pytest
from src import scoring as sc
from src import data_processing as dp
import polars as pl

from .conftest import MockArgs

def test_compare_methylation_pattern(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    
    """
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    motifs_scored_in_bins = motifs_scored_in_bins \
        .filter(~pl.col("bin_contig").str.contains("unbinned"))
    
    
    
    # Define the corresponding choices for each condition
    choices = [
        0,  # bin motif is methylated, contig motif is methylated
        1,  # bin motif is methylated, contig motif is not methylated
        1,  # bin motif is not methylated, contig motif is methylated
        0,  # bin motif is not methylated, contig motif is not methylated
        0,  # bin motif is methylated, contig motif is not observed
        0,  # bin motif is not methylated, contig motif is not observed
    ]
    
    motif_binary_compare = sc.compare_methylation_pattern_multiprocessed(
        motifs_scored_in_bins,
        args
    )
    
    
    contig_bin_comparison_score = sc.compare_methylation_pattern(motif_binary_compare, choices)
    
    contig_bin_comparison_score = contig_bin_comparison_score.to_pandas()
    
    assert contig_bin_comparison_score is not None
    assert set(contig_bin_comparison_score["bin"].unique()) == {'b1', 'b2', 'b3', 'b4'}
    assert len(contig_bin_comparison_score["contig"].unique()) == 14
    assert contig_bin_comparison_score[(contig_bin_comparison_score["bin"] == "b2") & (contig_bin_comparison_score["contig"] == "contig_4")]["binary_methylation_missmatch_score"].values[0] == 0
    
