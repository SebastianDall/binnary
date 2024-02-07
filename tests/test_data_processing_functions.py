import pytest
from src import data_processing
from .conftest import MockArgs
import pandas as pd

def test_feature_with_loaded_data(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]
    assembly_stats = loaded_data["assembly_stats"]
    # assembly_file = loaded_data["assembly_file"]

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
    # assert assembly_file is not None
    # assert len(assembly_file) == 16



def test_prepare_bin_motif_binary(loaded_data):
    """
    GIVEN loaded_data
    WHEN prepare_bin_motifs_binary is called
    THEN assert that the output contains only the expected columns
    """
    args = MockArgs()
    
    bin_motif_binary = data_processing.prepare_bin_motifs_binary(loaded_data["bin_motifs"], args)
    
    assert bin_motif_binary is not None
    assert bin_motif_binary.columns.tolist() == ["bin", "motif_mod", "mean_methylation", "methylation_binary"]
    assert bin_motif_binary[(bin_motif_binary["bin"] == "b3") & (bin_motif_binary["motif_mod"] == "m6_a")]["methylation_binary"].values[0] == 1
    
    # Assert that there are 4 motifs in bin 1
    assert bin_motif_binary[bin_motif_binary["bin"] == "b3"].shape[0] == 2
    


def test_motifs_scored_in_bins(loaded_data):
    # Access the loaded data directly if returned as a dictionary
    motifs_scored = loaded_data["motifs_scored"]
    bin_motifs = loaded_data["bin_motifs"]
    contig_bins = loaded_data["contig_bins"]
    assembly_stats = loaded_data["assembly_stats"]

    # Step 1 create bin_motif_binary
    args = MockArgs()
    bin_motif_binary = data_processing.prepare_bin_motifs_binary(bin_motifs, args)

    # Step 2: create motifs_scored_in_bins
    motifs_scored_in_bins = data_processing.prepare_motifs_scored_in_bins(
        motifs_scored,
        bin_motif_binary,
        contig_bins,
        assembly_stats,
    )

    # motifs_scored_in_bins
    assert sorted(bin_motif_binary["motif_mod"].unique().tolist()) == sorted(
        motifs_scored_in_bins["motif_mod"].unique().tolist()
    )
    assert motifs_scored_in_bins.shape[1] == 13  # number of columns
    # Assert contig 1 belongs to bin 1
    assert (
        motifs_scored_in_bins.loc[
            motifs_scored_in_bins["contig"] == "contig_1", "bin"
        ].unique()
        == "b1"
    )


def test_contig_motif_binary_function(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data and motifs_scored_in_bins_and_bin_motifs
    WHEN construct_contig_motif_binary is called
    THEN assert that the output contains only the expected columns, and that there are 2 unique motifs and 16 unique contigs
    """
    bin_motif_binary = motifs_scored_in_bins_and_bin_motifs["bin_motif_binary"]
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    contig_motif_binary = data_processing.construct_contig_motif_binary(motifs_scored_in_bins, bin_motif_binary["motif_mod"].unique(), args.mean_methylation_cutoff-0.15, args.n_motif_cutoff)
    
    assert contig_motif_binary is not None
    assert contig_motif_binary.columns.tolist() == ["bin_compare", "motif_mod", "methylation_binary_compare"]
    
    # Assert there are two unique motifs in contig_motif_binary
    assert contig_motif_binary["motif_mod"].nunique() == 2
    
    # Assert that there are 2 bins in contig_motif_binary
    assert contig_motif_binary["bin_compare"].nunique() == 16


def test_bin_motifs_from_motifs_scored_in_bins(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    GIVEN loaded_data and motifs_scored_in_bins_and_bin_motifs
    WHEN construct_bin_motifs_from_motifs_scored_in_bins is called
    THEN assert that the output contains only the expected columns
    """
    bin_motif_binary = motifs_scored_in_bins_and_bin_motifs["bin_motif_binary"]
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    bin_motifs_from_motifs_scored_in_bins = data_processing.construct_bin_motifs_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        bin_motif_binary["motif_mod"].unique(),
        args
    )
    
    print(bin_motifs_from_motifs_scored_in_bins.columns.tolist())
    
    assert bin_motifs_from_motifs_scored_in_bins is not None
    assert bin_motifs_from_motifs_scored_in_bins.columns.tolist() == ['bin', 'motif_mod', 'n_mod', 'n_nomod', 'mean_methylation', 'methylation_binary']
    
    
    
    
def test_compare_methylation_pattern(loaded_data, motifs_scored_in_bins_and_bin_motifs):
    """
    
    """
    bin_motif_binary = motifs_scored_in_bins_and_bin_motifs["bin_motif_binary"]
    motifs_scored_in_bins = motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"]
    
    args = MockArgs()
    
    motifs_of_interest = bin_motif_binary["motif_mod"].unique()
    
    contig_motif_binary = data_processing.construct_contig_motif_binary(
        motifs_scored_in_bins, 
        motifs_of_interest,
        args.mean_methylation_cutoff-0.15, 
        args.n_motif_cutoff
    )
    
    bin_motifs_from_motifs_scored_in_bins = data_processing.construct_bin_motifs_from_motifs_scored_in_bins(
        motifs_scored_in_bins,
        motifs_of_interest,
        args
    )
    
     # Merge bin_motifs_from_motifs_scored_in_bins and contig_motif_binary    
    motif_binary_compare = pd.merge(
        bin_motifs_from_motifs_scored_in_bins, contig_motif_binary, on="motif_mod"
    )
    
    # Define the corresponding choices for each condition
    choices = [
        0,  # bin motif is methylated, contig motif is methylated
        1,  # bin motif is methylated, contig motif is not methylated
        1,  # bin motif is not methylated, contig motif is methylated
        0,  # bin motif is not methylated, contig motif is not methylated
        0,  # bin motif is methylated, contig motif is not observed
        0,  # bin motif is not methylated, contig motif is not observed
    ]
    
    contig_bin_comparison_score = data_processing.compare_methylation_pattern(motif_binary_compare, choices)
    
    print("\n")
    print(contig_bin_comparison_score)
    
    assert contig_bin_comparison_score is not None
    assert set(contig_bin_comparison_score["bin"].unique()) == {"b2", "b3"}
    assert len(contig_bin_comparison_score["contig"].unique()) == 16
    assert contig_bin_comparison_score[(contig_bin_comparison_score["bin"] == "b2") & (contig_bin_comparison_score["contig"] == "contig_14")]["binary_methylation_missmatch_score"].values[0] == 0
    
