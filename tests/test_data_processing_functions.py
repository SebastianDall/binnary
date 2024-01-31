import pytest
from src import data_processing


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


def test_motifs_scored_in_bins(loaded_data):
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
        motifs_scored,
        bin_motifs,
        contig_bins,
        assembly_stats,
    )

    # bin_motifs
    ## Assert there are only the m1, m2, m3 in motif_mod
    assert sorted(bin_motifs["motif_mod"].unique()) == sorted(
        ["m1_a", "m2_a", "m3_a", "m6_a"]
    )

    # motifs_scored_in_bins
    print(motifs_scored_in_bins)
    assert sorted(bin_motifs["motif_mod"].unique().tolist()) == sorted(
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
