# import pytest
# from src import include_contigs
# from .conftest import MockArgs

# def test_include_contigs(loaded_data, motifs_scored_in_bins_and_bin_motifs):
#     """
#     GIVEN loaded_data
#     WHEN include_contigs is called
#     THEN assert that the output contains only contig 3 and the expected columns
#     """
#     args = MockArgs()
#     print("\n")
#     print(motifs_scored_in_bins_and_bin_motifs["bin_motif_binary"])
#     new_bin_contigs = include_contigs.include_contigs(
#         motifs_scored_in_bins_and_bin_motifs["motifs_scored_in_bins"],
#         motifs_scored_in_bins_and_bin_motifs["bin_motif_binary"]["motif_mod"].unique(),
#         args
#     )
#     print("\n")
#     print(new_bin_contigs)
#     assert new_bin_contigs is not None
#     assert new_bin_contigs["bin"].unique().tolist() == ["b2"]
#     assert new_bin_contigs["contig"].unique().tolist() == ["contig_14"]
    
    