import pytest
from src.cli_parser import get_parser


def test_argparser_detect_contamination_with_specific_data(mocker):
    # Mock sys.argv with the specific command line you provided
    mocker.patch(
        "sys.argv",
        [
            "./binnary.py",  # Assuming 'binnary.py' is a typo and should be 'binary.py' or similar
            "detect_contamination",
            "--motifs_scored",
            "data/motifs-scored.tsv",
            "--bin_motifs",
            "data/bin-motifs.tsv",
            "--contig_bins",
            "data/bins.tsv",
            "--assembly_stats",
            "data/assembly_info.txt",
            # "--assembly_file", "data/assembly_file.fasta",
            "--out",
            "output.tsv",
        ],
    )

    # Parse the arguments
    parser = get_parser()
    args = parser.parse_args()

    # Assertions to verify that all arguments are correctly parsed
    assert args.command == "detect_contamination"
    assert args.motifs_scored == "data/motifs-scored.tsv"
    assert args.bin_motifs == "data/bin-motifs.tsv"
    assert args.contig_bins == "data/bins.tsv"
    assert args.assembly_stats == "data/assembly_info.txt"
    # assert args.assembly_file == "data/assembly_file.fasta"
    assert args.out == "output.tsv"
