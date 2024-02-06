# binnairy

Binary methylation pattern binning

## Install

### Conda Environment

In the devcontainer, run the following commands to create a conda environment with the required dependencies:

```bash
micromamba create -n nanomotif -c conda-forge python=3.9
micromamba activate nanomotif
#python -m pip install nanomotif
```

## Test with mock data

Test data was made by hand and the mock assembly file with the script in `src/utilities/assembly_file_creator.py`.

To run a test, run the following command:

```bash
./binnary.py detect_contamination --motifs_scored data/motifs-scored.tsv \
    --bin_motifs data/bin-motifs.tsv \
    --contig_bins data/bins.tsv \
    --assembly_stats data/assembly_info.txt \
    --out output.tsv
```

## Test 2
    

```bash
./binnary.py detect_contamination \
    --motifs_scored data/PaPr00000216MP/nanomotif/motifs-scored.tsv \
    --bin_motifs data/PaPr00000216MP/nanomotif/bin-motifs.tsv \
    --contig_bins data/PaPr00000216MP/bins.tsv \
    --assembly_stats data/PaPr00000216MP/flye/assembly_info.txt \
    --out binnary-PaPr00000216MP.tsv
```
