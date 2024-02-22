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
    --out output.tsv
```

```bash
./binnary.py include_contigs --motifs_scored data/motifs-scored.tsv \
    --bin_motifs data/bin-motifs.tsv \
    --contig_bins data/bins.tsv \
    --out output.tsv
```

## Test 2

```bash
./binnary.py detect_contamination \
    --motifs_scored data/PaPr00000216MP_nm_0.1.18/nanomotif/motifs-scored.tsv \
    --bin_motifs data/PaPr00000216MP_nm_0.1.18/nanomotif/bin-motifs.tsv \
    --contig_bins data/PaPr00000216MP_nm_0.1.18/bins.tsv \
    --n_motif_contig_cutoff 10 \
    --threads 4 \
    --out output/PaPr00000216MP_nm_0.1.18_test
```

```bash
./binnary.py include_contigs \
    --threads 4 \
    --motifs_scored data/PaPr00000216MP_nm_0.1.18/nanomotif/motifs-scored.tsv \
    --bin_motifs data/PaPr00000216MP_nm_0.1.18/nanomotif/bin-motifs.tsv \
    --contig_bins data/PaPr00000216MP_nm_0.1.18/bins.tsv \
    --run_detect_contamination \
    --assembly_file data/PaPr00000216MP_nm_0.1.18/flye/assembly_file.fasta \
    --write_bins \
    --out output/PaPr00000216MP_nm_0.1.18_test
```
