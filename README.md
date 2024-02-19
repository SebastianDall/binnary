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

```bash
./binnary.py include_contigs --motifs_scored data/motifs-scored.tsv \
    --bin_motifs data/bin-motifs.tsv \
    --contig_bins data/bins.tsv \
    --assembly_stats data/assembly_info.txt \
    --out output.tsv
```


## Test 2
    

```bash
./binnary.py detect_contamination \
    --motifs_scored data/PaPr00000216MP_nm_0.1.18/nanomotif/motifs-scored.tsv \
    --bin_motifs data/PaPr00000216MP_nm_0.1.18/nanomotif/bin-motifs.tsv \
    --contig_bins data/PaPr00000216MP_nm_0.1.18/bins.tsv \
    --assembly_stats data/PaPr00000216MP_nm_0.1.18/flye/assembly_info.txt \
    --n_motif_contig_cutoff 10 \
    --out output/PaPr00000216MP_nm_0.1.18
```
--mean_methylation_cutoff 0.45 \
```bash
./binnary.py include_contigs \
    --motifs_scored data/PaPr00000216MP_nm_0.1.18/nanomotif/motifs-scored.tsv \
    --bin_motifs data/PaPr00000216MP_nm_0.1.18/nanomotif/bin-motifs.tsv \
    --contig_bins data/PaPr00000216MP_nm_0.1.18/bins.tsv \
    --assembly_stats data/PaPr00000216MP_nm_0.1.18/flye/assembly_info.txt \
    --run_detect_contamination \
    --out output/PaPr00000216MP_nm_0.1.18
```
