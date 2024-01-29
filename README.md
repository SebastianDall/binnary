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
python3 binnary.py --motifs_scored data/motifs-scored.tsv \
    --bin_motifs data/bin-motifs.tsv \
    --contig_bins data/bins.tsv \
    --assembly_stats data/assembly_info.txt \
    --assembly_file data/assembly_file.fasta \
    --out output.tsv
```

This should produce the following output:

```bash
matching_bins	contig	current_bin
b1	contig_10	b1
b1	contig_1	b1
b1	contig_2	b1
b1	contig_3	b1
b2	contig_4	b2
b2	contig_5	b2
b3	contig_7	b3
b3	contig_8	b3
b3	contig_9	b3
b4	contig_16	b4
b1:0.05074253606808246|b3:0.05144902395625766	contig_11	b1
b1|b3	contig_15	unbinned
```
