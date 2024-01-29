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
```bash
python3 binnary.py --motifs_scored data/motifs-scored.tsv \
    --bin_motifs data/bin-motifs.tsv \
    --contig_bins data/bins.tsv \
    --assembly_stats data/assembly_info.txt \
    --assembly_file data/assembly_file.fasta \
    --out output.tsv
```
