from .data_processing import prepare_motifs_scored_in_bins



def perform_analysis(motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file, args):
    """
    Performs the core analysis of the tool. This function is called from the main entry point.
    """
    # Find relevant motifs in bins and contigs
    motifs_scored_in_bins = prepare_motifs_scored_in_bins(motifs_scored, bin_motifs, contig_bins, assembly_stats)
    
    # Calculate belonging score
    belonging_score = calculate_belonging_score(motifs_scored_in_bins, args)
    
    
    return belonging_score





def calculate_belonging_score(motifs_scored_in_bins, args):
    """
    Finds the closest match between contig and bin based on methylation pattern.
    """
    # Step 1: Group by bin and motif_mod and calculate mean methylation
    bin_motif_binary = motifs_scored_in_bins[motifs_scored_in_bins['bin'] != "unbinned"].groupby(['bin', 'motif_mod'])['mean'].mean().reset_index(name='mean_methylation')

    # Step 2: Convert mean methylation values to binary
    bin_motif_binary['methylation_binary'] = (bin_motif_binary['mean_methylation'] >= args.mean_methylation_cutoff).astype(int)

    # Step 4: (Optional) Merge the binary methylation back to the original DataFrame if needed
    # This step is optional and depends on whether you need these values associated with the original DataFrame rows
    # motifs_scored_in_bins_with_binary = motifs_scored_in_bins.merge(bin_motif_binary[['bin', 'motif_mod', 'methylation_binary']], on=['bin', 'motif_mod'], how='left')
    
    return bin_motif_binary


