from .data_processing import prepare_motifs_scored_in_bins



def perform_analysis(motifs_scored, bin_motifs, contig_bins, assembly_stats, assembly_file, args):
    """
    Performs the core analysis of the tool. This function is called from the main entry point.
    """
    motifs_scored_in_bins = prepare_motifs_scored_in_bins(motifs_scored, bin_motifs, contig_bins, assembly_stats)
    
    return motifs_scored_in_bins