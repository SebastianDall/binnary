import random

random.seed(1234)  # Set the random seed for reproducibility

def generate_random_sequence(length):
    """Generates a random DNA sequence of specified length."""
    return ''.join(random.choice('ACGT') for _ in range(length))

def write_fasta_file(num_contigs, file_name, seq_length=1000, line_width=60):
    """Writes a FASTA file with a specified number of contigs."""
    with open(file_name, 'w') as f:
        for i in range(1, num_contigs + 1):
            # Write the header
            f.write(f">contig_{i}\n")
            
            # Generate and write the sequence, wrapped to line_width
            sequence = generate_random_sequence(seq_length)
            for j in range(0, len(sequence), line_width):
                f.write(sequence[j:j+line_width] + "\n")

# Parameters
num_contigs = 15  # Number of contigs to generate
fasta_file = 'data/assembly_file.fasta'  # Output file name
contig_length = 500  # Length of each contig

# Generate the FASTA file
write_fasta_file(num_contigs, fasta_file, seq_length=contig_length)

print(f"Generated {num_contigs} contigs in '{fasta_file}'")
