
def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in "GCgc")
    total_length = len(sequence)
    return (gc_count / total_length) * 100 if total_length > 0 else 0

def parse_fasta(fasta_file):
    """Simple FASTA parser to extract sequences."""
    sequences = []
    current_sequence = ""
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
            else:
                current_sequence += line.strip()
        if current_sequence:
            sequences.append(current_sequence)
    return sequences

def calculate_assembly_stats(fasta_file):
    contig_lengths = []
    total_gc_count = 0
    total_bases = 0

    # Parse the FASTA file manually
    sequences = parse_fasta(fasta_file)
    for seq in sequences:
        length = len(seq)
        contig_lengths.append(length)
        total_gc_count += sum(1 for base in seq if base in "GCgc")
        total_bases += length

    # Sort contig lengths in descending order for N50 calculation
    contig_lengths.sort(reverse=True)

    # Calculate N50
    cumulative_length = 0
    n50 = 0
    for length in contig_lengths:
        cumulative_length += length
        if cumulative_length >= total_bases / 2:
            n50 = length
            break

    # Calculate GC content
    gc_content = (total_gc_count / total_bases) * 100 if total_bases > 0 else 0

    # Results
    stats = {
        "Total assembly size": total_bases,
        "Number of contigs": len(contig_lengths),
        "N50": n50,
        "GC content (%)": round(gc_content, 2)
    }

    return stats