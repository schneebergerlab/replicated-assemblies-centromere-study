from Bio import SeqIO
import sys
from collections import Counter

def read_external_sequence(external_fasta_file):
    """Read the first sequence from the external FASTA file."""
    with open(external_fasta_file, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            return str(record.seq)
    raise ValueError("No sequences found in the external FASTA file.")

def analyze_sequences(fasta_file, external_sequence):
    seq_counter = Counter()
    total_count = 0

    with open(fasta_file, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            seq = str(record.seq)
            seq_counter[seq] += 1
            total_count += 1
    
    match_count = seq_counter[external_sequence]
    match_ratio = match_count / total_count if total_count > 0 else 0

    del seq_counter[external_sequence]
    most_common_seq, most_common_count = seq_counter.most_common(1)[0] if seq_counter else ("", 0)
    most_common_ratio = most_common_count / total_count if total_count > 0 else 0

    return total_count, match_count, match_ratio, most_common_seq, most_common_count, most_common_ratio

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <fasta1_file> <fasta2_file>")
        print("Output format: ")
        print("fasta1_file\ttotal_count\tmatch_count\tmatch_ratio\tmost_common_seq\tmost_common_count\tmost_common_ratio")
        sys.exit(1)

    fasta1_file = sys.argv[1]
    fasta2_file = sys.argv[2]

    external_sequence = read_external_sequence(fasta2_file)
    total_count, match_count, match_ratio, most_common_seq, most_common_count, most_common_ratio = analyze_sequences(fasta1_file, external_sequence)

    # Output results to one line
    print(f"{fasta1_file}\t{total_count}\t{match_count}\t{match_ratio:.4f}\t{most_common_seq}\t{most_common_count}\t{most_common_ratio:.4f}")

