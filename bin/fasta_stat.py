import sys
from Bio import SeqIO

def fasta_statistics(fasta_file):
    lengths = []
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            lengths.append(len(record.seq))

    total_sequences = len(lengths)
    total_length = sum(lengths)
    average_length = total_length / total_sequences if total_sequences > 0 else 0

    #print(f"Filename\tTotal Length\tTotal Sequences\tAverage Length")
    print(f"{fasta_file}\t{total_length}\t{total_sequences}\t{average_length:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    fasta_statistics(fasta_file)

