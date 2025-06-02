from Bio import SeqIO
import Levenshtein
import sys

def read_external_sequence(external_fasta_file):
    """Read the first sequence from the external FASTA file."""
    with open(external_fasta_file, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            return str(record.seq)
    raise ValueError("No sequences found in the external FASTA file.")

def calculate_edit_distance(fasta_file, external_sequence, output_file):
    distances = []
    with open(fasta_file, "r") as fasta, open(output_file, "w") as output:
        output.write("Sequence_Name\tEdit_Distance\n")
        for record in SeqIO.parse(fasta, "fasta"):
            seq_name = record.id
            seq = str(record.seq)
            distance = Levenshtein.distance(seq, external_sequence)
            distances.append(distance)
            output.write(f"{seq_name}\t{distance}\n")
    
    # Calculate the average edit distance
    average_distance = sum(distances) / len(distances) if distances else 0
    
    # Print the average distance to the console
    print(f"{fasta_file}\t{average_distance:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <fasta_file> <external_fasta_file> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    external_fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    external_sequence = read_external_sequence(external_fasta_file)
    calculate_edit_distance(fasta_file, external_sequence, output_file)

