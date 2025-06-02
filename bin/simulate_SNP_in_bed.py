import random
import sys
from Bio import SeqIO
import pandas as pd

def read_bed_file(bed_file):
    """Read and parse the BED file, returning a dictionary of allowed mutation ranges for each sequence."""
    bed_data = pd.read_csv(bed_file, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"], dtype={"chr": str, "start": int, "end": int})
    mutation_ranges = {}
    for _, row in bed_data.iterrows():
        if row["chr"] not in mutation_ranges:
            mutation_ranges[row["chr"]] = []
        mutation_ranges[row["chr"]].append((row["start"], row["end"]))
    return mutation_ranges

def mutate_sequence(sequence, ranges, snp_count):
    """Generate SNP mutations for a given sequence within specific ranges."""
    nucleotides = ['A', 'T', 'C', 'G']
    mutations = []
    mutated_sequence = list(sequence)

    possible_positions = []
    for start, end in ranges:
        possible_positions.extend(range(start, end))

    if len(possible_positions) < snp_count:
        raise ValueError("Not enough positions available to introduce the requested number of SNPs.")

    for _ in range(snp_count):
        position = random.choice(possible_positions)
        original_allele = mutated_sequence[position]
        possible_alleles = [n for n in nucleotides if n != original_allele]
        new_allele = random.choice(possible_alleles)
        mutated_sequence[position] = new_allele
        mutations.append((position + 1, original_allele, new_allele))  # Adjusting position to be 1-based
        possible_positions.remove(position)  # Remove this position to prevent duplicate mutations

    return "".join(mutated_sequence), mutations

def write_fasta(sequences, output_fasta):
    """Write the sequences to a fasta file."""
    with open(output_fasta, 'w') as f:
        for record_id, sequence in sequences.items():
            f.write(f">{record_id}\n")
            f.write(sequence + "\n")

def write_mutations(mutations, output_mutations):
    """Write the mutations to a file."""
    with open(output_mutations, 'w') as f:
        f.write("SequenceID\tPosition\tOriginal\tMutated\n")
        for seq_id, pos, original, mutated in mutations:
            f.write(f"{seq_id}\t{pos}\t{original}\t{mutated}\n")

def main(input_fasta, bed_file, snp_count, output_fasta, output_mutations):
    # Parse the input fasta file
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(input_fasta, "fasta")}
    
    # Parse the BED file
    mutation_ranges = read_bed_file(bed_file)

    mutations = []
    mutated_sequences = sequences.copy()

    for _ in range(snp_count):
        seq_id = random.choice(list(mutated_sequences.keys()))

        if seq_id not in mutation_ranges:
            print(f"No mutation ranges found for sequence {seq_id}. Skipping.")
            continue

        sequence = mutated_sequences[seq_id]
        mutated_sequence, snp_mutations = mutate_sequence(sequence, mutation_ranges[seq_id], 1)

        for pos, original, mutated in snp_mutations:
            mutations.append((seq_id, pos, original, mutated))

        mutated_sequences[seq_id] = mutated_sequence

    # Write the output files
    write_fasta(mutated_sequences, output_fasta)
    write_mutations(mutations, output_mutations)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py input.fasta input.bed snp_count output.fasta output_mutations.txt")
    else:
        input_fasta = sys.argv[1]
        bed_file = sys.argv[2]
        snp_count = int(sys.argv[3])
        output_fasta = sys.argv[4]
        output_mutations = sys.argv[5]
        main(input_fasta, bed_file, snp_count, output_fasta, output_mutations)

