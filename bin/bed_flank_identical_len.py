import sys
from Bio import SeqIO
import argparse

def read_positions_file(positions_file):
    """
    Read positions file which contains 6 columns: chrom1, pos1s, pos1e, chrom2, pos2s, pos2e.
    Returns a list of tuples.
    """
    positions = []
    with open(positions_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) != 6:
                continue
            chrom1 = parts[0]
            pos1s = int(parts[1])
            pos1e = int(parts[2])
            chrom2 = parts[3]
            pos2s = int(parts[4])
            pos2e = int(parts[5])
            positions.append((chrom1, pos1s, pos1e, chrom2, pos2s, pos2e))
    return positions

def read_fasta_file(fasta_file):
    """
    Read FASTA file and return a dictionary where keys are chromosome names and values are sequences.
    """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def compare_sequences(seq1, seq2):
    """
    Compare two sequences and return the number of matching base pairs.
    """
    min_len = min(len(seq1), len(seq2))
    matches = 0
    for i in range(min_len):
        if seq1[i] == seq2[i]:
            matches += 1
        else:
            break
    return matches

def find_matching_sequences(positions, sequences, bed_output):
    """
    For each pair of positions, find the number of matching base pairs upstream and downstream.
    """
    results = []
    for chrom1, pos1s, pos1e, chrom2, pos2s, pos2e in positions:
        if chrom1 in sequences and chrom2 in sequences:
            seq1 = sequences[chrom1]
            seq2 = sequences[chrom2]

            # Compare upstream sequences
            upstream1 = seq1[:pos1s - 1][::-1]  # Reverse sequence for upstream comparison
            upstream2 = seq2[:pos2s - 1][::-1]
            upstream_matches = compare_sequences(upstream1, upstream2)

            # Compare downstream sequences
            downstream1 = seq1[pos1e:]
            downstream2 = seq2[pos2e:]
            downstream_matches = compare_sequences(downstream1, downstream2)

            if bed_output:
                new_pos1s = pos1s - upstream_matches
                new_pos1e = pos1e + downstream_matches
                new_pos2s = pos2s - upstream_matches
                new_pos2e = pos2e + downstream_matches
                results.append((chrom1, new_pos1s, new_pos1e, chrom2, new_pos2s, new_pos2e))
            else:
                results.append((chrom1, pos1s, pos1e, chrom2, pos2s, pos2e, upstream_matches, downstream_matches))
        else:
            print(f"Chromosome {chrom1} or {chrom2} not found in FASTA file.")
    return results

def write_results(results, output_file, bed_output):
    """
    Write results to an output file.
    """
    with open(output_file, 'w') as f:
        for result in results:
            if bed_output:
                f.write("\t".join(map(str, result)) + "\n")
            else:
                f.write("\t".join(map(str, result)) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare sequences around positions in a FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Input FASTA file")
    parser.add_argument("positions_file", type=str, help="Positions file")
    parser.add_argument("output_file", type=str, help="Output file")
    parser.add_argument("--bed", action="store_true", help="Output positions in BED format")

    args = parser.parse_args()

    positions_file = args.positions_file
    fasta_file = args.fasta_file
    output_file = args.output_file
    bed_output = args.bed

    positions = read_positions_file(positions_file)
    sequences = read_fasta_file(fasta_file)
    results = find_matching_sequences(positions, sequences, bed_output)
    write_results(results, output_file, bed_output)

    #print(f"Results saved successfully to {output_file}")

