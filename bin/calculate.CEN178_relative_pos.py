import sys
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq  # Used for handling reverse complement sequences

def read_sequence(file_name):
    """Read sequence file and store each chromosome's sequence in a dictionary."""
    sequences = {}
    with open(file_name, "r") as file:
        current_chromosome = None
        sequence_lines = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_chromosome:
                    sequences[current_chromosome] = "".join(sequence_lines)
                current_chromosome = line[1:]  # Use the header line as chromosome name
                sequence_lines = []
            else:
                sequence_lines.append(line)
        if current_chromosome:  # Save the last chromosome sequence
            sequences[current_chromosome] = "".join(sequence_lines)
    return sequences

def read_pos_file(file_path):
    """Read and parse BED file, returning a DataFrame with chromosome and start positions."""
    return pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"], dtype={"chr": str, "start": int, "end": int})

def read_input_positions(file_path):
    """Read input file containing chromosome and position information."""
    return pd.read_csv(file_path, sep="\t", header=None, names=["chr", "position"], dtype={"chr": str, "position": int})

def get_sequence(sequences, chromosome, start, end):
    """Retrieve a substring from the sequence of the specified chromosome starting at 'start' (inclusive) and ending at 'end' (exclusive)."""
    if chromosome not in sequences:
        print(f"Chromosome {chromosome} not found in the sequence data.")
        sys.exit(1)
    return sequences[chromosome][start:end]

def pairwise_alignment(seq1, seq2):
    """Perform global alignment with affine gap penalties between two sequences."""
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -10, -1)
    return alignments[0]

def find_pairwise_points(align, pos, reverse_complement=False):
    """Find the corresponding position in CEN178 based on the input position in the aligned sequence."""
    pos_pairwise = -1  # Default to -1, indicating not found
    align_seq1, align_seq2 = align[0], align[1]
    align_pos1, align_pos2 = 0, 0

    if reverse_complement:
        pos = len(align_seq1.replace("-", "")) - pos + 1  # Adjust position for reverse complement sequences

    for i in range(len(align_seq1)):
        if align_seq1[i] != '-':
            align_pos1 += 1
        if align_seq2[i] != '-':
            align_pos2 += 1
        if align_pos1 == pos:
            pos_pairwise = align_pos2
            break

    return pos_pairwise

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py input.fa input_positions.txt path_to_unit_data_file output_file.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    input_positions_file = sys.argv[2]
    unit_data_file = sys.argv[3]
    output_file = sys.argv[4]

    # Fixed CEN178 sequence
    cen178_sequence = "AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG"

    # Read input sequence, unit_data file, and input positions file
    original_sequence = read_sequence(input_file)
    unit_data = read_pos_file(unit_data_file)
    input_positions = read_input_positions(input_positions_file)

    # Open the output file for writing
    with open(output_file, "w") as output:
        output.write("chr\tposition\tCEN178_position\tstrand\n")  # Write header line

        # Loop over each input position
        for _, row in input_positions.iterrows():
            chromosome = row["chr"]
            position = row["position"]

            # Filter unit_data to get only rows with the same chromosome
            unit_data_chr = unit_data[unit_data["chr"] == chromosome]

            if unit_data_chr.empty:
                print(f"No data found for chromosome {chromosome} in the unit data file.")
                continue

            # Get the unit coordinates that contain the input position
            if unit_data_chr[unit_data_chr["start"] > position].empty:
                print(f"No unit found that contains position {position} on chromosome {chromosome}")
                continue
            
            idx_unit_start = unit_data_chr[unit_data_chr["start"] < position].tail(1)["start"].values[0]
            idx_unit_end = unit_data_chr[unit_data_chr["end"] >= position].head(1)["end"].values[0]

            # Extract unit sequence
            unit_sequence = get_sequence(original_sequence, chromosome, idx_unit_start, idx_unit_end)

            # Check if the unit sequence is reverse complement to CEN178
            unit_seq_obj = Seq(unit_sequence)
            reverse_complement_unit_sequence = str(unit_seq_obj.reverse_complement())

            if pairwise2.align.globalms(unit_sequence, cen178_sequence, 2, -1, -10, -1)[0][2] >= pairwise2.align.globalms(reverse_complement_unit_sequence, cen178_sequence, 2, -1, -10, -1)[0][2]:
                # If direct alignment score is higher
                align = pairwise_alignment(unit_sequence, cen178_sequence)
                reverse_complement = False
            else:
                # If reverse complement alignment score is higher
                align = pairwise_alignment(reverse_complement_unit_sequence, cen178_sequence)
                reverse_complement = True

            # Calculate the corresponding position in CEN178
            position_in_cen178 = find_pairwise_points(align, position - idx_unit_start, reverse_complement)

            if position_in_cen178 == -1:
                print(f"Position {position} on chromosome {chromosome} could not be mapped to CEN178")
            else:
                strand = "-" if reverse_complement else "+"
                output.write(f"{chromosome}\t{position}\t{position_in_cen178}\t{strand}\n")

