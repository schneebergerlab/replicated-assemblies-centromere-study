import random
import numpy as np
import sys
import re
import pandas as pd
from Bio import pairwise2

def read_sequence(file_name):
    with open(file_name, "r") as file:
        sequence = file.read().strip()
    return sequence

def read_pos_file(file_path):
    """Read and parse a BED file, returning a DataFrame containing start positions with integer type."""
    return pd.read_csv(file_path, sep="\t", header=None, usecols=[1], names=["start"], dtype={"start": int})

def adjust_pos_coordinates(pos1, pos2):
    """Adjust pos1 coordinates based on pos2 instructions (INS/DEL)."""
    adjusted_pos1_temp1 = pos1

    for _, ins_del, p1, p2 in pos2:
        adjusted_pos1_temp2 = []
        for pos in adjusted_pos1_temp1:
            if ins_del == "DEL":
                if pos < p1:
                    adjusted_pos1_temp2.append(pos)  # no change to the positions smaller than p1
                elif p1 <= pos < p2:
                    continue  # delete positions between p1 and p2
                else:  # p2 <= pos
                    adjusted_pos1_temp2.append(pos - (p2 - p1))
            elif ins_del == "INS":
                if pos < p1:
                    adjusted_pos1_temp2.append(pos)  # no change to the positions smaller than p1
                elif p1 <= pos < p2:
                    adjusted_pos1_temp2.append(pos)  # insert/duplicate the positions between p1 and p2
                    adjusted_pos1_temp2.append(pos + (p2 - p1))
                else:  # p2 <= pos
                    adjusted_pos1_temp2.append(pos + (p2 - p1))
        adjusted_pos1_temp1 = adjusted_pos1_temp2
    return adjusted_pos1_temp1

def get_sequence(sequence, start, end):
    """Retrieve the sequence from the given sequence string, starting from 'start' position (0-based inclusive) to 'end' position (0-based exclusive)."""
    return "".join(sequence[start:end])  # Ensure the result is a string

def pairwise_alignment(seq1, seq2):
    """Perform pairwise alignment between two sequences using global alignment with affine gap penalties."""
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -10, -1)
    return alignments[0]

def find_pairwise_points(align, pos):
    pos_pairwise = -1  # Default to -1 if not found
    align_seq1, align_seq2 = align[0], align[1]
    align_pos1, align_pos2 = 0, 0

    for i in range(len(align_seq1)):
        if align_seq1[i] != '-':
            align_pos1 += 1
        if align_seq2[i] != '-':
            align_pos2 += 1
        if align_pos1 == pos:
            pos_pairwise = align_pos2
            break

    return pos_pairwise

def introduce_mutations(sequence, generation, num_generations, unit_data):
    mutated_sequence = list(sequence)
    mutation_records = []
    adjusted_pos = unit_data["start"].tolist()  # Corrected to use column name "start"

    for _ in range(num_generations):
        generation += 1

        actual_snps = 0
        num_snps = np.random.poisson(0.1)
        while actual_snps < num_snps:
            idx = random.randint(0, len(mutated_sequence) - 1)

            bases = ['A', 'T', 'C', 'G']
            bases.remove(mutated_sequence[idx])
            mutated_base = random.choice(bases)
            actual_snps += 1
            mutation_records.append((generation, "SNP", idx + 1, mutated_sequence[idx], mutated_base, "-"))
            mutated_sequence[idx] = mutated_base

        actual_indels = 0
        num_indels = np.random.poisson(0.5)
        while actual_indels < num_indels:
            idx = random.randint(0, len(mutated_sequence) - 1)
            indel_type = random.choice(["INS", "DEL"])

            copy_num = np.random.poisson(7.6)
            idx_unit_start = unit_data[unit_data["start"] <= idx].tail(1)["start"].values[0]
            idx_unit_end = unit_data[unit_data["start"] > idx].head(1)["start"].values[0]

            if pd.isna(idx_unit_start) or pd.isna(idx_unit_end):
                continue
            idx_unit_seq = get_sequence(mutated_sequence, idx_unit_start, idx_unit_end)

            if unit_data[unit_data["start"] > idx].head(copy_num).empty:
                continue
            idx_pairwise_unit_start = unit_data[unit_data["start"] > idx].head(copy_num).tail(1)["start"].values[0]

            if unit_data[unit_data["start"] > idx].head(copy_num + 1).empty:
                continue
            idx_pairwise_unit_end = unit_data[unit_data["start"] > idx].head(copy_num + 1).tail(1)["start"].values[0]

            if pd.isna(idx_pairwise_unit_start) or pd.isna(idx_pairwise_unit_end):
                continue
            idx_pairwise_unit_seq = get_sequence(mutated_sequence, idx_pairwise_unit_start, idx_pairwise_unit_end)

            if len(idx_unit_seq) == 0 or len(idx_pairwise_unit_seq) == 0:
                continue

            try:
                align = pairwise_alignment(idx_unit_seq, idx_pairwise_unit_seq)
            except IndexError:
                continue

            idx_pairwise = find_pairwise_points(align, idx - idx_unit_start)

            if idx_pairwise == -1:  # Skip if pairwise point was not found
                continue

            idx_pairwise_abs = int(idx_pairwise_unit_start) + idx_pairwise

            if indel_type == "INS":
                insertion_sequence = "".join(mutated_sequence[idx:idx_pairwise_abs])
                mutation_records.append((generation, indel_type, idx, mutated_sequence[idx - 1], "".join(mutated_sequence[idx - 1:idx_pairwise_abs]), copy_num))
                mutated_sequence[idx:idx] = list(insertion_sequence)
                actual_indels += 1
            elif indel_type == "DEL":
                mutation_records.append((generation, indel_type, idx, "".join(mutated_sequence[idx - 1:idx_pairwise_abs]), "".join(mutated_sequence[idx - 1]), copy_num))
                del mutated_sequence[idx:idx_pairwise_abs]
                actual_indels += 1

            # Update adjusted_pos
            indel_records = [(generation, indel_type, idx, idx_pairwise_abs)]
            adjusted_pos = sorted(adjust_pos_coordinates(adjusted_pos, indel_records))
            unit_data = pd.DataFrame(adjusted_pos, columns=["start"])

        actual_conversions = 0
        num_conversions = np.random.poisson(1)
        while actual_conversions < num_conversions:
            start = random.randint(0, len(mutated_sequence) - 1)
            conversion_size = np.random.poisson(20)
            end = start + conversion_size

            start_unit_start = unit_data[unit_data["start"] <= start].tail(1)["start"].values[0]
            if unit_data[unit_data["start"] > start].empty:
                continue
            start_unit_end = unit_data[unit_data["start"] > start].head(1)["start"].values[0]
            end_unit_start = unit_data[unit_data["start"] <= end].tail(1)["start"].values[0]
            if unit_data[unit_data["start"] > end].empty:
                continue
            end_unit_end = unit_data[unit_data["start"] > end].head(1)["start"].values[0]

            if pd.isna(start_unit_start) or pd.isna(start_unit_end) or pd.isna(end_unit_start) or pd.isna(end_unit_end):
                continue
            start_unit_seq = get_sequence(mutated_sequence, start_unit_start, start_unit_end)
            end_unit_seq = get_sequence(mutated_sequence, end_unit_start, end_unit_end)

            if unit_data[unit_data["start"] > start].head(1).empty:
                continue
            start_pairwise_unit_start = unit_data[unit_data["start"] > start].head(1).tail(1)["start"].values[0]

            if unit_data[unit_data["start"] > start].head(1 + 1).empty:
                continue
            start_pairwise_unit_end = unit_data[unit_data["start"] > start].head(1 + 1).tail(1)["start"].values[0]

            if unit_data[unit_data["start"] > end].head(1).empty:
                continue
            end_pairwise_unit_start = unit_data[unit_data["start"] > end].head(1).tail(1)["start"].values[0]

            if unit_data[unit_data["start"] > end].head(1 + 1).empty:
                continue
            end_pairwise_unit_end = unit_data[unit_data["start"] > end].head(1 + 1).tail(1)["start"].values[0]

            if pd.isna(start_pairwise_unit_start) or pd.isna(start_pairwise_unit_end) or pd.isna(end_pairwise_unit_start) or pd.isna(end_pairwise_unit_end):
                continue
            start_pairwise_unit_seq = get_sequence(mutated_sequence, start_pairwise_unit_start, start_pairwise_unit_end)
            end_pairwise_unit_seq = get_sequence(mutated_sequence, end_pairwise_unit_start, end_pairwise_unit_end)

            if len(start_unit_seq) == 0 or len(end_unit_seq) == 0 or len(start_pairwise_unit_seq) == 0 or len(end_pairwise_unit_seq) == 0:
                continue

            try:
                align1 = pairwise_alignment(start_unit_seq, start_pairwise_unit_seq)
                align2 = pairwise_alignment(end_unit_seq, end_pairwise_unit_seq)
            except IndexError:
                continue

            start_pairwise = find_pairwise_points(align1, start - start_unit_start)
            end_pairwise = find_pairwise_points(align2, end - end_unit_start)

            if start_pairwise == -1 or end_pairwise == -1:  # Skip if pairwise point was not found
                continue

            start_pairwise_abs = int(start_pairwise_unit_start) + start_pairwise
            end_pairwise_abs = int(end_pairwise_unit_start) + end_pairwise

            donor_sequence = "".join(mutated_sequence[start:end])
            receipt_sequence = "".join(mutated_sequence[start_pairwise_abs:end_pairwise_abs])
            if donor_sequence == receipt_sequence:
                conversion_out = "Identical"
            elif len(donor_sequence) == len(receipt_sequence):
                conversion_out = "SNP"
            else:
                conversion_out = "INDEL"
            mutation_records.append((generation, "Conversion", start_pairwise_abs, receipt_sequence, donor_sequence, conversion_out))

            actual_conversions += 1

            # Update adjusted_pos
            conversion_indel_records = []
            if conversion_out == "INDEL":
                if start_unit_start == end_unit_start and start_pairwise_unit_start == end_pairwise_unit_start:
                    if len(donor_sequence) > len(receipt_sequence): #INS
                        conversion_indel_records = [(generation, "INS", start_pairwise_abs, start_pairwise_abs + len(donor_sequence) - len(receipt_sequence))]
                    else: #DEL
                        conversion_indel_records = [(generation, "DEL", end_pairwise_abs + len(receipt_sequence) - len(donor_sequence), end_pairwise_abs)]
                else:
                    if (start_unit_end - start) > (start_pairwise_unit_end - start_pairwise_abs): #INS
                        conversion_indel_records = [(generation, "INS", start_pairwise_abs, start_pairwise_abs + (start_unit_end - start) - (start_pairwise_unit_end - start_pairwise_abs))]
                    elif (start_unit_end - start) < (start_pairwise_unit_end - start_pairwise_abs): #DEL
                        conversion_indel_records = [(generation, "DEL", start_pairwise_abs + (start_unit_end - start) - (start_pairwise_unit_end - start_pairwise_abs), start_pairwise_abs)]
                    if (end - end_unit_start) > (end_pairwise_abs - end_pairwise_unit_start): #INS
                        conversion_indel_records = [(generation, "INS", end_pairwise_abs, end_pairwise_abs + (end - end_unit_start) - (end_pairwise_abs - end_pairwise_unit_start))]
                    elif (end - end_unit_start) < (end_pairwise_abs - end_pairwise_unit_start): #DEL
                        conversion_indel_records = [(generation, "DEL", end_pairwise_abs + (end - end_unit_start) - (end_pairwise_abs - end_pairwise_unit_start), end_pairwise_abs)]

                adjusted_pos = sorted(adjust_pos_coordinates(adjusted_pos, conversion_indel_records))
                unit_data = pd.DataFrame(adjusted_pos, columns=["start"])

    return "".join(mutated_sequence), mutation_records, adjusted_pos

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python simulate.py input.fa num_generations path_to_unit_data_file output.fa mutation_record.txt adjusted_pos_output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    match = re.match(r'1\.fasta/(\d+)?generation\.out\.fa', input_file)
    generation = int(match.group(1)) if match else 0
    num_generations = int(sys.argv[2])
    unit_data_file = sys.argv[3]
    output_file = sys.argv[4]
    record_output = sys.argv[5]
    adjusted_pos_output = sys.argv[6]

    original_sequence = read_sequence(input_file)

    # Read the unit_data file
    unit_data = read_pos_file(unit_data_file)

    mutated_sequence, mutation_records, adjusted_pos = introduce_mutations(original_sequence, generation, num_generations, unit_data)

    with open(output_file, "w") as seq_f:
        seq_f.write(mutated_sequence + '\n')

    with open(record_output, "w") as record_f:
        for record in mutation_records:
            generation, type, idx, ref, mut, copy_num = record
            record_f.write(f"{generation}, {type}, {idx}, {ref}, {mut}, {copy_num}\n")

    with open(adjusted_pos_output, "w") as adjusted_f:
        for pos in sorted(adjusted_pos): # Ensure posiitons are sorted before writing to the file
            adjusted_f.write(f"centro_{generation}gen\t{pos}\n")

    print("Mutated sequence written to", output_file)
    print("Mutation records written to", record_output)
    print("Adjusted positions written to", adjusted_pos_output)

