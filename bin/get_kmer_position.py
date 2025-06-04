#!/usr/bin/env python
import argparse
import re
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_input_file(input_file):
    """Parse the input file and return the chromosome, position and kmer list."""
    with open(input_file, 'r') as f:
        line = f.readline().strip()
        chrom_pos, kmers = line.split('\t')
        chrom, pos = chrom_pos.split(',')
        pos = int(pos)
        kmer_list = kmers.split('; ')
    return chrom, pos, kmer_list

def find_kmer_positions(fasta_sequences, kmer):
    """Find all positions of a kmer in the sequences."""
    kmer_positions = []
    for chrom, sequence in fasta_sequences.items():
        for m in re.finditer(f'(?={kmer})', sequence):
            kmer_positions.append((chrom, m.start()))
    return kmer_positions

def process_kmer(kmer, chrom, pos, fasta_sequences):
    """Process a single kmer and return the BED lines."""
    kmer_positions = find_kmer_positions(fasta_sequences, kmer)
    kmer_length = len(kmer)
    bed_lines = []
    for kmer_chrom, kmer_pos in kmer_positions:
        if kmer_chrom == chrom and pos >= kmer_pos and pos <= kmer_pos + kmer_length:
            original_position_line = f'{kmer_chrom}\t{kmer_pos}\t{kmer_pos + kmer_length}\t'
            for other_chrom, other_pos in kmer_positions:
                if other_chrom != chrom or pos < other_pos or pos > other_pos + kmer_length:
                    bed_lines.append(f'{original_position_line}{other_chrom}\t{other_pos}\t{other_pos + kmer_length}')
    return bed_lines

def write_bed_files(chrom, pos, kmer_list, fasta_sequences, output_prefix):
    """Write the BED files for each kmer."""
    with ThreadPoolExecutor() as executor:
        futures = []
        for kmer_count, kmer in enumerate(kmer_list, start=1):
            future = executor.submit(process_kmer, kmer, chrom, pos, fasta_sequences)
            futures.append((future, kmer_count))
        
        for future, kmer_count in futures:
            bed_lines = future.result()
            output_file = f"{output_prefix}_{kmer_count}.bed"
            with open(output_file, 'w') as out_f:
                out_f.write('\n'.join(bed_lines) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Identify kmer positions in a FASTA file and output in BED format")
    parser.add_argument("input_file", help="Input file with positions and kmer sequences")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("output_prefix", help="Output file prefix for BED files")
    args = parser.parse_args()

    chrom, pos, kmer_list = parse_input_file(args.input_file)
    fasta_sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.fasta, "fasta")}
    write_bed_files(chrom, pos, kmer_list, fasta_sequences, args.output_prefix)

    print(f"BED files have been written with prefix '{args.output_prefix}'")

if __name__ == "__main__":
    main()

