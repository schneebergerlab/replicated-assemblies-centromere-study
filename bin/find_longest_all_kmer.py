import sys
import re
import time
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed

def read_fasta_file(fasta_file):
    """Reads a FASTA file and returns a dictionary with sequence IDs as keys and sequences as values."""
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict

def read_positions_file(positions_file):
    """Reads a positions file and returns a list of (chromosome, position) tuples. Converts 1-based to 0-based positions."""
    positions = []
    with open(positions_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:  # Ensure there are at least two parts (chrom and pos)
                chrom = parts[0]
                pos = int(parts[1]) - 1  # Convert 1-based to 0-based
                positions.append((chrom, pos))
            else:
                print(f"Skipping invalid line: {line.strip()}")
    return positions

def find_kmer_positions(fasta_dict, chrom, kmer, exclude_pos):
    """Finds all positions of a given kmer in the sequence, excluding the specified position."""
    positions = []
    for ch, seq in fasta_dict.items():
        if ch == chrom:
            for match in re.finditer(f'(?=({re.escape(kmer)}))', seq):
                start = match.start(1)
                if start <= exclude_pos < start + len(kmer):
                    continue
                positions.append((ch, start))
    return positions

def find_all_valid_kmers(chrom, seq, pos, fasta_dict, max_k, min_k):
    """Finds all kmers from a given position that have at least one hit in non-original positions."""
    exclude_pos = pos
    valid_kmers = []
    for k in range(max_k, min_k - 1, -1):
        kmer_start = max(0, pos - k + 1)
        kmer_end = min(len(seq), pos + k)
        for i in range(kmer_start, kmer_end - k + 1):
            kmer = seq[i:i + k]
            kmer_positions = find_kmer_positions(fasta_dict, chrom, kmer, exclude_pos)
            if kmer_positions:
                valid_kmers.append(kmer)
    return valid_kmers

def process_position(chrom, pos, fasta_dict, max_k, min_k):
    """Processes a single position to find all valid kmers with at least one non-original position hit."""
    if chrom in fasta_dict:
        seq = fasta_dict[chrom]
        valid_kmers = find_all_valid_kmers(chrom, seq, pos, fasta_dict, max_k, min_k)
        if valid_kmers:
            result = "; ".join(valid_kmers)
        else:
            result = "No_valid_kmers_found"
        return f"{chrom},{pos + 1}\t{result}"  # Convert back to 1-based for output
    else:
        return f"{chrom},{pos + 1}\tChromosome not found in FASTA file"  # Convert back to 1-based for output

def process_files(fasta_file, positions_file, output_file, max_k, min_k):
    """Main function to process files and output results."""
    start_time = time.time()  # Record start time

    fasta_dict = read_fasta_file(fasta_file)
    positions = read_positions_file(positions_file)

    results = []
    with ThreadPoolExecutor(max_workers=20) as executor:  # Adjust max_workers as needed
        futures = []
        for chrom, pos in positions:
            futures.append(executor.submit(process_position, chrom, pos, fasta_dict, max_k, min_k))
        for future in as_completed(futures):
            results.append(future.result())

    with open(output_file, 'w') as out_file:
        for result in results:
            out_file.write(result + '\n')

    end_time = time.time()  # Record end time
    print(f"Execution finished in {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <fasta_file> <positions_file> <output_file> <max_k> <min_k>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    positions_file = sys.argv[2]
    output_file = sys.argv[3]
    max_k = int(sys.argv[4])
    min_k = int(sys.argv[5])

    process_files(fasta_file, positions_file, output_file, max_k, min_k)

