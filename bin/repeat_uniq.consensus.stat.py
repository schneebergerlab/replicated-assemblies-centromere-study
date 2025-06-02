import sys
from collections import Counter

def parse_fasta(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        seq = ''
        for line in file:
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                    seq = ''
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return sequences

def calculate_statistics(sequences):
    num_sequences = len(sequences)
    sequence_length = len(sequences[0])
    
    counts = []
    for i in range(sequence_length):
        column = [seq[i].upper() for seq in sequences]  # Convert to uppercase
        count = Counter(column)
        counts.append(count)
    
    return counts, num_sequences

def generate_output(counts, num_sequences, output_file):
    with open(output_file, 'w') as file:
        header = ["Position", "A_count", "T_count", "G_count", "C_count", "-_count", 
                  "Total", "A_freq", "T_freq", "G_freq", "C_freq", "-_freq", "Consensus"]
        file.write('\t'.join(header) + '\n')
        
        for i, count in enumerate(counts):
            position = i + 1
            A_count = count.get('A', 0)
            T_count = count.get('T', 0)
            G_count = count.get('G', 0)
            C_count = count.get('C', 0)
            gap_count = count.get('-', 0)
            
            total_count = A_count + T_count + G_count + C_count + gap_count
            
            A_freq = A_count / num_sequences
            T_freq = T_count / num_sequences
            G_freq = G_count / num_sequences
            C_freq = C_count / num_sequences
            gap_freq = gap_count / num_sequences
            
            freqs = {
                'A': A_freq,
                'T': T_freq,
                'G': G_freq,
                'C': C_freq,
                '-': gap_freq
            }
            
            consensus = max(freqs, key=freqs.get)
            if freqs[consensus] < 0.5:
                consensus = '-'
            
            file.write(f"{position}\t{A_count}\t{T_count}\t{G_count}\t{C_count}\t{gap_count}\t{total_count}\t"
                       f"{A_freq:.2f}\t{T_freq:.2f}\t{G_freq:.2f}\t{C_freq:.2f}\t{gap_freq:.2f}\t"
                       f"{consensus}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_fasta output_txt")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_txt = sys.argv[2]

    sequences = parse_fasta(input_fasta)
    counts, num_sequences = calculate_statistics(sequences)
    generate_output(counts, num_sequences, output_txt)

