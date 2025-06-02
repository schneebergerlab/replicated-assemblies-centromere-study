import sys

def extract_consensus(input_file, output_file):
    consensus_sequence = []

    with open(input_file, 'r') as file:
        # Skip the header line
        next(file)
        
        for line in file:
            parts = line.strip().split('\t')
            consensus = parts[-1]
            
            if consensus != '-':
                consensus_sequence.append(consensus)

    # Join the consensus sequence into a single string
    consensus_sequence = ''.join(consensus_sequence)

    # Calculate the length of the consensus sequence
    consensus_length = len(consensus_sequence)

    # Write the consensus sequence and its length to the output file
    with open(output_file, 'w') as file:
        file.write(f'{consensus_sequence}\t{consensus_length}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_txt output_txt")
        sys.exit(1)

    input_txt = sys.argv[1]
    output_txt = sys.argv[2]

    extract_consensus(input_txt, output_txt)

