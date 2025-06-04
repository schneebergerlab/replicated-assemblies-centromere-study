import argparse
from Bio import SeqIO
# import matplotlib.pyplot as plt  # <-- Plotting turned off

def read_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        identifier = record.id
        sequence = str(record.seq)
        sequences[identifier] = sequence
    return sequences

def generate_word_hash(sequence, word_length):
    word_hash = {}
    for i in range(len(sequence) - word_length + 1):
        word = sequence[i:i+word_length]
        if word in word_hash:
            word_hash[word].append(i)
        else:
            word_hash[word] = [i]
    return word_hash

def count_mismatches(word1, word2):
    return sum(c1 != c2 for c1, c2 in zip(word1, word2))

def find_matching_words(seq1, seq2, word_length, max_mismatches):
    matches = []
    word_hash_seq1 = generate_word_hash(seq1, word_length)
    for i in range(len(seq2) - word_length + 1):
        word = seq2[i:i+word_length]
        if word in word_hash_seq1:
            for match_pos in word_hash_seq1[word]:
                mismatches = count_mismatches(word, seq1[match_pos:match_pos+word_length])
                if mismatches <= max_mismatches:
                    matches.append((match_pos, i, word, '+', mismatches))
    return matches

def find_complement_matches(seq1, seq2, word_length, max_mismatches):
    complement_matches = []
    word_hash_seq1 = generate_word_hash(seq1, word_length)
    for i in range(len(seq2) - word_length + 1):
        word = seq2[i:i+word_length]
        complement_word = reverse_complement(word)
        if complement_word in word_hash_seq1:
            for match_pos in word_hash_seq1[complement_word]:
                mismatches = count_mismatches(complement_word, seq1[match_pos:match_pos+word_length])
                if mismatches <= max_mismatches:
                    complement_matches.append((match_pos, i, complement_word, '-', mismatches))
    return complement_matches

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = ''.join(complement_dict.get(base, base) for base in reversed(sequence))
    return reverse_complement_seq

def write_matches_to_file(matches, complement_matches, output_file):
    with open(output_file, 'w') as file:
        file.write("Pos_Seq1\tPos_Seq2\tMatched_Word\tStrand\tMismatch_num\n")
        for match in matches:
            file.write(f"{match[0]}\t{match[1]}\t{match[2]}\t{match[3]}\t{match[4]}\n")
        for match in complement_matches:
            file.write(f"{match[0]}\t{match[1]}\t{match[2]}\t{match[3]}\t{match[4]}\n")

# def plot_dot_plot(matches, complement_matches, word_length, labels, sequences1, sequences2, output_file):
#     fig, ax = plt.subplots()
#     len_seq1 = len(list(sequences1.values())[0])
#     len_seq2 = len(list(sequences2.values())[0])
#     for match in matches:
#         x1, y1 = match[0] + 1, match[1] + 1
#         x2, y2 = x1 + word_length - 1, y1 + word_length - 1
#         ax.plot([x1, x2], [y1, y2], color='blue')
#     for match in complement_matches:
#         x1, y2 = match[0] + 1, match[1] + 1
#         x2, y1 = x1 + word_length - 1, y2 + word_length - 1
#         ax.plot([x1, x2], [y1, y2], color='red')
#     plt.xlim(0, len_seq1 + 1)
#     plt.ylim(0, len_seq2 + 1)
#     plt.xticks(range(1, len_seq1 + 1))
#     plt.yticks(range(1, len_seq2 + 1))
#     plt.title(f"Dot Plot (Word Length = {word_length})")
#     plt.xlabel(labels[0])
#     plt.ylabel(labels[1])
#     pdf_output_file = output_file.rsplit('.', 1)[0] + '.pdf'
#     plt.savefig(pdf_output_file)
#     plt.close()

def main():
    parser = argparse.ArgumentParser(description="Word-based alignment script")
    parser.add_argument("-l", "--word_length", type=int, default=3, help="Length of matching word (default: 3)")
    parser.add_argument("-i", "--max_mismatches", type=int, default=0, help="Maximum allowed mismatches in matching word (default: 0)")
    parser.add_argument("fasta1", type=str, help="First FASTA file")
    parser.add_argument("fasta2", type=str, help="Second FASTA file")
    parser.add_argument("-o", "--output_file", type=str, default="output.txt", help="Output file name (default: output.txt)")
    args = parser.parse_args()

    sequences1 = read_fasta(args.fasta1)
    sequences2 = read_fasta(args.fasta2)

    if not sequences1 or not sequences2:
        print("Error: At least one of the input FASTA files is empty.")
        return

    labels1 = list(sequences1.keys())
    labels2 = list(sequences2.keys())

    matches = find_matching_words(list(sequences1.values())[0], list(sequences2.values())[0], args.word_length, args.max_mismatches)
    complement_matches = find_complement_matches(list(sequences1.values())[0], list(sequences2.values())[0], args.word_length, args.max_mismatches)

    if matches or complement_matches:
        write_matches_to_file(matches, complement_matches, args.output_file)
        # plot_dot_plot(matches, complement_matches, args.word_length, [labels1[0], labels2[0]], sequences1, sequences2, args.output_file)
    else:
        print("No matches found.")

if __name__ == "__main__":
    main()
