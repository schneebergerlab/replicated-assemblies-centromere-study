import csv
import argparse

def read_bed_file(bed_file):
    bed_dict = {}
    with open(bed_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 3:
                continue  # Skip rows with less than 3 columns
            chr_bed, start, end = row[:3]  # Only take the first three columns
            if chr_bed not in bed_dict:
                bed_dict[chr_bed] = []
            bed_dict[chr_bed].append((int(start), int(end)))
    return bed_dict

def find_bed_range(chr_bed, pos, bed_dict):
    if chr_bed not in bed_dict:
        return None
    for start, end in bed_dict[chr_bed]:
        if start <= pos <= end:
            return start, end
    return None

def process_files(bed_file, txt_file, output_file):
    bed_dict = read_bed_file(bed_file)
    
    with open(txt_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = [row for row in reader if len(row) >= 9]  # Ensure rows have at least 9 columns

    # Sort rows by the absolute value of distance
    rows.sort(key=lambda x: abs(int(x[7])))

    for row in rows:
        chr1, start1, end1, chr2, start2, end2, size, distance, same_chr = row
        start1 = int(start1)
        start2 = int(start2)
        distance = int(distance)
        
        if same_chr == "YES":
            bed_range1 = find_bed_range(chr1, start1, bed_dict)
            bed_range2 = find_bed_range(chr2, start2, bed_dict)

            if bed_range1 and bed_range2:
                relative_pos1 = (start1 - bed_range1[0]) / (bed_range1[1] - bed_range1[0])
                relative_pos2 = (start2 - bed_range2[0]) / (bed_range2[1] - bed_range2[0])

                if abs(relative_pos1 - relative_pos2) <= 0.02:  # Difference less than or equal to 2%
                    with open(output_file, 'w') as out_f:
                        writer = csv.writer(out_f, delimiter='\t')
                        writer.writerow(row)
                    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process bed and txt files.")
    parser.add_argument("bed_file", help="Path to the input bed file.")
    parser.add_argument("txt_file", help="Path to the input txt file.")
    parser.add_argument("output_file", help="Path to the output file.")
    
    args = parser.parse_args()
    
    process_files(args.bed_file, args.txt_file, args.output_file)

