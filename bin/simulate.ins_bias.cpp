#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <random>
#include <algorithm>
#include <cmath>
#include <regex>
#include <map>

using namespace std;

// Random number generator
random_device rd;
mt19937 gen(rd());

// Poisson distribution helper
int poisson_sample(double lambda) {
    poisson_distribution<int> dist(lambda);
    return dist(gen);
}

// Geometric distribution helper
int geometric_sample(double p) {
    geometric_distribution<int> dist(p);
    return dist(gen) + 1; // +1 because geometric_distribution starts from 0
}

// Uniform random integer
int random_int(int min, int max) {
    uniform_int_distribution<int> dist(min, max);
    return dist(gen);
}

// Structure for mutation records
struct MutationRecord {
    int generation;
    string type;
    int idx;
    string ref;
    string mut;
    string copy_num;
};

// Structure for indel records
struct IndelRecord {
    int generation;
    string type;
    int p1;
    int p2;
};

// Read sequence from file
string read_sequence(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    
    string sequence;
    string line;
    while (getline(file, line)) {
        sequence += line;
    }
    
    // Remove whitespace
    sequence.erase(remove_if(sequence.begin(), sequence.end(), ::isspace), sequence.end());
    
    file.close();
    return sequence;
}

// Read position file (BED format)
vector<int> read_pos_file(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    
    vector<int> positions;
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string col1;
        int start;
        
        if (iss >> col1 >> start) {
            positions.push_back(start);
        }
    }
    
    file.close();
    return positions;
}

// Adjust positions based on INS/DEL operations
vector<int> adjust_pos_coordinates(const vector<int>& pos1, const vector<IndelRecord>& pos2) {
    vector<int> adjusted = pos1;
    
    for (const auto& record : pos2) {
        vector<int> temp;
        
        for (int pos : adjusted) {
            if (record.type == "DEL") {
                if (pos < record.p1) {
                    temp.push_back(pos);
                } else if (pos >= record.p2) {
                    temp.push_back(pos - (record.p2 - record.p1));
                }
                // Skip positions between p1 and p2
            } else if (record.type == "INS") {
                if (pos < record.p1) {
                    temp.push_back(pos);
                } else if (pos >= record.p1 && pos < record.p2) {
                    temp.push_back(pos);
                    temp.push_back(pos + (record.p2 - record.p1));
                } else {
                    temp.push_back(pos + (record.p2 - record.p1));
                }
            }
        }
        
        adjusted = temp;
    }
    
    return adjusted;
}

// Get subsequence
string get_sequence(const string& sequence, int start, int end) {
    if (start < 0 || end > sequence.length() || start >= end) {
        return "";
    }
    return sequence.substr(start, end - start);
}

// Simple pairwise alignment (NW algorithm)
pair<string, string> pairwise_alignment(const string& seq1, const string& seq2) {
    int match = 2;
    int mismatch = -1;
    int gap_open = -10;
    int gap_extend = -1;
    
    int len1 = seq1.length();
    int len2 = seq2.length();
    
    // DP matrix
    vector<vector<int>> score(len1 + 1, vector<int>(len2 + 1, 0));
    
    // Initialize
    for (int i = 0; i <= len1; i++) {
        score[i][0] = gap_open + i * gap_extend;
    }
    for (int j = 0; j <= len2; j++) {
        score[0][j] = gap_open + j * gap_extend;
    }
    
    // Fill matrix
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int match_score = score[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch);
            int delete_score = score[i-1][j] + gap_extend;
            int insert_score = score[i][j-1] + gap_extend;
            
            score[i][j] = max({match_score, delete_score, insert_score});
        }
    }
    
    // Traceback
    string align1, align2;
    int i = len1, j = len2;
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && score[i][j] == score[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)) {
            align1 = seq1[i-1] + align1;
            align2 = seq2[j-1] + align2;
            i--;
            j--;
        } else if (i > 0 && score[i][j] == score[i-1][j] + gap_extend) {
            align1 = seq1[i-1] + align1;
            align2 = '-' + align2;
            i--;
        } else {
            align1 = '-' + align1;
            align2 = seq2[j-1] + align2;
            j--;
        }
    }
    
    return {align1, align2};
}

// Find pairwise points
int find_pairwise_points(const pair<string, string>& align, int pos) {
    const string& align_seq1 = align.first;
    const string& align_seq2 = align.second;
    
    int align_pos1 = -1;
    int align_pos2 = -1;
    
    for (size_t i = 0; i < align_seq1.length(); i++) {
        if (align_seq1[i] != '-') {
            align_pos1++;
        }
        if (align_seq2[i] != '-') {
            align_pos2++;
        }
        
        if (align_pos1 == pos) {
            if (align_seq2[i] == '-') {
                return -1;
            }
            return align_pos2;
        }
    }
    
    return -1;
}

// Find value in sorted vector (lower bound)
int find_lower_bound(const vector<int>& vec, int value) {
    auto it = lower_bound(vec.begin(), vec.end(), value);
    if (it != vec.end() && *it ==value) {
	    return *it;
    }

    if (it == vec.begin()) return -1;
    return *(--it);
}

// Find value in sorted vector (upper bound)
int find_upper_bound(const vector<int>& vec, int value) {
    auto it = upper_bound(vec.begin(), vec.end(), value);
    if (it == vec.end()) return -1;
    return *it;
}

// Get nth element from iterator
int get_nth_from_position(const vector<int>& vec, int start_value, int n) {
    auto it = upper_bound(vec.begin(), vec.end(), start_value);
    if (it == vec.end()) return -1;
    
    for (int i = 0; i < n && it != vec.end(); i++, it++);
    
    if (it == vec.end()) return -1;
    return *it;
}

// Main mutation function
tuple<string, vector<MutationRecord>, vector<int>> introduce_mutations(
    const string& sequence,
    int generation,
    int num_generations,
    vector<int> unit_data)
{
    string mutated_sequence = sequence;
    vector<MutationRecord> mutation_records;
    vector<int> adjusted_pos = unit_data;

    for (int gen = 0; gen < num_generations; gen++) {

        if (mutated_sequence.length() < 50000) {
            cerr << "Sequence too short (<50kb), stopping simulation at generation "
                 << generation << endl;
            break;
        }
	else if (mutated_sequence.length() > 500000000) {
		cerr << "Sequence too long (>100Mb), stopping simulation at generation "
			<< generation << endl;
            break;
        }

        generation++;

        // ========== SNP ==========
        int num_snps = poisson_sample(mutated_sequence.length() * 1.97e-8);
        int snp_count = 0, attempt = 0, max_attempt = num_snps * 10;

        while (snp_count < num_snps && attempt < max_attempt) {
            attempt++;

            int idx = random_int(0, mutated_sequence.length() - 1);

            vector<char> bases = {'A','T','C','G'};
            bases.erase(remove(bases.begin(), bases.end(),
                               mutated_sequence[idx]), bases.end());

            char mutated_base =
                bases[random_int(0, bases.size() - 1)];

            MutationRecord record;
            record.generation = generation;
            record.type = "SNP";
            record.idx = idx + 1;
            record.ref = string(1, mutated_sequence[idx]);
            record.mut = string(1, mutated_base);
            record.copy_num = "-";

            mutation_records.push_back(record);
            mutated_sequence[idx] = mutated_base;

            snp_count++;
        }

        // ========== INSERTIONS (tandem duplication) ==========
        int num_insertions = poisson_sample(mutated_sequence.length() * 1.15e-8);

	int ins_count = 0; attempt = 0; max_attempt = num_insertions * 50;

        while (ins_count < num_insertions && attempt < max_attempt) {

            attempt++;

            int idx = random_int(0, mutated_sequence.length() - 1);
            int copy_num = max(1, poisson_sample(10.51));

            // ----- current unit -----
            int idx_unit_start = find_lower_bound(adjusted_pos, idx);
            int idx_unit_end   = find_upper_bound(adjusted_pos, idx);

            if (idx_unit_start == -1 || idx_unit_end == -1) continue;

            // ----- donor unit -----
            int donor_unit_start =
                get_nth_from_position(adjusted_pos, idx, copy_num - 1);

            int donor_unit_end =
                get_nth_from_position(adjusted_pos, idx, copy_num);

            if (donor_unit_start == -1 || donor_unit_end == -1) continue;

            // ----- alignment -----
            int offset = idx - idx_unit_start;

            string current_unit_seq =
                get_sequence(mutated_sequence,
                             idx_unit_start,
                             idx_unit_end);

            string donor_unit_seq =
                get_sequence(mutated_sequence,
                             donor_unit_start,
                             donor_unit_end);

            if (current_unit_seq.empty() || donor_unit_seq.empty())
                continue;

            auto align =
                pairwise_alignment(current_unit_seq,
                                   donor_unit_seq);

            int mapped_offset =
                find_pairwise_points(align, offset);

            if (mapped_offset == -1)
                continue;

            int idx_pairwise_abs =
                donor_unit_start + mapped_offset;

            if (idx_pairwise_abs <= idx)
                continue;

	    size_t ins_size = idx_pairwise_abs - idx;

	    if (mutated_sequence.length() + ins_size > 500000000)
		    continue;

            // ----- duplication sequence -----
            string insertion_sequence =
                mutated_sequence.substr(idx,
                                        idx_pairwise_abs - idx);

            if (insertion_sequence.empty())
                continue;

            // ----- mutation record -----
            string ref_base, mut_sequence;

            if (idx == 0) {
                ref_base = "";
                mut_sequence = insertion_sequence;
            } else {
                ref_base = string(1, mutated_sequence[idx - 1]);
                mut_sequence =
                    mutated_sequence.substr(idx - 1,
                                            idx_pairwise_abs - (idx - 1));
            }

            MutationRecord record;
            record.generation = generation;
            record.type = "INS";
            record.idx = idx;
            record.ref = ref_base;
            record.mut = mut_sequence;
            record.copy_num = to_string(copy_num);
            mutation_records.push_back(record);

            // ----- apply insertion -----
            mutated_sequence.insert(idx, insertion_sequence);

            // ----- adjust coordinates -----
            vector<IndelRecord> indel_records =
                {{generation, "INS", idx, idx_pairwise_abs}};

            adjusted_pos =
                adjust_pos_coordinates(adjusted_pos, indel_records);

            sort(adjusted_pos.begin(), adjusted_pos.end());

            // ----- success -----
            ins_count++;
        }

        // ========== DELETIONS ==========
        int num_deletions = poisson_sample(mutated_sequence.length() * 6.93e-9);
        int del_count = 0; attempt = 0; max_attempt = num_deletions * 100;

        while (del_count < num_deletions && attempt < max_attempt) {
            attempt++;

            int idx = random_int(0, mutated_sequence.length() - 1);
            int copy_num = max(1, poisson_sample(5.19));

            int idx_unit_start = find_lower_bound(adjusted_pos, idx);
            int idx_unit_end   = find_upper_bound(adjusted_pos, idx);
            if (idx_unit_start == -1 || idx_unit_end == -1) continue;

            int target_unit_start =
                get_nth_from_position(adjusted_pos, idx, copy_num - 1);
            int target_unit_end =
                get_nth_from_position(adjusted_pos, idx, copy_num);
            if (target_unit_start == -1 || target_unit_end == -1) continue;

            int offset = idx - idx_unit_start;

            string current_unit_seq =
                get_sequence(mutated_sequence, idx_unit_start, idx_unit_end);
            string target_unit_seq =
                get_sequence(mutated_sequence, target_unit_start, target_unit_end);

            if (current_unit_seq.empty() || target_unit_seq.empty())
                continue;

            auto align =
                pairwise_alignment(current_unit_seq, target_unit_seq);

            int mapped_offset =
                find_pairwise_points(align, offset);
            if (mapped_offset == -1) continue;

            int idx_pairwise_abs =
                target_unit_start + mapped_offset;
            if (idx_pairwise_abs <= idx) continue;

	    size_t del_size = idx_pairwise_abs - idx;

	    if (mutated_sequence.length() - del_size < 1000000)
		    continue;

            string ref_sequence, alt_base;

            if (idx == 0) {
                ref_sequence =
                    mutated_sequence.substr(0, idx_pairwise_abs);
                alt_base = "";
            } else {
                ref_sequence =
                    mutated_sequence.substr(idx - 1,
                                            idx_pairwise_abs - (idx - 1));
                alt_base = string(1, mutated_sequence[idx - 1]);
            }

            MutationRecord record;
            record.generation = generation;
            record.type = "DEL";
            record.idx = idx;
            record.ref = ref_sequence;
            record.mut = alt_base;
            record.copy_num = to_string(copy_num);
            mutation_records.push_back(record);

            mutated_sequence.erase(idx, idx_pairwise_abs - idx);

            vector<IndelRecord> indel_records =
                {{generation, "DEL", idx, idx_pairwise_abs}};

            adjusted_pos =
                adjust_pos_coordinates(adjusted_pos, indel_records);

            sort(adjusted_pos.begin(), adjusted_pos.end());

            del_count++;
        }

        // ========== LARGE DELETIONS ==========
        //int num_large_deletions = poisson_sample(0.004);
	int num_large_deletions = 0;
        int ldel_count = 0; attempt = 0; max_attempt = num_large_deletions * 200;

        while (ldel_count < num_large_deletions && attempt < max_attempt) {
            attempt++;

            int idx = random_int(0, mutated_sequence.length() - 1);
            int copy_num = max(1, poisson_sample(264.0));

            int idx_unit_start = find_lower_bound(adjusted_pos, idx);
            int idx_unit_end   = find_upper_bound(adjusted_pos, idx);
            if (idx_unit_start == -1 || idx_unit_end == -1) continue;

            int target_unit_start =
                get_nth_from_position(adjusted_pos, idx, copy_num - 1);
            int target_unit_end =
                get_nth_from_position(adjusted_pos, idx, copy_num);
            if (target_unit_start == -1 || target_unit_end == -1) continue;

            int offset = idx - idx_unit_start;

            string current_unit_seq =
                get_sequence(mutated_sequence, idx_unit_start, idx_unit_end);
            string target_unit_seq =
                get_sequence(mutated_sequence, target_unit_start, target_unit_end);

            if (current_unit_seq.empty() || target_unit_seq.empty())
                continue;

            auto align =
                pairwise_alignment(current_unit_seq, target_unit_seq);

            int mapped_offset =
                find_pairwise_points(align, offset);
            if (mapped_offset == -1) continue;

            int idx_pairwise_abs =
                target_unit_start + mapped_offset;
            if (idx_pairwise_abs <= idx) continue;

	    size_t del_size = idx_pairwise_abs - idx;

	    if (mutated_sequence.length() - del_size < 1000000)
		    continue;
	    
	    string ref_sequence, alt_base;
	    
	    if (idx == 0) {
		    ref_sequence =
			    mutated_sequence.substr(0, idx_pairwise_abs);
		    alt_base = "";
	    } else {
		    ref_sequence =
			    mutated_sequence.substr(idx - 1,
					    idx_pairwise_abs - (idx - 1));
		    
		    alt_base = string(1, mutated_sequence[idx - 1]);
	    }

            MutationRecord record;
            record.generation = generation;
            record.type = "Large_DEL";
            record.idx = idx;
            record.ref = ref_sequence;
            record.mut = alt_base;
            record.copy_num = to_string(copy_num);
            mutation_records.push_back(record);

            mutated_sequence.erase(idx, idx_pairwise_abs - idx);

            vector<IndelRecord> indel_records =
                {{generation, "DEL", idx, idx_pairwise_abs}};

            adjusted_pos =
                adjust_pos_coordinates(adjusted_pos, indel_records);

            sort(adjusted_pos.begin(), adjusted_pos.end());

            ldel_count++;
        }

	// ========== CONVERSION ==========
	
	int num_conversions = poisson_sample(mutated_sequence.length() * 3.39e-8);

	int conv_count = 0; attempt = 0; max_attempt = num_conversions * 100;

	while (conv_count < num_conversions && attempt < max_attempt) {
		attempt++;

        	int start = random_int(0, mutated_sequence.length() - 1);
        	int conversion_size = geometric_sample(0.05);
        	int end = start + conversion_size;

       		if (end >= mutated_sequence.length()) continue;

        	// ----- donor start unit -----
        	int start_unit_start = find_lower_bound(adjusted_pos, start);
        	int start_unit_end   = find_upper_bound(adjusted_pos, start);

        	// ----- donor end unit -----
        	int end_unit_start = find_lower_bound(adjusted_pos, end);
        	int end_unit_end   = find_upper_bound(adjusted_pos, end);

        	if (start_unit_start == -1 || start_unit_end == -1 ||
				end_unit_start   == -1 || end_unit_end   == -1)
			continue;

        	// ----- receptor units (next unit) -----
        	int start_pairwise_unit_start =
			get_nth_from_position(adjusted_pos, start, 0);

        	int start_pairwise_unit_end =
			get_nth_from_position(adjusted_pos, start, 1);
		
		int end_pairwise_unit_start =
			get_nth_from_position(adjusted_pos, end, 0);

        	int end_pairwise_unit_end =
			get_nth_from_position(adjusted_pos, end, 1);

        	if (start_pairwise_unit_start == -1 || start_pairwise_unit_end == -1 ||
				end_pairwise_unit_start   == -1 || end_pairwise_unit_end   == -1)
			continue;

        	// ----- sequences -----
        	string start_unit_seq =
			get_sequence(mutated_sequence,
					start_unit_start,
					start_unit_end);

        	string end_unit_seq =
			get_sequence(mutated_sequence,
					end_unit_start,
					end_unit_end);

        string start_pairwise_unit_seq =
            get_sequence(mutated_sequence,
                         start_pairwise_unit_start,
                         start_pairwise_unit_end);

        string end_pairwise_unit_seq =
            get_sequence(mutated_sequence,
                         end_pairwise_unit_start,
                         end_pairwise_unit_end);

        if (start_unit_seq.empty() || end_unit_seq.empty() ||
            start_pairwise_unit_seq.empty() || end_pairwise_unit_seq.empty())
            continue;

        // ----- alignment -----
        auto align1 =
            pairwise_alignment(start_unit_seq,
                               start_pairwise_unit_seq);

        auto align2 =
            pairwise_alignment(end_unit_seq,
                               end_pairwise_unit_seq);

        int start_offset = start - start_unit_start;
        int end_offset   = end   - end_unit_start;

        int start_pairwise =
            find_pairwise_points(align1, start_offset);

        int end_pairwise =
            find_pairwise_points(align2, end_offset);

        if (start_pairwise == -1 || end_pairwise == -1)
            continue;

        int start_pairwise_abs =
            start_pairwise_unit_start + start_pairwise;

        int end_pairwise_abs =
            end_pairwise_unit_start + end_pairwise;

        if (end_pairwise_abs <= start_pairwise_abs)
            continue;

        // ----- donor / receptor sequences -----
        string donor_sequence =
            mutated_sequence.substr(start, end - start);

        string receptor_sequence =
            mutated_sequence.substr(start_pairwise_abs,
                                    end_pairwise_abs - start_pairwise_abs);

        if (donor_sequence.empty() || receptor_sequence.empty())
            continue;

        if (donor_sequence.length() != receptor_sequence.length())
            continue;

        // ----- record conversion type -----
        string conversion_out;

        if (donor_sequence == receptor_sequence)
            conversion_out = "Identical";
        else if (donor_sequence.length() == receptor_sequence.length())
            conversion_out = "SNP";
        else
            conversion_out = "INDEL";

        MutationRecord record;
        record.generation = generation;
        record.type = "Conversion";
        record.idx = start_pairwise_abs;
        record.ref = receptor_sequence;
        record.mut = donor_sequence;
        record.copy_num = conversion_out;
        mutation_records.push_back(record);

        // ===== APPLY CONVERSION =====
        int donor_len    = donor_sequence.length();
        int receptor_len = receptor_sequence.length();

        mutated_sequence.replace(
            start_pairwise_abs,
            receptor_len,
            donor_sequence);

        // ===== ADJUST COORDINATES IF LENGTH CHANGED =====
        if (donor_len != receptor_len) {

            vector<IndelRecord> indel_records;

            if (donor_len > receptor_len) {

                indel_records.push_back(
                    {generation, "INS",
                     start_pairwise_abs,
                     start_pairwise_abs + donor_len});

            } else {

                indel_records.push_back(
                    {generation, "DEL",
                     start_pairwise_abs,
                     start_pairwise_abs + receptor_len});
            }

            adjusted_pos =
                adjust_pos_coordinates(adjusted_pos, indel_records);

            sort(adjusted_pos.begin(), adjusted_pos.end());
        }

        conv_count++;
    }

    //cout<<"Gen "<<generation<<" CONV "<<conv_count<<"/"<<num_conversions<<" att="<<attempt<<endl;
    
    }

    return make_tuple(mutated_sequence, mutation_records, adjusted_pos);
}

int main(int argc, char* argv[]) {
    if (argc != 7) {
        cerr << "Usage: " << argv[0] << " input.fa num_generations path_to_unit_data_file output.fa mutation_record.txt adjusted_pos_output.txt" << endl;
        return 1;
    }
    
    string input_file = argv[1];
    int num_generations = stoi(argv[2]);
    string unit_data_file = argv[3];
    string output_file = argv[4];
    string record_output = argv[5];
    string adjusted_pos_output = argv[6];
    
    // Extract generation from filename
    int generation = 0;
    regex pattern(R"(1\.fasta/(\d+)?generation\.out\.fa)");
    smatch match;
    if (regex_match(input_file, match, pattern)) {
        if (match[1].matched) {
            generation = stoi(match[1]);
        }
    }
    
    // Read input
    string original_sequence = read_sequence(input_file);
    vector<int> unit_data = read_pos_file(unit_data_file);
    
    // Perform mutations
    auto result = introduce_mutations(original_sequence, generation, num_generations, unit_data);
    string mutated_sequence = get<0>(result);
    vector<MutationRecord> mutation_records = get<1>(result);
    vector<int> adjusted_pos = get<2>(result);
    
    // Write output sequence
    ofstream seq_f(output_file);
    if (!seq_f.is_open()) {
        cerr << "Error: Cannot write to " << output_file << endl;
        return 1;
    }
    seq_f << mutated_sequence << "\n";
    seq_f.close();
    
    // Write mutation records
    ofstream record_f(record_output);
    if (!record_f.is_open()) {
        cerr << "Error: Cannot write to " << record_output << endl;
        return 1;
    }
    for (const auto& record : mutation_records) {
	//if(record.type != "SNP") {
	//	cout << "Non-SNP detected this run: " 
	//	<< record.type 
	//	<< " at generation " 
	//	<< record.generation << endl;
	//}

        record_f << record.generation << ", " << record.type << ", " 
                 << record.idx << ", " << record.ref << ", " 
                 << record.mut << ", " << record.copy_num << "\n";
    }
    record_f.close();
    
    // Write adjusted positions
    ofstream adjusted_f(adjusted_pos_output);
    if (!adjusted_f.is_open()) {
        cerr << "Error: Cannot write to " << adjusted_pos_output << endl;
        return 1;
    }
    sort(adjusted_pos.begin(), adjusted_pos.end());
    for (int pos : adjusted_pos) {
        adjusted_f << "centro_" << generation << "gen\t" << pos << "\n";
    }
    adjusted_f.close();
    
    cout << "Mutated sequence written to " << output_file << endl;
    cout << "Mutation records written to " << record_output << endl;
    cout << "Adjusted positions written to " << adjusted_pos_output << endl;
    
    return 0;
}
