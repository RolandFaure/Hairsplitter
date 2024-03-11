#include "separate_reads.h"
#include "cluster_graph.h"
#include "tools.h"

#include "input_output.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <chrono> //for time measurement
#include <iostream>
#include <fstream> //for reading files
#include <sstream> //for parsing strings
#include <string>
#include <omp.h> //for efficient parallelization
#include <map>
#include <set>
#include <cmath> //for sqrt and pow
#include <algorithm> //for sort and std::find
#include <ctime>

using namespace std::chrono;
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::atoi;
using std::make_pair;
using std::pair;
using std::set;
using robin_hood::unordered_map;
using std::max;
using std::min;

/**
 * @brief Parse a column file and store the SNPs for each contig
 * 
 * @param file File to parse
 * @param snps Associates each contig with its SNPs
 * @param index_of_names Associates each contig with its index in snps
 * @param name_of_contigs Associates index of contigs with their name
 * @param numberOfReads Associates each contig with the number of reads aligned on it
 */
void parse_column_file(
    std::string file, 
    std::vector<std::vector<Column>>  &snps, 
    std::unordered_map<string, int> &index_of_names,
    std::unordered_map<int, string> &name_of_contigs,
    std::vector<std::vector<string>> &names_of_reads,
    std::vector<long int> &length_of_contigs,
    std::vector<std::vector<std::pair<int,int>>> &readLimits,
    std::vector<int>& numberOfReads,
    int max_coverage){

    std::ifstream infile(file);
    std::string line;
    int nbreadshere = 0;
    names_of_reads = vector<vector<string>>(0);
    bool numbers = false; //are snps encoded as numbers or as chars ?
    bool firstsnpline = true;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        std::string line_type;
        iss >> line_type;
        if (line_type == "CONTIG"){
            index_of_names[line] = snps.size();
            name_of_contigs[snps.size()] = line;
            //parse the length of the contig, which is the third field of the line
            string name;
            string length;
            iss >> name >> length;
            length_of_contigs.push_back(std::atoi(length.c_str()));
            snps.push_back(std::vector<Column>());
            numberOfReads.push_back(0);
            names_of_reads.push_back(vector<string>(0));
            readLimits.push_back(vector<pair<int,int>>(0));
        }
        else if (line_type == "SNPS"){
            string pos;
            string ref_base_string;
            string second_frequent_base_string;
            char ref_base;
            char second_frequent_base;
            string content;
            iss >> pos >> ref_base_string >> second_frequent_base_string;

            if (firstsnpline && (!std::isalpha(ref_base_string[0]) && ref_base_string[0] != '-')){
                numbers = true;
            }
            if (numbers){
                ref_base = (char) std::stoi(ref_base_string);
                second_frequent_base = (char) std::stoi(second_frequent_base_string);
            }
            else{
                ref_base = ref_base_string[0];//the string should be of length 1
                second_frequent_base = second_frequent_base_string[0];
            }
            firstsnpline = false;
            //get the content and the indices
            string readsIdx;
            iss >> readsIdx >> content;
            //content is a comma-separated list of the bases at this position represented by ints or chars: remove the commas and convert to char
            string new_content = "";
            string integer = ""; //the variant can be either an integer encoding a char or directly a char
            for (auto c : content){
                if (c == ','){
                    if (integer == " "){
                        new_content += " ";
                    }
                    else if (numbers){
                        new_content += (unsigned char) std::stoi(integer);
                    }
                    else{
                        new_content += integer; //in this case "integer" is actually already a char
                    }
                    integer = "";
                }
                else{
                    integer += c;
                }
            }
            content = new_content;

            iss >> readsIdx;
            string readIdx = "";
            vector<int> readIdxs;
            for (auto c : readsIdx){
                if (c == ','){
                    readIdxs.push_back(std::atoi(readIdx.c_str()));
                    readIdx = "";
                }
                else{
                    readIdx += c;
                }
            }

            Column snp;
            snp.pos = std::atoi(pos.c_str());
            snp.ref_base = ref_base;
            snp.second_base = second_frequent_base;
            int cov = 0;
            for (int n = 0; n < content.size(); n++){
                if (content[n] != ' ' && cov < max_coverage){ //ignore the reads after this
                    snp.content.push_back(content[n]);
                    snp.readIdxs.push_back(readIdxs[n]);
                }   
                if (content[n] != ' ' && readIdxs[n] >= 0){
                    cov++;
                }  
            }
            snps[snps.size()-1].push_back(snp);                

        }
        else if (line_type == "READ"){
            numberOfReads[numberOfReads.size()-1]++;
            names_of_reads[names_of_reads.size()-1].push_back(line);
            string name, startRead, endRead, startContig, endContig;
            iss >> name >> startRead >> endRead >> startContig >> endContig;
            // cout << "LINE sep reasd qslkjmflmsqdj " << line << endl;
            try {
                readLimits[readLimits.size()-1].push_back(make_pair(std::stoi(startContig), std::stoi(endContig)));
            }
            catch (const std::invalid_argument& ia) {
                cout << "error in parsing read limits" << endl;
                cout << "line : " << line << endl;
                exit(1);
            }
        }

    }
    infile.close();
}

/**
 * @brief Compute pairwise similarities and differences between reads
 * 
 * @param snps All the snps in the contig
 * @param sims_and_diffs Resulting matrix of similarities and differences
 */
void list_similarities_and_differences_between_reads(
    std::vector<Column> &snps, 
    vector<vector<pair<int,int>>> &sims_and_diffs){

    // set<int> debug_interesting_reads = {-2, -1};
    // for (auto r = 0 ; r < sims_and_diffs.size() ; r++){
    //     debug_interesting_reads.insert(r);
    // }
    // std::unordered_map<int, string> debug_strings;
    // for (auto r : debug_interesting_reads){
    //     debug_strings[r] = "";
    // }

    for (Column snp : snps){

        cout << "suspssdcect position : " << snp.pos << endl;
        // if (position <200000){
        //     cout<< "ldjflqmj" << endl;
        //     print_snp(snps[position], mask);
        // }

        //what bases occur in which cluster ?
        robin_hood::unordered_map<unsigned char,int> bases_in_total;
        for (auto r = 0 ; r < snp.readIdxs.size() ; r++){
            unsigned char base = snp.content[r];
            if (bases_in_total.find(base) == bases_in_total.end()){
                bases_in_total[base] = 0;
            }
            bases_in_total[base]++;
        }
        //find the second most frequent base
        unsigned char second_most_frequent_base = ' ';
        unsigned char most_frequent_base = ' ';
        int second_most_frequent_base_count = 0;
        int most_frequent_base_count = 0;
        for (auto b : bases_in_total){
            if (b.second > most_frequent_base_count){
                second_most_frequent_base = most_frequent_base;
                second_most_frequent_base_count = most_frequent_base_count;
                most_frequent_base_count = b.second;
                most_frequent_base = b.first;
            }
            else if (b.second > second_most_frequent_base_count){
                second_most_frequent_base = b.first;
                second_most_frequent_base_count = b.second;
            }
        }

        // debug_strings[-2] += second_most_frequent_base;
        // debug_strings[-1] += most_frequent_base;
        

        int idx1 = 0;
        for (int read1 = 0 ; read1 < sims_and_diffs.size() ; read1 ++){
            while (idx1 < snp.readIdxs.size() && snp.readIdxs[idx1] < read1){
                idx1++;
            }
            if (idx1 >= snp.readIdxs.size() || snp.readIdxs[idx1] > read1){
                // if (debug_interesting_reads.find(read1) != debug_interesting_reads.end()){
                //     debug_strings[read1] += " ";
                // }
                continue;
            }
            // else{
            //     if (debug_interesting_reads.find(read1) != debug_interesting_reads.end()){
            //         if (bases_in_total[snps[position].content[idx1]] >= second_most_frequent_base_count){
            //             // debug_strings[read1] += snps[position].content[idx1];
            //             if (snps[position].content[idx1] == most_frequent_base){
            //                 debug_strings[read1] += snps[position].content[idx1];
            //             }
            //             else if (snps[position].content[idx1] == second_most_frequent_base){
            //                 debug_strings[read1] += snps[position].content[idx1];
            //             }
            //             else{
            //                 debug_strings[read1] += " ";
            //             }
            //         }
            //         else{
            //             debug_strings[read1] += " ";
            //         }
            //     }
            // }
            unsigned char base1 = snp.content[idx1];

            int idx2 = idx1+1;
            for (int read2 = read1+1 ; read2 < sims_and_diffs.size() ; read2 ++){
                while (idx2 < snp.readIdxs.size() && snp.readIdxs[idx2] < read2){
                    idx2++;
                }
                if (idx2 >= snp.readIdxs.size() || snp.readIdxs[idx2] > read2){
                    continue;
                }
                unsigned char base2 = snp.content[idx2];

                // if (bases_in_total[base1] > second_most_frequent_base_count && bases_in_total[base2] >= second_most_frequent_base_count){
                //     cout << "posodddiicie ";
                //     print_snp(snps[position], mask);
                //     cout << "baddese1 " << base1 << " " << bases_in_total[base1] << " base2 " << base2 << " " << bases_in_total[base2] << endl;
                // }

                // cout << "qijdioqddsp " << bases_in_total[base1] << " " << second_most_frequent_base_count << " " << bases_in_total[base2] << " " << base1  << " "<< base2 << endl;

                if (bases_in_total[base1] >= second_most_frequent_base_count && bases_in_total[base2] >= second_most_frequent_base_count && base1 != base2){
                    sims_and_diffs[read1][read2].second++;
                    sims_and_diffs[read2][read1].second++;
                }
                else if (bases_in_total[base1] >= second_most_frequent_base_count && bases_in_total[base2] >= second_most_frequent_base_count && base1 == base2){
                    sims_and_diffs[read1][read2].first++;
                    sims_and_diffs[read2][read1].first++;
                    if (base1 == second_most_frequent_base){ //this is very strong signal
                        sims_and_diffs[read1][read2].first+= 2;
                        sims_and_diffs[read2][read1].first+= 2;
                    }
                }
            }
        }
    }
}

/**
 * @brief Compute pairwise similarities and differences between reads
 * 
 * @param snps All the snps in the contig
 * @param sims_and_diffs Resulting matrix of similarities and differences
 */
void list_similarities_and_differences_between_reads2(
    std::vector<Column> &snps, 
    vector<vector<pair<int,int>>> &sims_and_diffs){

    //convert snps to a vector of strings
    vector<string> snps_by_reads (sims_and_diffs.size(), string(snps.size(), ' '));
    vector<char> dominant_char (snps.size(), ' ');
    vector<char> second_char (snps.size(), ' ');
    for (int s = 0 ; s < snps.size() ; s++){
        dominant_char[s] = snps[s].ref_base;
        second_char[s] = snps[s].second_base;
        Column snp = snps[s];
        for (auto r = 0 ; r < snp.readIdxs.size() ; r++){
            snps_by_reads[snp.readIdxs[r]][s] = snp.content[r];
        }
    }

    //now list the similarities and differences
    for (auto r1 = 0 ; r1 < snps_by_reads.size() ; r1++){
        for (auto r2 = r1+1 ; r2 < snps_by_reads.size() ; r2++){
            for (auto pos = 0 ; pos < snps_by_reads[r1].size() ; pos++){
                if (snps_by_reads[r1][pos] != ' ' && snps_by_reads[r2][pos] != ' '){
                    if (snps_by_reads[r1][pos] == snps_by_reads[r2][pos] && snps_by_reads[r1][pos] == dominant_char[pos]){
                        sims_and_diffs[r1][r2].first++;
                        sims_and_diffs[r2][r1].first++;
                    }
                    else if (snps_by_reads[r1][pos] == snps_by_reads[r2][pos] && snps_by_reads[r1][pos] == second_char[pos]){
                        sims_and_diffs[r1][r2].first+= 3;
                        sims_and_diffs[r2][r1].first+= 3;
                    }
                    else if (snps_by_reads[r1][pos] != snps_by_reads[r2][pos] && snps_by_reads[r1][pos] == dominant_char[pos] && snps_by_reads[r2][pos] == second_char[pos]){
                        sims_and_diffs[r1][r2].second++;
                        sims_and_diffs[r2][r1].second++;
                    }
                    else if (snps_by_reads[r1][pos] != snps_by_reads[r2][pos] && snps_by_reads[r1][pos] == second_char[pos] && snps_by_reads[r2][pos] == dominant_char[pos]){
                        sims_and_diffs[r1][r2].second++;
                        sims_and_diffs[r2][r1].second++;
                    }
                }
            }
        }
    }
}

/**
 * @brief Compute pairwise similarities and differences between reads using matrix algebra using the Eigen library
 * 
 * @param snps All the snps in the contig
 * @param similarity sparse matrix numbering the number of identical loci between reads
 * @param difference sparse matrix numbering the number of different loci between reads
 */
void list_similarities_and_differences_between_reads3(
    std::vector<Column> &snps, 
    Eigen::SparseMatrix<int>& similarity,
    Eigen::SparseMatrix<int>& difference){

    //generate two triplet lists
    vector<Eigen::Triplet<int>> triplet_list_ref;
    vector<Eigen::Triplet<int>> triplet_list_alt;
    
    //now fill the matrix with the snps
    for (int s = 0 ; s < snps.size() ; s++){
        Column snp = snps[s];
        for (auto r = 0 ; r < snp.readIdxs.size() ; r++){
            if (snp.content[r] == snp.ref_base){
                triplet_list_ref.push_back(Eigen::Triplet<int>(snp.readIdxs[r], s, 1));
            }
            else if (snp.content[r] == snp.second_base){
                triplet_list_alt.push_back(Eigen::Triplet<int>(snp.readIdxs[r], s, 1));
            }
        }
    }

    //now use two linear algebra to formulate the distance if terms of matrix multiplication
    //build 4 matrices from the subsampled matrix

    //generate matrix_1s from triplet_list_alt
    Eigen::SparseMatrix<int> matrix_1s(similarity.rows(), snps.size());
    matrix_1s.setFromTriplets(triplet_list_alt.begin(), triplet_list_alt.end());

    //second matrix : transpose of matrix_1
    Eigen::SparseMatrix<int> matrix_2s = matrix_1s.transpose();

    //generate matrix_3s from triplet_list_ref
    Eigen::SparseMatrix<int> matrix_3s(similarity.rows(), snps.size()); 
    matrix_3s.setFromTriplets(triplet_list_ref.begin(), triplet_list_ref.end());

    //fourth matrix : transpose of matrix_3
    Eigen::SparseMatrix<int> matrix_4s = matrix_3s.transpose();

    //multiplication of matrices to obtain the similarity and difference matrices
    similarity = 3*matrix_1s * matrix_2s + matrix_3s * matrix_4s; //the 3 is to give more weight to the second most frequent base, this is more exceptionnal
    difference = matrix_1s * matrix_4s + matrix_3s * matrix_2s;

    //set the diagonals to 0
    for (int i = 0 ; i < similarity.outerSize() ; ++i){
        for (Eigen::SparseMatrix<int>::InnerIterator it(similarity, i); it; ++it){
            if (it.row() == it.col()){
                it.valueRef() = 0;
            }
        }
    }

    for (int i = 0 ; i < difference.outerSize() ; ++i){
        for (Eigen::SparseMatrix<int>::InnerIterator it(difference, i); it; ++it){
            if (it.row() == it.col()){
                it.valueRef() = 0;
            }
        }
    }
}

/**
 * @brief Create a graph linking very similar reads
 * 
 * @param mask Indicates which reads span the whole window - process only those
 * @param snps 
 * @param chunk chunk of the contig that is being processed
 * @param sizeOfWindow size of chunks
 * @param adjacency_matrix Result, set of pairs of coordinates (containing only ones)
 * @param errorRate Error rate of the reads
 */
void create_read_graph(
    vector <bool> &mask,
    // std::vector<Column> &snps, 
    int chunk, 
    int sizeOfWindow,
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs,
    vector<vector<int>> &adjacency_matrix,
    float &errorRate){

    set<int> listOfGroups;
    int max_cluster = 0;

    for (int read1 = 0 ; read1 < mask.size() ; read1 ++){
        if (mask[read1]){

            vector <float> distance_with_other_reads (mask.size(), 0);
            int max_compat = 5; //to remove reads that match on few positions, see how much compatibility you can find at most
            for (auto r = 0 ; r < distance_with_other_reads.size() ; r++){
                if (mask[r] && r != read1 && sims_and_diffs[read1][r].first > 0){
                    float diffs = max(0,sims_and_diffs[read1][r].second-1); //-1 to make sure that 1 difference does not make such a difference 
                    distance_with_other_reads[r] = 1 - diffs / float(sims_and_diffs[read1][r].first+sims_and_diffs[read1][r].second);
                    if (sims_and_diffs[read1][r].first > max_compat){
                        max_compat = sims_and_diffs[read1][r].first;
                    }
                }
            }

            for (auto r = 0 ; r < distance_with_other_reads.size() ; r++){
                if (mask[r] && r != read1 && sims_and_diffs[read1][r].first + sims_and_diffs[read1][r].second < max(5.0,0.7*max_compat) ){
                    distance_with_other_reads[r] = 0;
                }
            }


            vector<pair<int, float>> smallest;
            for (int r = 0 ; r < distance_with_other_reads.size() ; r++){
                smallest.push_back(make_pair(r,distance_with_other_reads[r]));
            }
            //sort smallest by distance in decreasing order
            sort(smallest.begin(), smallest.end(), [](const pair<int, float>& a, const pair<int, float>& b) {
                return a.second > b.second;
            });

            int nb_of_neighbors = 0;
            float distance_threshold_below_which_two_reads_are_considered_different = 1 - errorRate*2;
            float distance_threshold_above_which_two_reads_should_be_linked= 1 ;
            if (smallest.size() > 1){
                distance_threshold_above_which_two_reads_should_be_linked = smallest[0].second - (smallest[0].second - smallest[1].second)*3;
            }
            if (distance_threshold_above_which_two_reads_should_be_linked == 1){ //if there are actually identical reads, you still want to tolerate at least some errors, or you will end up making strictly clonal clusters
                //take the first five non-1 distance of smallest
                int idx = 0;
                while (idx < smallest.size() && smallest[idx].second == 1){
                    idx+=1;
                }
                if (idx < smallest.size()){
                    idx = min(idx+4, (int) smallest.size()-1);
                    distance_threshold_above_which_two_reads_should_be_linked = smallest[idx].second;
                }
            }

            for (auto neighbor : smallest){
                if (neighbor.second > distance_threshold_below_which_two_reads_are_considered_different 
                    && (nb_of_neighbors < 5 || neighbor.second == 1 || neighbor.second >= distance_threshold_above_which_two_reads_should_be_linked)
                    && mask[neighbor.first]){
                    nb_of_neighbors++;
                    
                    adjacency_matrix[read1][neighbor.first] = 1;
                    adjacency_matrix[neighbor.first][read1] = 1;
                }
            }
        }
    }

    // cout << "adjacency mawtrix 25 : " << endl;
    // for (auto r = 0 ; r < adjacency_matrix[19].size() ; r++){
    //     cout << r << " - " << adjacency_matrix[19][r] << " " << compatibilities[19][r] << " " << imcompatibilities[19][r] << endl;
    // }
    // cout << endl;

}
 
 
 /**
  * @brief Create a read graph, using less memory but much more time than list_similarities_and_differences and create_read_graph
  * 
  * @param snps 
  * @param mask 
  * @param chunk 
  * @param sizeOfWindow 
  * @param adjacency_matrix 
  * @param errorRate 
  */
void create_read_graph_low_memory(
    std::vector<Column> &snps,
    std::vector <bool> &mask,
    int chunk,
    int sizeOfWindow,
    std::vector<std::vector<int>> &neighbor_list_low_memory, //containing the ordered list of neighbors for each read
    float &errorRate){

    //first convert snps to a sparse matrix
    Eigen::MatrixXi matrix_1(mask.size(), snps.size());
    //fill the matrix with -1
    matrix_1.setConstant(-1);
    int idx_snp;
    for (Column snp : snps){
        for (auto r = 0 ; r < snp.readIdxs.size() ; r++){
            if (snp.content[r] == snp.ref_base){
                matrix_1(snp.readIdxs[r], idx_snp) = 0;
            }
            else if (snp.content[r] == snp.second_base){
                matrix_1(snp.readIdxs[r], idx_snp) = 1;
            }
        }
        idx_snp++;
    }

    //go through the reads and to see which are full of -1 and add them to mask
    vector<bool> mask_extend (mask.size(), false);
    for (auto r = 0 ; r < mask.size() ; r++){
        if (mask[r]){
            bool all_minus_one = true;
            for (auto s = 0 ; s < snps.size() ; s++){
                if (matrix_1(r,s) != -1){
                    all_minus_one = false;
                    break;
                }
            }
            if (all_minus_one == false){
                mask_extend[r] = true;
            }
        }
    }

    for (int read1=0 ; read1 < mask.size() ; read1++){
        if (mask_extend[read1]){
            vector<float> distance_with_other_reads (mask.size(), 0);
            vector<int> similarity_with_other_reads (mask.size(), 0);
            vector<int> difference_with_other_reads (mask.size(), 0); 
            int max_compat = 5;
            for (int read2 = 0 ; read2 < mask.size() ; read2++){
                if (mask_extend[read2] && read1 != read2){
                    //compare the sims and diffs between read1 and read2
                    int nb_similar = 0;
                    int nb_different = 0;
                    //compute the number of cells where read1 and read2 have 1
                    for (int pos = 0 ; pos < snps.size() ; pos++){

                        if (matrix_1(read1, pos) == 1 && matrix_1(read2, pos) == 1){
                            nb_similar+=3;
                        }
                        else if (matrix_1(read1, pos) == 0 && matrix_1(read2, pos) == 0){
                            nb_similar++;
                        }
                        else if (matrix_1(read1, pos) != -1 && matrix_1(read2, pos) != -1){
                            nb_different++;
                        }

                    }

                    distance_with_other_reads[read2] = 1 - max(0, nb_different-1) / float(nb_different+nb_similar); //-1 to make sure that just 1 difference does not make a difference
                    if (nb_similar > max_compat){
                        max_compat = nb_similar;
                    }
                    similarity_with_other_reads[read2] = nb_similar;
                    difference_with_other_reads[read2] = nb_different;
                }
            }

            for (auto r = 0 ; r < distance_with_other_reads.size() ; r++){
                if (mask[r] && r != read1 && similarity_with_other_reads[r] + difference_with_other_reads[r] < max(5.0,0.7*max_compat) ){
                    distance_with_other_reads[r] = 0;
                }
            }
            
            vector<pair<int, float>> smallest;
            for (int r = 0 ; r < distance_with_other_reads.size() ; r++){
                smallest.push_back(make_pair(r,distance_with_other_reads[r]));
            }

            //sort smallest by distance in decreasing order
            sort(smallest.begin(), smallest.end(), [](const pair<int, float>& a, const pair<int, float>& b) {
                return a.second > b.second;
            });

            int nb_of_neighbors = 0;
            float distance_threshold_below_which_two_reads_are_considered_different = 1 - errorRate*2;
            float distance_threshold_above_which_two_reads_should_be_linked= 1 ;
            if (smallest.size() > 1){
                distance_threshold_above_which_two_reads_should_be_linked = smallest[0].second - (smallest[0].second - smallest[1].second)*3;
            }
            if (distance_threshold_above_which_two_reads_should_be_linked == 1){ //if there are actually identical reads, you still want to tolerate at least some errors, or you will end up making strictly clonal clusters
                //take the first five non-1 distance of smallest
                int idx = 0;
                while (idx < smallest.size() && smallest[idx].second == 1){
                    idx+=1;
                }
                if (idx < smallest.size()){
                    idx = min(idx+4, (int) smallest.size()-1);
                    distance_threshold_above_which_two_reads_should_be_linked = smallest[idx].second;
                }
            }

            for (auto neighbor : smallest){
                if (neighbor.second > distance_threshold_below_which_two_reads_are_considered_different 
                    && (nb_of_neighbors < 5 || neighbor.second == 1 || neighbor.second >= distance_threshold_above_which_two_reads_should_be_linked)
                    && mask[neighbor.first]){
                    nb_of_neighbors++;

                    // if (read1 == 11645 || read1 == 12747 || read1 == 22677){ //11645 is the weird cluster, 12747 the perfect one and  22677 the hub
                    //     cout << "adding link " << read1 << " " << neighbor.first << " " << neighbor.second << endl;
                    // }

                    neighbor_list_low_memory[read1].push_back(neighbor.first);
                    neighbor_list_low_memory[neighbor.first].push_back(read1);
                }
            }
        }
    }

    //sort each entry of neighbor_list_low_memory and remove duplicates
    for (auto r = 0 ; r < neighbor_list_low_memory.size() ; r++){
        sort(neighbor_list_low_memory[r].begin(), neighbor_list_low_memory[r].end());
        //remove duplicate using unique
        neighbor_list_low_memory[r] = vector<int>(neighbor_list_low_memory[r].begin(), std::unique(neighbor_list_low_memory[r].begin(), neighbor_list_low_memory[r].end()));
    }
}

/**
 * @brief Create a graph linking very similar reads
 * 
 * @param mask Indicates which reads span the whole window - process only those
 * @param chunk chunk of the contig that is being processed
 * @param sizeOfWindow size of chunks
 * @param similarity input, matrix of similarities
 * @param difference input, matrix of differences
 * @param adjacency_matrix Result, sparse adjacency matrix
 * @param errorRate Error rate of the reads
 */
void create_read_graph_matrix(
    std::vector <bool> &mask,
    int chunk,
    int sizeOfWindow,
    Eigen::SparseMatrix<int>& similarity,
    Eigen::SparseMatrix<int>& difference,
    Eigen::SparseMatrix<int> &adjacency_matrix, //containing only the 1s
    float &errorRate){

    set<int> listOfGroups;
    int max_cluster = 0;

    vector<Eigen::Triplet<int>> triplet_list; //i,j, value to build the sparse adjacnceny matrix efficiently

    for (int read1 = 0 ; read1 < mask.size() ; read1 ++){
        if (mask[read1]){

            vector <float> distance_with_other_reads (mask.size(), 0);
            vector <int> sims (mask.size(), 0);
            vector <int> diffs (mask.size(), 0);
            int max_compat = 5; //to remove reads that match on few positions, see how much compatibility you can find at most

            //iterate through all the other reads

            for (Eigen::SparseMatrix<int>::InnerIterator it(similarity, read1); it; ++it){
                int read2 = it.row();
                if (mask[read2] && read1 != read2){
                    sims[read2] = it.value();
                }
            }

            for (Eigen::SparseMatrix<int>::InnerIterator it(difference, read1); it; ++it){
                int read2 = it.row();
                if (mask[read2] && read1 != read2){
                    diffs[read2] = it.value();
                }
            }

            for (auto r = 0 ; r < mask.size() ; r++){
                if (mask[r] && r != read1 && sims[r]> 0){
                    float diff = max(0,diffs[r]-1); //-1 to make sure that 1 difference does not make such a difference 
                    distance_with_other_reads[r] = 1 - diff / float(sims[r]+diffs[r]);
                    if (sims[r] > max_compat){
                        max_compat = sims[r];
                    }
                }
            }

            for (auto r = 0 ; r < distance_with_other_reads.size() ; r++){
                if (mask[r] && r != read1 && sims[r] + diffs[r] < max(5.0,0.7*max_compat) ){
                    distance_with_other_reads[r] = 0;
                }
            }


            vector<pair<int, float>> smallest;
            for (int r = 0 ; r < distance_with_other_reads.size() ; r++){
                smallest.push_back(make_pair(r,distance_with_other_reads[r]));
            }
            //sort smallest by distance in decreasing order
            sort(smallest.begin(), smallest.end(), [](const pair<int, float>& a, const pair<int, float>& b) {
                return a.second > b.second;
            });

            int nb_of_neighbors = 0;
            float distance_threshold_below_which_two_reads_are_considered_different = 1 - errorRate*2;
            float distance_threshold_above_which_two_reads_should_be_linked= 1 ;
            if (smallest.size() > 1){
                distance_threshold_above_which_two_reads_should_be_linked = smallest[0].second - (smallest[0].second - smallest[1].second)*3;
            }
            if (distance_threshold_above_which_two_reads_should_be_linked == 1){ //if there are actually identical reads, you still want to tolerate at least some errors, or you will end up making strictly clonal clusters
                //take the first five non-1 distance of smallest
                int idx = 0;
                while (idx < smallest.size() && smallest[idx].second == 1){
                    idx+=1;
                }
                if (idx < smallest.size()){
                    idx = min(idx+4, (int) smallest.size()-1);
                    distance_threshold_above_which_two_reads_should_be_linked = smallest[idx].second;
                }
            }

            for (auto neighbor : smallest){
                if (neighbor.second > distance_threshold_below_which_two_reads_are_considered_different 
                    && (nb_of_neighbors < 5 || neighbor.second == 1 || neighbor.second >= distance_threshold_above_which_two_reads_should_be_linked)
                    && mask[neighbor.first]){
                    nb_of_neighbors++;
                    
                    triplet_list.push_back(Eigen::Triplet<int>(read1,neighbor.first,1));
                    triplet_list.push_back(Eigen::Triplet<int>(neighbor.first,read1,1));
                }
            }
        }
    }

    //set matrix from triplet list, do not fill duplicates
    adjacency_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());

    //go through the matrix and set all the positive values to 1
    for (int k = 0; k < adjacency_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency_matrix, k); it; ++it) {
            it.valueRef() = 1;
        }
    }
}



/**
 * @brief Merge the clusters that are too close to each other
 * 
 * @param localClusters initial clusters
 * @param adjacency_matrix graph
 * @param mask reads to be considered
 * @return std::vector<int> resulting clusters
 */
std::vector<int> merge_clusterings(std::vector<std::vector<int>> &localClusters,
    std::vector<std::vector<int>> &neighbor_list, 
    Eigen::SparseMatrix<int> &adjacency_matrix, bool low_memory, std::vector <bool> &mask){

    if (localClusters.size() == 0){
        return vector<int>(neighbor_list.size(), 0);
    }

    vector<double> clusters_aggregated(localClusters[0].size(), 0);
    for (auto i = 0 ; i < localClusters.size() ; i++){
        for (auto j = 0 ; j < localClusters[i].size() ; j++){
            clusters_aggregated[j] += localClusters[i][j]*std::pow(2.0,i);
        }
    }

    std::unordered_map<double, int> clusters_aggregated_map;
    vector<int> cluster_aggregated_ints;
    int index = 0;
    for (auto i = 0 ; i < clusters_aggregated.size() ; i++){
        if (clusters_aggregated_map.find(clusters_aggregated[i]) == clusters_aggregated_map.end()){
            clusters_aggregated_map[clusters_aggregated[i]] = index;
            cluster_aggregated_ints.push_back(index);
            index++;
        }
        else{
            cluster_aggregated_ints.push_back(clusters_aggregated_map[clusters_aggregated[i]]);
        }
    }

    // put -2 for the reads that are not in the mask
    for (auto i = 0 ; i < cluster_aggregated_ints.size() ; i++){
        if (!mask[i]){
            cluster_aggregated_ints[i] = -2;
        }
    }

    vector<int> re_clustered;
    if (low_memory){
        re_clustered = chinese_whispers(neighbor_list, cluster_aggregated_ints, mask);
    }
    else{
        re_clustered = chinese_whispers_high_memory(adjacency_matrix, cluster_aggregated_ints, mask);
    }

    return re_clustered;
}

/**
 * @brief Aggregate many clusterings into one, merging the clusters that are too close to each other
 * 
 * @param snps 
 * @param localClusters 
 * @param strengthened_adjacency_matrix 
 * @param mask_at_this_position 
 * @param haplotypes 
 * @param errorRate 
 */
void finalize_clustering( 
    std::vector<Column> &snps,  
    vector<vector<int>> &localClusters, 
    std::vector<std::vector<int>> &strengthened_neighbor_list,
    Eigen::SparseMatrix<int> &adjacency_matrix,
    bool low_memory, 
    vector<bool> &mask_at_this_position,
    vector<int> &haplotypes,
    float errorRate,
    int posstart,
    int posend){

    if (localClusters.size() == 0){
        for (auto r = 0 ; r < mask_at_this_position.size() ; r++){
            if (mask_at_this_position[r]){
                haplotypes[r] = -1;
            }
            else{
                haplotypes[r] = -2;
            }
        }
        return;
    }
 
    vector<int> new_clusters = merge_clusterings(localClusters, strengthened_neighbor_list, adjacency_matrix, low_memory, mask_at_this_position);
    vector<int> clusteredReads = new_clusters;

    // find all clusters that are too small and flag them as a -1 cluster, so they can be rescued later
    unordered_map<int, int> clusterSizes;
    for (auto r = 0 ; r < clusteredReads.size() ; r++){
        if (mask_at_this_position[r] == false){
            clusteredReads[r] = -2;
        }
        else{
            clusterSizes[clusteredReads[r]] += 1;
        }
    }
    for (auto r = 0 ; r < clusteredReads.size() ; r++){

        float minSizeOfCluster = 5;

        if (clusterSizes[clusteredReads[r]] < minSizeOfCluster && clusteredReads[r] != -2){
            clusteredReads[r] = -1;
        }
    }

    vector<int> mergedHaplotypes = clusteredReads;
    //re-number the clusters
    unordered_map<int, int> clusterToHaplotype;
    int haplotype = 0;
    for (auto r = 0 ; r < mergedHaplotypes.size() ; r++){
        if (mergedHaplotypes[r] > -1){
            if (clusterToHaplotype.find(mergedHaplotypes[r]) == clusterToHaplotype.end()){
                clusterToHaplotype[mergedHaplotypes[r]] = haplotype;
                haplotype++;
            }
            mergedHaplotypes[r] = clusterToHaplotype[mergedHaplotypes[r]];
        }
    }

    //clustered reads represent subsets of reads that have never been separated in no partition
    //however, some haplotypes may have been separated in some partitions (especially if there are many haplotypes)
    //however, they should not be actually separated in snps: merge them

    if (low_memory){
        // cout << "lendidi sep reads " << strengthened_neighbor_list.size() << " " << mergedHaplotypes.size() << endl;
        // for (auto i : mergedHaplotypes){
        //     cout << i << " ";
        // }  
        // cout << endl;
        haplotypes = chinese_whispers(strengthened_neighbor_list, mergedHaplotypes, mask_at_this_position);
    }
    else{
        haplotypes = chinese_whispers_high_memory(adjacency_matrix, mergedHaplotypes, mask_at_this_position);
    }

    //re-number the haplotypes
    unordered_map<int, int> haplotypeToIndex;
    haplotypeToIndex[-1] = -1; 
    if (snps.size() == 0){
        haplotypeToIndex[-1] = 0; 
    }
    haplotypeToIndex[-2] = -2;
    int index_h = 0;
    for (auto h : haplotypes){
        if (haplotypeToIndex.find(h) == haplotypeToIndex.end()){
            haplotypeToIndex[h] = index_h;
            index_h += 1;
        }
    }
    for (auto r = 0 ; r < haplotypes.size() ; r++){
        haplotypes[r] = haplotypeToIndex[haplotypes[r]];
    }

    merge_close_clusters(strengthened_neighbor_list, adjacency_matrix, low_memory, haplotypes, mask_at_this_position); //by disturbing the CW clustering

    haplotypes = merge_wrongly_split_haplotypes(haplotypes, snps, strengthened_neighbor_list, adjacency_matrix, low_memory, posstart, posend); //by looking at the actual reads
}

/**
 * 
 * @brief Checks if the separation of some clusters is never justified by snps and if so, merges them
 * 
 * @param clusteredReads First separation of the reads
 * @param snps Vector of columns containing the snps
 * @param chunk To look only at the interesting snps
 * @param suspectPostitions To look only at the interesting snps
 * @param sizeOfWindow To look only at the interesting snps
 * @return std::vector<int> 
 */

std::vector<int> merge_wrongly_split_haplotypes(
    std::vector<int> &clusteredReads, 
    std::vector<Column> &snps,  
    std::vector<std::vector<int>> &neighbor_list,
    Eigen::SparseMatrix<int> &adjacency_matrix,
    bool low_memory,
    int posstart,
    int posend){

    set<int> listOfGroups;
    int max_cluster = 0;
    unordered_map<int, int> indexOfGroups;
    unordered_map<int,int> sizeOfGroups;

    int index = 0;
    for (auto read = 0 ; read < clusteredReads.size() ; read++){
        if (clusteredReads[read] > -1){
            listOfGroups.emplace(clusteredReads[read]);
            if (indexOfGroups.find(clusteredReads[read]) == indexOfGroups.end()){
                indexOfGroups[clusteredReads[read]] = index;
                index++;
                sizeOfGroups[clusteredReads[read]] = 0;
            }
            sizeOfGroups[clusteredReads[read]] += 1;
        }
    }
    vector<vector<int>> imcompatibilities (listOfGroups.size(), vector<int> (listOfGroups.size(), 0));
    vector<vector<int>> pos_of_last_incompatibilities (listOfGroups.size(), vector<int> (listOfGroups.size(), -10));

    if (listOfGroups.size() <= 1){
        vector <int> one_cluster (clusteredReads.size(), 0);
        for (auto r = 0 ; r < clusteredReads.size() ; r++){
            if (clusteredReads[r] == -2){
                one_cluster[r] = -2;
            }
        }
        return one_cluster;
    }

    //vector associating to each read how many times it is excluded from each cluster
    vector < unordered_map<int, int> > with_which_clusters_goes_this_read (clusteredReads.size());
    vector < int > on_how_many_snps_is_this_read_defined (clusteredReads.size(), 0);
    //fill the vector with 0s
    for (auto read = 0 ; read < clusteredReads.size() ; read++){
        for (auto group : listOfGroups){
            with_which_clusters_goes_this_read[read][group] = 0;
        }
    }

    for (auto snp : snps){
        if (snp.pos >= posstart && snp.pos < posend){

            //what bases occur in which cluster ?
            unordered_map<int,unsigned char> cluster_to_majority_base;
            unordered_map<int,unordered_map<unsigned char,int>> bases_in_each_cluster;
            unordered_map<int,int> nb_bases_in_each_cluster;
            unordered_map<unsigned char,int> bases_in_total;
            for (auto r = 0 ; r < snp.readIdxs.size() ; r++){
                int read = snp.readIdxs[r];
                unsigned char base = snp.content[r];
                int cluster = clusteredReads[read];
                if (cluster > -1){
                    if (bases_in_each_cluster.find(cluster) == bases_in_each_cluster.end()){
                        bases_in_each_cluster[cluster] = unordered_map<unsigned char,int>();
                    }
                    if (bases_in_each_cluster[cluster].find(base) == bases_in_each_cluster[cluster].end()){
                        bases_in_each_cluster[cluster][base] = 0;
                    }
                    if (bases_in_total.find(base) == bases_in_total.end()){
                        bases_in_total[base] = 0;
                    }
                    bases_in_total[base]++;
                    bases_in_each_cluster[cluster][base]++;
                    nb_bases_in_each_cluster[cluster]++;
                }
            }
            
            //count the majority base of each cluster at this position
            set <unsigned char> maxbases;
            for (auto cluster : bases_in_each_cluster){
                int secondMax = 0;
                int max = 0; //at leash 60% the bases need to agree
                char maxBase = ' ';
                for (auto base : cluster.second){
                    if (base.second >= max){
                        maxBase = base.first;
                        secondMax = max;
                        max = base.second;
                    }
                    else if (base.second > secondMax){
                        secondMax = base.second;
                    }
                }
                if (secondMax*2 > max || nb_bases_in_each_cluster[cluster.first]*0.5 > max){
                    maxBase = ' ';
                    // cout << "secondMax: " << secondMax << " max: " << max << endl;
                    // cout << "bases_in_each_cluster[cluster.first] :" << endl;
                    // for (auto base : bases_in_each_cluster[cluster.first]){
                    //     cout << base.first << " " << base.second << endl;
                    // }
                }
                cluster_to_majority_base[cluster.first] = maxBase;
                if (maxBase != ' '){
                    maxbases.emplace(maxBase);
                }
            }
            if (maxbases.size() <= 1){
                continue;
            }

            // if (chunk == 14){
            //     cout << "chqhiuà " << position << endl;
            //     vector<bool> no_mask(clusteredReads.size(), true);
            //     print_snp(snps[position],no_mask);
            // }

            for (auto group1 : listOfGroups){
                for (auto group2 : listOfGroups){
                    if (cluster_to_majority_base[group1] != ' ' && cluster_to_majority_base[group2] != ' ' && group1 > group2){
                        if (cluster_to_majority_base[group1] != cluster_to_majority_base[group2] && snp.pos - pos_of_last_incompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] > 10){ //>10 to make sure it's not two close snps

                            // if (imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] <= 1){
                            //     cout << "fqdsttewlj sep reads are incompatible " << group1 << " " << group2 << " " << snp.pos << " " << "ACGT-"[cluster_to_majority_base[group1]%5]
                            //          << " " << "ACGT-"[cluster_to_majority_base[group2]%5] << endl;
                            // }

                            imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] += 1;
                            imcompatibilities[indexOfGroups[group2]][indexOfGroups[group1]] += 1;

                            pos_of_last_incompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] = snp.pos;
                            pos_of_last_incompatibilities[indexOfGroups[group2]][indexOfGroups[group1]] = snp.pos;

                            // if (group1 == 7 && group2 == 94){
                            //     cout << "incompatibility btw " << group1 << " and " << group2  << "at position " << position << endl;
                            //     cout << cluster_to_majority_base[group1] << " " << cluster_to_majority_base[group2] << endl;
                            //     cout << "bases in each cluster group1 " << endl;
                            //     for (auto base : bases_in_each_cluster[group1]){
                            //         cout << base.first << " " << base.second << endl;
                            //     }
                            //     cout << "bases in each cluster group2 " << endl;
                            //     for (auto base : bases_in_each_cluster[group2]){
                            //         cout << base.first << " " << base.second << endl;
                            //     }
                            //     cout << "eecddxww " << endl;
                            //     int nr = 0;
                            //     for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
                            //         int read = snps[position].readIdxs[r];
                            //         while (nr < read){
                            //             if (clusteredReads[nr] >= 0){
                            //                 cout << " ";
                            //             }
                            //             nr += 1;
                            //         }
                            //         if (clusteredReads[nr] >= 0){
                            //             unsigned char c = snps[position].content[r];
                            //             if (c > 126){
                            //                 cout << (unsigned char) (c - 80);
                            //             }
                            //             else{
                            //                 cout << (unsigned char) c;
                            //             }
                            //         }
                            //         nr += 1;
                            //     }
                            //     cout << endl;
                            // }
                        }
                    }
                }
            }
        }
    }

    // cout << "imcompatibilities computed : " << endl;
    // for (auto group1 : listOfGroups){
    //     for (auto group2 : listOfGroups){
    //         cout << "incompatiebility btw " << group1 << " (" << sizeOfGroups[group1] << ") and " << group2 << " (" << sizeOfGroups[group2] << ") "  << imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] << endl;
    //     }
    //     cout << endl;
    // }

    // list all the distances between different clusters
    std::map < pair<int, int>, double > number_of_links_between_the_two_clusters;
    unordered_map <int,int> number_of_links_in_each_cluster;
    if (low_memory){
        for (auto read1 = 0 ; read1 < clusteredReads.size() ; read1++){
            for (auto read2 = 0 ; read2 < clusteredReads.size() ; read2++){
                bool linked = false;
                if (std::binary_search(neighbor_list[read1].begin(), neighbor_list[read1].end(),read2) == true){
                    linked = true;
                }
                if (linked){
                    int cluster1 = clusteredReads[read1];
                    int cluster2 = clusteredReads[read2];
                    if (cluster1 != cluster2){
                        if (number_of_links_between_the_two_clusters.find(make_pair(cluster1, cluster2)) == number_of_links_between_the_two_clusters.end()){
                            number_of_links_between_the_two_clusters[make_pair(cluster1, cluster2)] = 0;
                        }
                        number_of_links_between_the_two_clusters[make_pair(cluster1, cluster2)]++;
                        
                    }
                    if (number_of_links_in_each_cluster.find(cluster1) == number_of_links_in_each_cluster.end()){
                        number_of_links_in_each_cluster[cluster1] = 0;
                    }
                    number_of_links_in_each_cluster[cluster1]++;
                }
            }
        }
    }
    else{
        for (int k = 0; k < adjacency_matrix.outerSize(); ++k) {
            for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency_matrix, k); it; ++it) {
                int read1 = it.row();
                int read2 = it.col();

                int cluster1 = clusteredReads[read1];
                int cluster2 = clusteredReads[read2];
                if (cluster1 != cluster2){
                    if (number_of_links_between_the_two_clusters.find(make_pair(cluster1, cluster2)) == number_of_links_between_the_two_clusters.end()){
                        number_of_links_between_the_two_clusters[make_pair(cluster1, cluster2)] = 0;
                    }
                    number_of_links_between_the_two_clusters[make_pair(cluster1, cluster2)]++;
                    
                }
                if (number_of_links_in_each_cluster.find(cluster1) == number_of_links_in_each_cluster.end()){
                    number_of_links_in_each_cluster[cluster1] = 0;
                }
                number_of_links_in_each_cluster[cluster1]++;
            }
        }
    }

    for (auto link : number_of_links_between_the_two_clusters){
        // cout << "qfiodududu " << link.second << " " << number_of_links_in_each_cluster[link.first.first] << endl;
        number_of_links_between_the_two_clusters[link.first] = link.second / number_of_links_in_each_cluster[link.first.first];
    }
    // pairs of clusters by closeness
    vector < pair < pair<int, int>, double > > sorted_links;
    for (auto link : number_of_links_between_the_two_clusters){
        sorted_links.emplace_back(link);
    }
    sort(sorted_links.begin(), sorted_links.end(), [](const pair < pair<int, int>, double > &a, const pair < pair<int, int>, double > &b){
        return a.second > b.second;
    });

    //now that all the incompatibilities are computed, we can merge the clusters or not
    unordered_map <int, int> old_group_to_new_group;
    for (auto group : listOfGroups){
        old_group_to_new_group[group] = group;
    }
    old_group_to_new_group[-1] = -1;
    old_group_to_new_group[-2] = -2;

    for (auto pair_of_clusters : sorted_links){
        if (pair_of_clusters.second > 0.01){
            int cluster1 = pair_of_clusters.first.first;
            int cluster2 = pair_of_clusters.first.second;

            if (old_group_to_new_group[cluster1] == old_group_to_new_group[cluster2]){ //have already been merged through another cluster
                continue;
            }
            //check if there are no incompatibilities between the two clusters
            // cout << "looking if i can merge " << cluster1 << " and " << cluster2 << endl;
            bool incompatibility = false;
            for (auto group1 : listOfGroups){
                if (old_group_to_new_group[group1] == old_group_to_new_group[cluster1]){
                    for (auto group2 : listOfGroups){
                        if (old_group_to_new_group[group2] == old_group_to_new_group[cluster2]){
                            if (imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] > 1 ){ //if there is only one incompatible snp, it might be an error
                                // cout << "incompatiebility btw " << group1 << " and " << group2 << " " << imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] << endl;
                                // cout << "merging " << cluster1 << " and " << cluster2 << " would create an incompatibility" << endl;
                                incompatibility = true;
                            }
                        }
                    }
                }
            }
            if (!incompatibility){
                // cout << "merging wiydiu  " << cluster1 << " and " << cluster2 << endl;
                for (auto group2 : listOfGroups){
                    if (old_group_to_new_group[group2] == old_group_to_new_group[cluster2]){
                        old_group_to_new_group[group2] = old_group_to_new_group[cluster1];
                    }
                }
            }
        }
    }

    //re-number the clusters
    unordered_map <int, int> new_group_to_index;
    int new_group_index = 0;
    for (auto group : listOfGroups){
        if (new_group_to_index.find(old_group_to_new_group[group]) == new_group_to_index.end()){
            new_group_to_index[old_group_to_new_group[group]] = new_group_index;
            new_group_index++;
        }
    }
    for (auto group : listOfGroups){
        old_group_to_new_group[group] = new_group_to_index[old_group_to_new_group[group]];
    }

    // cout << "old_group_to_new_group ddrr :" << endl;
    // for (auto group : listOfGroups){
    //     cout << group << " " << old_group_to_new_group[group] << endl;
    // }
    // cout << "new_groups ddrr :" << endl;
    // for (auto new_group : new_groups){
    //     for (auto group : new_group){
    //         cout << group << " ";
    //     }
    //     cout << endl;
    // }

    vector <int> new_clustered_reads (clusteredReads.size(), -1);
    for (auto read = 0 ; read < clusteredReads.size() ; read++){
        new_clustered_reads[read] = old_group_to_new_group[clusteredReads[read]];
    }

    return new_clustered_reads;

}

int main(int argc, char *argv[]){

    if (argc < 8){
        if (argc==2 && (argv[1] == string("-h") || argv[1] == string("--help"))){
            cout << "Usage: ./separate_reads <columns> <num_threads> <error_rate> <meta> <low_memory> <rarest-strain-abundance> <outfile> <DEBUG>" << endl;
            return 0;
        }

        cout << "Usage: ./separate_reads <columns> <num_threads> <error_rate> <meta> <low_memory> <rarest-strain-abundance> <outfile> <DEBUG>" << endl;
        return 1;
    }

    string columns_file = argv[1];
    int num_threads = atoi(argv[2]);
    bool meta = bool(atoi(argv[3]));
    float errorRate = atof(argv[4]);
    bool debug = bool(atoi(argv[8]));
    string outfile = argv[7];
    bool low_memory = bool(atoi(argv[5]));

    float rarest_strain_abundance = atof(argv[6]);
    int max_coverage;
    if (rarest_strain_abundance == 0){
        max_coverage = 1000000000;
    }
    else{
        int max_coverage = 50 / atof(argv[6]); //with a coverage of 50 a strain can be well recovered, thus we can downsample to 50/abundance
    }

    //create empty output file
    ofstream out(outfile);
    out.close();

    std::unordered_map<string, int> index_of_names;
    std::unordered_map<int, string> name_of_contigs;
    vector<vector<string>> names_of_reads;
    std::vector<int> numberOfReads;
    std::vector<std::vector<Column>> snps_in;
    std::vector<long int> length_of_contigs;
    vector<vector<pair<int,int>>> readLimits;
    parse_column_file(columns_file, snps_in, index_of_names, name_of_contigs, names_of_reads, length_of_contigs, readLimits, numberOfReads, max_coverage);

    //to make a progress bar, compute the length of the assembly
    long int total_length = 0;
    for (auto l : length_of_contigs){
        total_length += l;
    }
 
    //choosing size of window: compute the mean length of the reads (and compute coverge)
    int numberOfReadsHere = 0;
    int sumLength = 0;
    int numberOfReadsAbove4000 = 0;
    vector<float> coverages (snps_in.size(), 0);
    int readIdx = 0;
    for (auto c : readLimits){
        for (auto r : c){
            numberOfReadsHere++;
            sumLength += r.second-r.first+1;
            coverages[readIdx] += r.second-r.first+1;
            if (r.second-r.first+1 > 4000){
                numberOfReadsAbove4000++;
            }
        }
        coverages[readIdx] /= length_of_contigs[readIdx];
        readIdx++;
    }
    double meanLength = sumLength / double(numberOfReadsHere);
    int sizeOfWindow = 2000;
    if (numberOfReadsAbove4000 < 20 && meanLength < 4000 && meanLength > 2000){
        sizeOfWindow = 1000;
    }
    else if (numberOfReadsAbove4000 < 20 && meanLength < 2000){
        sizeOfWindow = 500;
    }

    // estimate the haploid coverage as the min coverage of the 90% most covered contigs (weighted by length)
    vector<pair<float,int>> weighted_coverages;
    for (auto i = 0 ; i < coverages.size() ; i++){
        weighted_coverages.push_back(make_pair(coverages[i]*length_of_contigs[i], i));
    }
    sort(weighted_coverages.begin(), weighted_coverages.end(), [](const pair<float,int>& a, const pair<float,int>& b) {
        return a.first > b.first;
    });
    float haploid_coverage = 0;
    float total_length_covered = 0;
    for (auto i = 0 ; i < weighted_coverages.size() ; i++){
        total_length_covered += length_of_contigs[weighted_coverages[i].second];
        if (total_length_covered > 0.9*total_length){
            haploid_coverage = coverages[weighted_coverages[i].second];
            break;
        }
    }
    cout << "HAPOLOID COVERAGE " << haploid_coverage << endl;


    //to make the progress bar
    int total_computed_length = 0; 
    auto time_0 = high_resolution_clock::now();
    auto last_time = time_0;

    //separate the reads on each contig parralelly
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (auto n = 0 ; n < snps_in.size() ; n++){

        // if (name_of_contigs[n].substr(7,9) != "edge_2877"){
        //     cout << "skipping " << name_of_contigs[n] << endl;
        //     continue;
        // }

        bool low_memory_now = low_memory;
        if (coverages[n] > 5000){
            low_memory_now = true;
        }

        auto snps = snps_in[n];
        int numberOfReadsHere = numberOfReads[n];
        if (snps.size() == 0){
            continue;
        }

        #pragma omp critical (cout)
        {
            time_t now = time(NULL);
            cout << "separating reads on contig " << name_of_contigs[n] << " [" << std::ctime(&now) <<"]\n";
        }

        vector<vector<pair<int,int>>> sims_and_diffs;
        
        Eigen::SparseMatrix<int> similarity(numberOfReadsHere, numberOfReadsHere);
        Eigen::SparseMatrix<int> difference(numberOfReadsHere, numberOfReadsHere);
        if (!low_memory_now){
            list_similarities_and_differences_between_reads3(snps, similarity, difference);
        }

        // cout << "similrities and differences computed " << numberOfReadsHere << endl;
        // for (auto i : sims_and_diffs[0]){
        //     cout << i.first << " " << i.second << endl;
        // } 

        vector<pair<pair<int,int>, vector<int>>> threadedReads; //associates to each window the reads that are in it
        int suspectPostitionIdx = 0;
        int chunk = -1;
        int upperBound;
        while ((chunk+1)*sizeOfWindow + 100 <= length_of_contigs[n]){
            // if (chunk*sizeOfWindow != 7000){
            //     cout << "csksdlk " << endl;
            //     chunk++;
            //     continue;
            // }
            chunk++;
            upperBound = (chunk+1)*sizeOfWindow;
            if ((chunk+1)*sizeOfWindow + 100 > length_of_contigs[n]){ //to avoid having extremely short terminal windows
                upperBound = length_of_contigs[n]+1;
            }

            if (suspectPostitionIdx >= snps.size() || snps[suspectPostitionIdx].pos > upperBound-1){
                //no snp in this window, just add the reads, with -2 for the reads that are not here and 0 for the reads that are here
                vector<int> readsHere (numberOfReadsHere, -2);
                for (auto r = 0 ; r < readLimits[n].size() ; r++){
                    int pointLeft = chunk*sizeOfWindow;
                    int pointRight = min(upperBound-1, int(length_of_contigs[n]));
                    int middlePoint = (pointLeft+pointRight)/2;
                    if (middlePoint < 500){
                        middlePoint = min(500, int(length_of_contigs[n]/2));
                    }
                    if (middlePoint > int(length_of_contigs[n])-500){
                        middlePoint = max(int(length_of_contigs[n]/2), int(length_of_contigs[n])-500);
                    }

                    if (readLimits[n][r].first <= middlePoint && readLimits[n][r].second >= middlePoint){
                        readsHere[r] = 0;
                    }
                }
                
                threadedReads.push_back(make_pair(make_pair(chunk*sizeOfWindow, min(upperBound-1, int(length_of_contigs[n]))), readsHere));
                continue;
            }

            //only consider reads that span the whole length of the window: create a mask
            vector<bool> mask_at_this_position (numberOfReadsHere, false); 
            for (auto r = 0 ; r < snps[suspectPostitionIdx].readIdxs.size() ; r++){
                mask_at_this_position[snps[suspectPostitionIdx].readIdxs[r]] = true;
            }

            while (suspectPostitionIdx < snps.size() && snps[suspectPostitionIdx].pos < upperBound-1){
                suspectPostitionIdx++;
            } //suspectPostitionIdx is now the first position after the window
            if (suspectPostitionIdx > 0){
                suspectPostitionIdx--; //suspectPostitionIdx is now the last position in the window
            }
            int idxmask = 0;
            for (auto r = 0 ; r < snps[suspectPostitionIdx].readIdxs.size() ; r++){
                while (idxmask < snps[suspectPostitionIdx].readIdxs[r]){
                    mask_at_this_position[idxmask] = false;
                    idxmask++;
                }
                idxmask++;
            }
            suspectPostitionIdx++; //suspectPostitionIdx is now the first position after the window

            vector<vector<int>> neighbor_list_low_memory (numberOfReadsHere, vector<int> (0));
            vector<vector<int>> neighbor_list_low_memory_strengthened (numberOfReadsHere, vector<int> (0));
            vector<vector<int>> adjacency_matrix_high_memory;
            vector<vector<int>> strengthened_adjacency_matrix_high_memory;
            Eigen::SparseMatrix<int> adjacency_matrix (numberOfReadsHere, numberOfReadsHere);

            if (!low_memory_now){
                create_read_graph_matrix(mask_at_this_position, chunk, sizeOfWindow, similarity, difference, adjacency_matrix, errorRate);
            }
            else{
                for (auto v : neighbor_list_low_memory){
                    v.reserve(20);
                }
                // cout << "he qd q lll" << endl;
                create_read_graph_low_memory(snps, mask_at_this_position, chunk, sizeOfWindow, neighbor_list_low_memory, errorRate);

                // cout << "ididinnd q" << endl;
                // neighbor_list_low_memory_strengthened = strengthen_adjacency_matrix(neighbor_list_low_memory, numberOfReadsHere);
                neighbor_list_low_memory_strengthened = neighbor_list_low_memory;
            }
            // cout << "ociojccood" << endl;

            vector<int> clustersStart (numberOfReadsHere, 0);
            for (auto r = 0 ; r < numberOfReadsHere ; r++){
                clustersStart[r] = r;
            }

            vector<vector<int>> allclusters_debug;
            vector<int> clusteredReads1;
            vector<int> clusteredReads1_2;
            if (low_memory_now){
                clusteredReads1 = chinese_whispers(neighbor_list_low_memory_strengthened, clustersStart, mask_at_this_position);
            }
            else{
                clusteredReads1 = chinese_whispers_high_memory(adjacency_matrix, clustersStart, mask_at_this_position);
            }
            allclusters_debug.push_back(clusteredReads1);
            vector<vector<int>> localClusters = {};

            // cout << "here are all the interesting positions" << endl;
            // for (auto p : interestingPositions){
            //     cout << p << " ";
            // }
            // cout << endl;
            // cout << "runnignd cw again and again" << endl;
            int lastpos = -10;
            for (auto snp : snps){
                if (snp.pos >= chunk*sizeOfWindow && snp.pos < chunk*sizeOfWindow + sizeOfWindow && snp.pos > lastpos+10){
                    lastpos = snp.pos;
                    // cout << "in dldjk position " << position << " : " << endl;

                    unordered_map<unsigned char, int> charToIndex;
                    vector<int> clustersStart2 (numberOfReadsHere, 0);
                    for (auto r = 0 ; r < numberOfReadsHere ; r++){
                        clustersStart2[r] = r;
                    }
                    for (auto r = 0 ; r < snp.content.size() ; r++){
                        if (mask_at_this_position[snp.readIdxs[r]] == true){
                            if (charToIndex.find(snp.content[r]) == charToIndex.end()){
                                charToIndex[snp.content[r]] = snp.readIdxs[r];
                            }
                            int indexHere = charToIndex[snp.content[r]];
                            clustersStart2[snp.readIdxs[r]] = indexHere;
                        }
                    }
                    
                    vector<int> clusteredReads_local;// = chinese_whispers(strengthened_adjacency_matrix, clustersStart2, mask_at_this_position); 
                    if (low_memory_now){
                        clusteredReads_local = chinese_whispers(neighbor_list_low_memory_strengthened, clustersStart2, mask_at_this_position);
                    }
                    else{
                        clusteredReads_local = chinese_whispers_high_memory(adjacency_matrix, clustersStart2, mask_at_this_position);
                    }
                    localClusters.push_back(clusteredReads_local);
                    
                    allclusters_debug.push_back(clustersStart2);
                    allclusters_debug.push_back(clusteredReads_local);
                }         
            }

            vector<int> haplotypes(numberOfReadsHere, -2);
            finalize_clustering(snps, localClusters, neighbor_list_low_memory_strengthened, adjacency_matrix, low_memory, mask_at_this_position, haplotypes, errorRate, chunk*sizeOfWindow, chunk*sizeOfWindow + sizeOfWindow);
            // cout << "outputting graph hs/tmp/graph_" <<  std::to_string(chunk*sizeOfWindow) +".gdf" << endl;
            // if (low_memory){
            //     outputGraph_low_memory(neighbor_list_low_memory_strengthened, haplotypes, "hs/tmp/graph_"+std::to_string(chunk*sizeOfWindow)+".gdf");
            // }
            // else{
            //     outputGraph(adjacency_matrix_high_memory, haplotypes, "hs/tmp/graph_"+std::to_string(chunk*sizeOfWindow)+".gdf");
            //     vector<bool> nomask (haplotypes.size(), true);
            //     outputGraph_several_clusterings(adjacency_matrix_high_memory, allclusters_debug, nomask, "hs/tmp/graph_"+std::to_string(chunk*sizeOfWindow)+"_all.gdf");
            // }
            // exit(0);


            threadedReads.push_back(make_pair(make_pair(chunk*sizeOfWindow, min(upperBound-1, int(length_of_contigs[n]))), haplotypes));
        }
        //recursively go through all the windows once again, seeing if a window can be inspired by its neighbors

        //progress bar if time since last progress bar > 10 seconds
        if (duration_cast<seconds>(high_resolution_clock::now() - last_time).count() > 10){
            last_time = high_resolution_clock::now();
            auto total_time_elapsed = duration_cast<seconds>(last_time - time_0).count();
            #pragma omp critical (cout)
            {
                if (length_of_contigs[n] < 0){
                    cout << "NEGATIVE LENGTH OF CONTIG " << name_of_contigs[n] << " " << length_of_contigs[n] << endl;
                }
                else{
                    total_computed_length += length_of_contigs[n];
                    cout << "Progress in separate_reads: " << 100*total_computed_length/total_length << "%. Estimated time left: " << total_time_elapsed*total_length/(total_computed_length) - total_time_elapsed << " seconds. \n";
                }
            }
        }

        //if the assembly is not meta, check that all haplotypes are abundant enough to correspond to a chromosome. Suppress the ones that are not
        if (!meta){
            for (auto r : threadedReads){
                //count each haplotype in threadedReads.second
                unordered_map<int, int> haplotypeToCount;
                for (auto h : r.second){
                    if (h >= 0){
                        if (haplotypeToCount.find(h) == haplotypeToCount.end()){
                            haplotypeToCount[h] = 0;
                        }
                        haplotypeToCount[h]++;
                    }
                }

                double coverage_limit = min(haploid_coverage*0.7, 25.0); //above this, the haplotype is trustworthy

                //count how many haplotypes are above the coverage limit
                int numberOfHaplotypesAboveLimit = 0;
                for (auto h : haplotypeToCount){
                    if (h.second > coverage_limit){
                        numberOfHaplotypesAboveLimit++;
                    }
                    else{
                        cout << "Dismissing haplotype " << h.first << " on contig " << name_of_contigs[n] << " between " << r.first.first << " and " << r.first.second << " because it has a coverage of " << h.second << " < " << coverage_limit << endl;
                    }
                }

                for (auto h = 0 ; h < r.second.size() ; h++){
                    if (haplotypeToCount[r.second[h]] < coverage_limit){
                        r.second[h] = -1;
                    }
                }

                if (numberOfHaplotypesAboveLimit <= 1){
                    for (auto h = 0 ; h < r.second.size() ; h++){
                        if (r.second[h] != -2){
                            r.second[h] = 0;
                        }
                    }
                }
            }
        }

        
        //append threadedReads to file
        #pragma omp critical
        {
            ofstream out(outfile, std::ios_base::app);
            out << name_of_contigs[n] << endl;
            for (auto r : names_of_reads[n]){
                out << r << "\n";
            }
            //now output all the groups
            vector<int> readsHere;
            vector<int> groups;
            for (auto r : threadedReads){
                readsHere.clear();
                groups.clear();
                out << "GROUP\t" << r.first.first << "\t" << r.first.second << "\t";
                for (auto h = 0 ; h < r.second.size() ; h++){
                    int group = r.second[h];
                    if (group != -2){
                        readsHere.push_back(h);
                        groups.push_back(group);
                    }
                }
                for (auto h = 0 ; h < readsHere.size() ; h++){
                    out << readsHere[h] << ",";
                }
                out << "\t";
                for (auto h = 0 ; h < groups.size() ; h++){
                    out << groups[h] << ",";
                }
                out << "\n";
            }
            out.close();
        }
    }

    return 0;
}










