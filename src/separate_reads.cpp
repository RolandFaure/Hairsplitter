#include "separate_reads.h"
#include "cluster_graph.h"

#include <iostream>
#include <fstream> //for reading files
#include <sstream> //for parsing strings
#include <string>
#include <omp.h> //for efficient parallelization
#include <map>
#include <set>
#include <cmath> //for sqrt and pow

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
    std::vector<int>& numberOfReads){

    std::ifstream infile(file);
    std::string line;
    int nbreadshere = 0;
    names_of_reads = vector<vector<string>>(0);
    bool numbers = false;
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
                ref_base = ref_base_string[0];//the string should be of length 1 anyway
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
                        new_content += (char) std::stoi(integer);
                    }
                    else{
                        new_content += integer; //in this cas "integer" is actually already a char
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
            for (int n = 0; n < content.size(); n++){
                if (content[n] != ' '){
                    snp.content.push_back(content[n]);
                    snp.readIdxs.push_back(readIdxs[n]);
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

        // cout << "suspssdcect position : " << position << endl;
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
            if (snp.readIdxs[idx1] > read1){
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
 * @brief Create a graph linking very similar reads
 * 
 * @param mask Indicates which reads span the whole window - process only those
 * @param snps 
 * @param chunk chunk of the contig that is being processed
 * @param sizeOfWindow size of chunks
 * @param adjacency_matrix Result
 * @param errorRate Error rate of the reads
 */
void create_read_graph(
    vector <bool> &mask,
    std::vector<Column> &snps, 
    int chunk, 
    int sizeOfWindow,
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs,
    std::vector< std::vector<int>> &adjacency_matrix,
    float &errorRate){

    set<int> listOfGroups;
    int max_cluster = 0;
    unordered_map<int, int> indexOfGroups;

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

            // if (read1 == 16 || true){
            //     cout << "ddqfhhe distance of 16 to other reads: " << endl;
            //     for (int r = 0 ; r < distance_with_other_reads.size() ; r++){
            //         if (mask[r] && r != read1){
            //             cout << r << " " << distance_with_other_reads[r] << " " << float(sims_and_diffs[read1][r].second) << " " << float(sims_and_diffs[read1][r].first) << endl;
            //         }
            //     }
            //     cout << endl;
            //     // cout << "sims and diffffsss between 133 and 43 : " << sims_and_diffs[read1][43].first << " " << sims_and_diffs[read1][43].second << " " << distance_with_other_reads[43]<< endl;
            //     // cout << "sims and diffffsss between 133 and 141 : " << sims_and_diffs[read1][141].first << " " << sims_and_diffs[read1][141].second << " " << distance_with_other_reads[141] << endl;
            // }

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
            
            // if (read1 == 1){
            //     cout << "readss1 : " << read1 << endl;
            //     for (auto r = 0 ; r < distance_with_other_reads.size() ; r ++){
            //         if (mask[r]){
            //             cout << r << " " << distance_with_other_reads[r] << ", ";
            //         }
            //     }
            //     cout << endl;
            // }

            int nb_of_neighbors = 0;
            float distance_threshold_below_which_two_reads_are_considered_different = 1 - errorRate*5;
            float distance_threshold_above_which_two_reads_should_be_linked= 1;
            if (smallest.size() > 1){
                distance_threshold_above_which_two_reads_should_be_linked = smallest[0].second - (smallest[0].second - smallest[1].second)*3;
            }

            for (auto neighbor : smallest){
                if (neighbor.second > distance_threshold_below_which_two_reads_are_considered_different 
                    && (nb_of_neighbors < 5 || neighbor.second == 1 || neighbor.second > distance_threshold_above_which_two_reads_should_be_linked)
                    && mask[neighbor.first]){
                    nb_of_neighbors++;
                    
                    adjacency_matrix[read1][neighbor.first] = 1;
                    adjacency_matrix[neighbor.first][read1] = 1;
                }
            }

            // for (auto read2 = 0 ; read2 < mask.size() ; read2 ++){
            //     if (mask[read2] && distance_with_other_reads[read2] >= threshold){
            //         adjacency_matrix[read1][read2] = 1;
            //         // adjacency_matrix[read2][read1] = 1;
            //     }
            // }
        }
    }

    // cout << "adjacency mawtrix 25 : " << endl;
    // for (auto r = 0 ; r < adjacency_matrix[19].size() ; r++){
    //     cout << r << " - " << adjacency_matrix[19][r] << " " << compatibilities[19][r] << " " << imcompatibilities[19][r] << endl;
    // }
    // cout << endl;

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
    std::vector< std::vector<int>> &adjacency_matrix, std::vector <bool> &mask){

    if (localClusters.size() == 0){
        return vector<int>(adjacency_matrix.size(), 0);
    }

    vector<double> clusters_aggregated(localClusters[0].size(), 0);
    for (auto i = 0 ; i < localClusters.size() ; i++){
        for (auto j = 0 ; j < localClusters[i].size() ; j++){
            clusters_aggregated[j] += localClusters[i][j]*std::pow(2.0,i);
        }
    }

    unordered_map<double, int> clusters_aggregated_map;
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

    auto re_clustered = chinese_whispers(adjacency_matrix, cluster_aggregated_ints, mask);

    return re_clustered;
}

int main(int argc, char *argv[]){
    if (argc < 5){
        cout << "Usage: ./separate_reads <columns> <num_threads> <error_rate> <DEBUG> <outfile> " << endl;
        
        return 1;
    }

    string columns_file = argv[1];
    int num_threads = atoi(argv[2]);
    float errorRate = atof(argv[3]);
    bool debug = bool(atoi(argv[4]));
    string outfile = argv[5];

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
    parse_column_file(columns_file, snps_in, index_of_names, name_of_contigs, names_of_reads, length_of_contigs, readLimits, numberOfReads);

    //choosing size of window: compute the mean length of the 1000 first reads
    int numberOfReadsHere = 0;
    int sumLength = 0;
    for (auto c : readLimits){
        for (auto r : c){
            numberOfReadsHere++;
            sumLength += r.second-r.first+1;
            if (numberOfReadsHere > 1000){
                break;
            }
        }
        if (numberOfReadsHere > 1000){
            break;
        }
    }
    double meanLength = sumLength / double(numberOfReadsHere);
    int sizeOfWindow = 2000;
    if (meanLength < 4000 && meanLength > 2000){
        sizeOfWindow = 1000;
    }
    else if (meanLength < 2000){
        sizeOfWindow = 500;
    }

    //separate the reads on each contig parralelly
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (auto n = 0 ; n < snps_in.size() ; n++){
        auto snps = snps_in[n];
        int numberOfReadsHere = numberOfReads[n];
        if (snps.size() == 0){
            continue;
        }

        #pragma omp critical (cout)
        {
            cout << "separating reads on contig " << name_of_contigs[n] << "\r";
        }

        vector<vector<pair<int,int>>> sims_and_diffs (numberOfReadsHere, vector<pair<int,int>>(numberOfReadsHere, make_pair(0,0)));

        list_similarities_and_differences_between_reads(snps, sims_and_diffs);

        // cout << "similrities and differences computed between reads 0 and other" << endl;
        // for (auto i : sims_and_diffs[0]){
        //     cout << i.first << " " << i.second << endl;
        // }
 

        vector<pair<pair<int,int>, vector<int>>> threadedReads;
        int suspectPostitionIdx = 0;
        int chunk = -1;
        while ((chunk+1)*sizeOfWindow-1 < length_of_contigs[n]){
            chunk++;

            if (suspectPostitionIdx >= snps.size() || snps[suspectPostitionIdx].pos > (chunk+1)*sizeOfWindow){
                //no snp in this window, just add the reads, with -2 for the reads that are not here and 0 for the reads that are here
                vector<int> readsHere (numberOfReadsHere, -2);
                for (auto r = 0 ; r < readLimits[n].size() ; r++){
                    int pointLeft = chunk*sizeOfWindow;
                    int pointRight = min((chunk+1)*sizeOfWindow-1, int(length_of_contigs[n]));
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
                
                threadedReads.push_back(make_pair(make_pair(chunk*sizeOfWindow, min((chunk+1)*sizeOfWindow-1, int(length_of_contigs[n]))), readsHere));
                continue;
            }
            //only consider reads that span the whole length of the window: create a mask
            vector<bool> mask_at_this_position (numberOfReadsHere, false); 
            for (auto r = 0 ; r < snps[suspectPostitionIdx].readIdxs.size() ; r++){
                mask_at_this_position[snps[suspectPostitionIdx].readIdxs[r]] = true;
            }

            while (suspectPostitionIdx < snps.size() && snps[suspectPostitionIdx].pos < (chunk+1)*sizeOfWindow){
                suspectPostitionIdx++;
            } //suspectPostitionIdx is now the first position after the window
            suspectPostitionIdx--; //suspectPostitionIdx is now the last position in the window

            int idxmask = 0;
            for (auto r = 0 ; r < snps[suspectPostitionIdx].readIdxs.size() ; r++){
                while (idxmask < snps[suspectPostitionIdx].readIdxs[r]){
                    mask_at_this_position[idxmask] = false;
                    idxmask++;
                }
                idxmask++;
            }

            suspectPostitionIdx++; //suspectPostitionIdx is now the first position after the window

            vector<vector<int>> adjacency_matrix (mask_at_this_position.size(), vector<int>(mask_at_this_position.size(), 0));
            create_read_graph(mask_at_this_position, snps, chunk, sizeOfWindow, sims_and_diffs, adjacency_matrix, errorRate);
            vector<int> clustersStart (adjacency_matrix.size(), 0);
            for (auto r = 0 ; r < adjacency_matrix.size() ; r++){
                clustersStart[r] = r;
            }
            auto strengthened_adjacency_matrix = strengthen_adjacency_matrix(adjacency_matrix);
            strengthened_adjacency_matrix = adjacency_matrix;

            vector<vector<int>> allclusters_debug;
            vector<int> clusteredReads1 = chinese_whispers(strengthened_adjacency_matrix, clustersStart, mask_at_this_position);
            allclusters_debug.push_back(clusteredReads1);
            vector<vector<int>> localClusters = {};

            // cout << "here are all the interesting positions" << endl;
            // for (auto p : interestingPositions){
            //     cout << p << " ";
            // }
            // cout << endl;

            for (auto snp : snps){
                if (snp.pos >= chunk*sizeOfWindow && snp.pos < chunk*sizeOfWindow + sizeOfWindow){
                    // cout << "in dldjk position " << position << " : " << endl;

                    unordered_map<unsigned char, int> charToIndex;
                    vector<int> clustersStart2 (adjacency_matrix.size(), 0);
                    for (auto r = 0 ; r < adjacency_matrix.size() ; r++){
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
                    
                    vector<int> clusteredReads_local = chinese_whispers(strengthened_adjacency_matrix, clustersStart2, mask_at_this_position); 
                    localClusters.push_back(clusteredReads_local);
                    
                    allclusters_debug.push_back(clustersStart2);
                    allclusters_debug.push_back(clusteredReads_local);
                }         
            }

            vector<int> new_clusters = merge_clusterings(localClusters, strengthened_adjacency_matrix, mask_at_this_position);
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

                float minSizeOfCluster = min(float(5.0),max(float(3.0), errorRate*100)); //with HiFi reads, we can find clusters of size 3

                if (clusterSizes[clusteredReads[r]] < minSizeOfCluster && clusteredReads[r] != -2){
                    clusteredReads[r] = -1;
                }
            }

            // vector<int> mergedHaplotypes = merge_wrongly_split_haplotypes(clusteredReads, snps, chunk, interestingPositions, adjacency_matrix, sizeOfWindow);
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

            // cout << "merged haploitypes : " << endl;
            // for (auto h : mergedHaplotypes){
            //     if (h > -1){
            //         cout << h << " ";
            //     }
            // }
            // cout << endl;

            //clustered reads represent subsets of reads that have never been separated in no partition
            //however, some haplotypes may have been separated in some partitions (especially if there are many haplotypes)
            //however, they should not be actually separated in snps: merge them

            vector<int> haplotypes = chinese_whispers(strengthened_adjacency_matrix, mergedHaplotypes, mask_at_this_position);

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

            merge_close_clusters(strengthened_adjacency_matrix, haplotypes, mask_at_this_position);
            if (debug){
                cout << "haploutypes sepreads : " << chunk << endl;
                int n = 0;
                for (auto h : haplotypes){
                    if (h > -1){
                        // cout << n << ":" <<h << " ";
                        cout << h;
                    }
                    // else{
                    //     cout << " ";
                    // }
                    n += 1;
                }
                cout << endl;
            }
            threadedReads.push_back(make_pair(make_pair(chunk*sizeOfWindow, min((chunk+1)*sizeOfWindow-1, int(length_of_contigs[n]))), haplotypes));
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










