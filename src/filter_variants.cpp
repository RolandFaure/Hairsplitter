#include "filter_variants.h"
#include "Partition.h"
#include "tools.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h> //for efficient parallelization
#include <cmath>
#include <set>
#include <cctype> //to check the type of a string (is it a number ?)

// #include "clipp.h" //library to build command line interfaces

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
using std::min;
using std::max;
using std::ofstream;
using std::ifstream;
using std::pair;
using std::make_pair;
// using namespace clipp;

/**
 * @brief Parse a column file and store the SNPs for each contig
 * 
 * @param file File to parse
 * @param snps Associates each contig with its SNPs
 * @param name_of_contigs Associates each contig with its index in snps
 * @param numberOfReads Associates each contig with the number of reads aligned on it
 */
void parse_column_file(
    std::string file, 
    std::vector<std::vector<Column>> &snps, 
    std::unordered_map<string, int>& name_of_contigs, 
    std::unordered_map<int, std::string>& name_of_contigs2,
    std::vector<int> &numberOfReads){

    std::ifstream infile(file);
    std::string line;
    bool numbers = false; //are the variants encoded as numbers or as letters ?
    bool firstsnpline = true;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        std::string line_type;
        iss >> line_type;
        if (line_type == "CONTIG"){
            name_of_contigs[line] = snps.size();
            string nameofcontig;
            iss >> nameofcontig;
            name_of_contigs2[snps.size()] = nameofcontig;
            snps.push_back(std::vector<Column>());
        }
        else if (line_type == "SNPS"){

            string pos;
            string ref_base_string;
            string second_frequent_base_string;
            char ref_base;
            char second_frequent_base;
            string idxs_str;
            string content;
            iss >> pos >> ref_base_string >> second_frequent_base_string >> idxs_str >> content;

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

            //content is a comma-separated list of the bases at this position represented by ints: remove the commas and convert to char
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

            //parse idx as a comma-separated list of integers
            vector<int> idxs;
            string idx = "";
            for (auto c : idxs_str){
                if (c == ','){
                    idxs.push_back(std::stoi(idx));
                    idx = "";
                }
                else{
                    idx += c;
                }
            }

            Column snp;
            snp.pos = std::atoi(pos.c_str());
            snp.ref_base = ref_base;
            snp.second_base = second_frequent_base;
            for (int n = 0; n < content.size(); n++){
                if (content[n] != ' '){
                    snp.content.push_back(content[n]);
                    snp.readIdxs.push_back(idxs[n]);
                }                
            }
            snps[snps.size()-1].push_back(snp);               
        }

    }
    infile.close();
}

/**
 * @brief Output the robust columns in a new column file
 * 
 * @param initial_column_file Use the original column file to get the contigs and the reads
 * @param new_snps 
 */
void output_new_column_file(std::string initial_column_file, std::vector<std::vector<Column>>  &new_snps, std::unordered_map<string, int> &name_of_contigs, string &output_file){

    ofstream out(output_file);
    std::string contig;
    std::ifstream infile(initial_column_file);
    std::string line;
    bool firstsnpline = false;
    while (std::getline(infile, line)){
        std::istringstream iss(line);
        std::string line_type;
        iss >> line_type;
        if (line_type == "CONTIG"){
            firstsnpline = false;
            contig = line;
            out << line << "\n";
        }
        else if (line_type == "READ"){
            out << line << "\n";
        }
        else if (line_type == "SNPS" && !firstsnpline){
            firstsnpline = true;
            //then output all the new snps
            for (auto snp : new_snps[name_of_contigs[contig]]){
                out << "SNPS\t" << snp.pos << "\t" << (int) snp.ref_base << "\t" << (int) snp.second_base << "\t";
                string content = "";
                string idxs = "";
                for (auto c  = 0 ; c < snp.content.size(); c++){
                    idxs += std::to_string(snp.readIdxs[c]) + ",";
                    content += std::to_string( (int) snp.content[c]) + ",";
                }
                out << idxs << "\t" << content << "\n";
            }
        }
    }
    infile.close();
    out.close();

}

/**
 * @brief Output the robust columns in a new vcf file
 * 
 * @param new_snps 
 * @param name_of_contigs 
 * @param output_file 
 */
void output_new_vcf_file(
    std::string &vcf_in,
    std::vector<std::vector<Column>>  &new_snps, 
    std::unordered_map<int, std::string>& name_of_contigs, 
    std::string &output_file){

    //create a set inventorying all the positions that are in the new snps
    std::unordered_map<string, std::set<int>> new_snp_positions;
    for (auto contig = 0 ; contig < new_snps.size() ; contig++){
        new_snp_positions[name_of_contigs[contig]] = std::set<int>();
        for (auto snp : new_snps[contig]){
            new_snp_positions[name_of_contigs[contig]].insert(snp.pos);
        }
    }

    //open the vcf file and write the header
    ofstream out(output_file);
    ifstream in(vcf_in);

    //go through all the lines and output only the ones which positions are in the new snps
    std::string line;
    while (std::getline(in, line)){
        std::istringstream iss(line);
        
        if (line[0] == '#'){
            out << line << "\n";
        }
        else { //this is a putative variant
            std::string contigname;
            iss >> contigname;
            int pos;
            iss >> pos;
            if (new_snp_positions.find(contigname) != new_snp_positions.end() && new_snp_positions[contigname].find(pos) != new_snp_positions[contigname].end()){
                out << line << "\n";
            }
        }
    }

    in.close();
    out.close();
}

/**
 * @brief Compute the distance from the partition to the column
 * 
 * @param par1 
 * @param par2 
 * @param ref corresponding reference base
 * @return distancePartition (n00, n01, n10, n11, phased, partition_to_augment)
 */
distancePartition distance(Partition &par1, Column &par2, char ref_base){

    /*
    when computing the distance, there is not 5 letters but 2 : the two alleles, which are the two most frequent letters
    */
    distancePartition res;
    res.augmented = true;
    vector <int> idxs1 = par1.getReads();
    vector<short> part1 = par1.getPartition();
    vector<float> confs1 = par1.getConfidence();
    vector<int> more1 = par1.getMore();
    vector<int> less1 = par1.getLess();

    vector <unsigned int> idxs2 = par2.readIdxs;
    vector <unsigned char> part2 = par2.content;
    
    float numberOfBases = 0;

    unordered_map<unsigned char, int> content2;

    auto n2 = 0;
    auto n1 = 0;
    for (auto r : idxs2){
        while (n1 < idxs1.size() && idxs1[n1] < idxs2[n2]){
            n1++;
        }
        if (n1 >= idxs1.size()){
            break;
        }
        if ( idxs1[n1] == idxs2[n2] && part1[n1] != -2){
            if (content2.find(part2[n2]) == content2.end()){
                content2[part2[n2]] = 0;
            }
            numberOfBases+=1;
            content2[part2[n2]] += 1;
        }
        n2++;
    }

    if (numberOfBases == 0){ //not comparable
        res.n00 = 0;
        res.n01 = 0;
        res.n10 = 0;
        res.n11 = 0;
        res.solid10 = 0;
        res.solid11 = 0;
        res.solid00 = 0;
        res.solid01 = 0;
        res.augmented = false;
        return res;
    }

    //determine first and second most frequent bases in par2
    
    unsigned char mostFrequent = ref_base;
    auto maxFrequence = content2[ref_base];
    //now find the second most frequent base
    unsigned char secondFrequent = ' ';
    int maxFrequence2 = -1;
    for (auto c : content2){
        if (ref_base != c.first) {
            if (c.second > maxFrequence2){
                secondFrequent = c.first;
                maxFrequence2 = c.second;
            }
        }
    }
    // cout << "mqodjkd " << mostFrequent << " " << secondFrequent << "\n";

    int matches00 = 0;
    int matches01 = 0;
    int matches10 = 0;
    int matches11 = 0;
    int solid11 = 0;
    int solid10 = 0;
    int solid01 = 0;
    int solid00 = 0;
    Column newPartition;

    newPartition.readIdxs = {};
    for (auto ci = 0 ; ci < par2.content.size() ; ci++){
        unsigned char c = par2.content[ci];
        if (c == mostFrequent){
            newPartition.readIdxs.push_back(par2.readIdxs[ci]);
            newPartition.content.push_back('A');

        }
        else if (c == secondFrequent){
            newPartition.readIdxs.push_back(par2.readIdxs[ci]);
            newPartition.content.push_back('a');
        }
        else{
            newPartition.readIdxs.push_back(par2.readIdxs[ci]);
            newPartition.content.push_back(' ');
        }
    }

    // cout << " idx1ss : " << endl;
    // for (auto i = 0 ; i < idxs1.size() ; i++){
    //     if (part1[i] != -2){
    //         cout << idxs1[i] << " ";
    //     }
    // }
    // cout << "\n";
    // cout << " idx2ss : " << endl;
    // for (auto i : par2.readIdxs){
    //     cout << i << " ";
    // }
    // cout << "\n";

    string debug1 = "";
    string debug2 = "";
    n1 = 0;
    auto it1 = idxs1.begin();
    n2 = 0;
    auto it2 = par2.readIdxs.begin();
    while (it1 != idxs1.end() && it2 != par2.readIdxs.end()){

        if (*it1 == *it2){
            ++it1;
            ++it2;

            if (par2.content[n2] == mostFrequent){

                if (part1[n1] == 1){
                    matches11 += 1;
                    debug1 += "1";
                    debug2 += par2.content[n2];
                    if (less1[n1] <= 1 && more1[n1] >= 3){
                        solid11 += 1;
                    }
                }
                else if (part1[n1] == -1){
                    matches01 += 1;
                    debug1 += "0";
                    debug2 += par2.content[n2];
                    if (less1[n1] <= 1 && more1[n1] >= 3){
                        solid01 += 1;
                    }
                }
            }
            else if (par2.content[n2] == secondFrequent){

                if (part1[n1] == 1){
                    matches10 += 1;
                    debug1 += "1";
                    debug2 += par2.content[n2];
                    if (less1[n1] <= 1 && more1[n1] >= 3){
                        solid10 += 1;
                    }
                }
                else if (part1[n1] == -1){
                    matches00 += 1;
                    debug1 += "0";
                    debug2 += par2.content[n2];
                    if (less1[n1] <= 1 && more1[n1] >= 3){
                        solid00 += 1;
                    }
                }
            }
            // cout << "nt1 " << n1 << " " << part1[n1] << " n2 " << n2 << " " << part2[n2] << " "<< *it1 <<  "\n";
            n1++;
            n2++;
        }
        else if (*it2 > *it1){
            ++it1;
            n1++;
        }
        else{
            ++it2;
            n2++;
        }
    }

    res.n00 = matches00;
    res.n01 = matches01;
    res.n10 = matches10;
    res.n11 = matches11;
    res.solid10 = solid10;
    res.solid11 = solid11;
    res.solid00 = solid00;
    res.solid01 = solid01;
    res.phased = 1;
    res.secondBase = secondFrequent;
    res.partition_to_augment = newPartition;
    //cout << "Computing..." << maxScore << " " << par1.size()-res.nonComparable << " " << res.nmismatch << endl;

    // cout << "fqiopudpc \n" << debug1 << "\n" << debug2 << "\n";

    return res;
}

/**
 * @brief comparing how close two partitions are
 * 
 * @param par1 
 * @param par2 
 * @param threshold_p number of highly diverging reads to consider the partitions different
 * @return distancePartition 
 */
distancePartition distance(Partition &par1, Partition &par2, int threshold_p){
    /*
    Two metrics are used to compare two partition : the chi to see if the two partition correlate when they are small
                                                    the p when the partitions are bigger and you can do probabilities on them
    */
    int numberOfComparableBases = 0;

    float chi = 0;

    vector<int> idx1 = par1.getReads();
    vector<int> idx2 = par2.getReads();

    vector<short> part1 = par1.getPartition();
    vector<short> part2 = par2.getPartition();

    vector<int> more1 = par1.getMore();
    vector<int> less1 = par1.getLess();

    vector<int> more2 = par2.getMore();
    vector<int> less2 = par2.getLess();

    int scores [2] = {0,0}; //the scores when directing mostFrequent on either mostfrequent2 or secondFrequent2
    short ndivergentPositions[2] = {0,0}; //number of positions where these two partitions could not have been so different by chance
    short nNotSurePositions[2] = {0,0}; //number of positions where the partitions disagree even though one of them was pretty sure
    //remember all types of matches for the chi square test
    int matches00[2] = {0,0};
    int matches01[2] = {0,0};
    int matches10[2] = {0,0};
    int matches11[2] = {0,0};

    int r1 = 0;
    int r2 = 0;
    while (r1 < idx1.size() && r2 < idx2.size()){
        if (idx1[r1]< idx2[r2]){
            r1++;
        }
        else if (idx2[r2] < idx1[r1]){
            r2++;
        }
        else if (more1[r1] > 1 && more2[r2] > 1){

            numberOfComparableBases += 1;

            float threshold1 = 0.5*(more1[r1]+less1[r1]) + 3*sqrt((more1[r1]+less1[r1])*0.5*(1-0.5)); //to check if we deviate significantly from the "random read", that is half of the time in each partition
            float threshold2 = 0.5*(more2[r2]+less2[r2]) + 3*sqrt((more2[r2]+less2[r2])*0.5*(1-0.5)); //to check if we deviate significantly from the "random read", that is half of the time in each partition

            if (part2[r2] == 1){

                if (part1[r1] == 1){
                    scores[0] += 1;
                    scores[1] -= 1;
                    matches11[0] += 1;
                    matches10[1] += 1;

                    //if both positions are certain, this may be bad
                    if (more1[r1] > threshold1 && more2[r2] > threshold2){
                        ndivergentPositions[1] += 1;
                    }
                    if (more1[r1] > threshold1 || more2[r2] > threshold2){
                        nNotSurePositions[1] += 1;
                    }
                }
                else if (part1[r1] == -1){
                    scores[0] -= 1;
                    scores[1] += 1;
                    matches01[0] += 1;
                    matches00[1] += 1;

                    //if both positions are certain, this may be bad
                    if (more1[r1] > threshold1 && more2[r2] > threshold2){
                        ndivergentPositions[0] += 1;
                    }
                    if (more1[r1] > threshold1 || more2[r2] > threshold2){
                        nNotSurePositions[0] += 1;
                    }
                }
            }
            else if (part2[r2] == -1){

                if (part1[r1] == 1){
                    scores[0] -= 1;
                    scores[1] += 1;
                    matches10[0] += 1;
                    matches11[1] += 1;

                    //if both positions are certain, this may be bad
                    if (more1[r1] > threshold1 && more2[r2] > threshold2){
                        ndivergentPositions[0] += 1;
                    }
                    if (more1[r1] > threshold1 || more2[r2] > threshold2){
                        nNotSurePositions[0] += 1;
                    }
                }
                else if (part1[r1] == -1){
                    scores[0] += 1;
                    scores[1] -= 1;
                    matches00[0] += 1;
                    matches01[1] += 1;

                    //if both positions are certain, this may be bad
                    if (more1[r1] > threshold1 && more2[r2] > threshold2){
                        ndivergentPositions[1] += 1;
                    }
                    if (more1[r1] > threshold1 || more2[r2] > threshold2){
                        nNotSurePositions[1] += 1;
                    }
                }
            }

            r1++;
            r2++;
        }
        else{
            r1++;
            r2++;
        }
    }


    distancePartition res;
    res.augmented = true;

    //check if there are too many unenxplainable positions
    if ((ndivergentPositions[0] >= threshold_p && ndivergentPositions[1] >= threshold_p) || (nNotSurePositions[0] >= 5 && nNotSurePositions[1] >= 5) || numberOfComparableBases == 0){
        // cout << "Should not merge those two partitions ! " << endl;
        res.augmented = false;
    }

    /*
    now there aren't too many unexplainable postions. 
    However, that could be due to little partitions on which we could not do stats. For those, do a chi-square
    */

    //now look at the best scores

    int maxScore = scores[0];
    int maxScoreIdx = 0; //can be either 0 or 1
    
    if (scores[1] > maxScore){
        maxScore = scores[1];
        maxScoreIdx = 1;
    }

    res.n00 = matches00[maxScoreIdx];
    res.n01 = matches01[maxScoreIdx];
    res.n10 = matches10[maxScoreIdx];
    res.n11 = matches11[maxScoreIdx];
    res.phased = -2*maxScoreIdx + 1; // worth -1 or 1

    return res ;
}

/**
 * @brief Compute the chi-square test with one degree of freedom
 * 
 * @param dis 
 * @return float 
 */
float computeChiSquare(distancePartition dis){

    int n = dis.n00 + dis.n01 + dis.n10 + dis.n11;
    if (n == 0){
        return 0;
    }
    float pmax1 = float(dis.n10+dis.n11)/n;
    float pmax2 = float(dis.n01+dis.n11)/n;
    //now a few exceptions when chisquare cannot be computed
    if (pmax1*(1-pmax1) == 0 && pmax2*(1-pmax2) == 0){
        return -1;
    }
    if (pmax1*pmax2*(1-pmax1)*(1-pmax2) == 0){ //if there is only one base in one partition, it can't be compared
        return 0;
    }

    // cout << "ps : " << pmax1 << ", " << pmax2 << endl;
    // cout << "expected/obtained : " << dis.n00 << "/"<<(1-pmax1)*(1-pmax2)*n << " ; " << dis.n01 << "/"<<(1-pmax1)*pmax2*n
    // << " ; " << dis.n10 << "/" << pmax1*(1-pmax2)*n << " ; " << dis.n11 << "/" << pmax1*pmax2*n << endl;
    //chi square test with 1 degree of freedom
    float res;
    res = pow((dis.n00-(1-pmax1)*(1-pmax2)*n),2)/((1-pmax1)*(1-pmax2)*n)
        + pow((dis.n01-(1-pmax1)*pmax2*n),2)/((1-pmax1)*pmax2*n)
        + pow((dis.n10-pmax1*(1-pmax2)*n),2)/(pmax1*(1-pmax2)*n)
        + pow((dis.n11-pmax1*pmax2*n),2)/(pmax1*pmax2*n);

    return res;
}

/**
 * @brief Keep only robust variants, i.e. variants that occur on several different positions in the contig
 * 
 * @param snps_in Input snps
 * @param snps_out Filtered snps
 * @param num_threads 
 */
void keep_only_robust_variants(
    std::vector<std::vector<Column>>  &snps_in, 
    std::vector<std::vector<Column>>  &snps_out, 
    std::unordered_map<string, int> &name_of_contigs,
    float mean_error,
    int num_threads){

    snps_out = vector<vector<Column>>(snps_in.size());
    omp_set_num_threads(num_threads);
    //iterate parallely over all contigs
    #pragma omp parallel for
    for (auto n = 0 ; n < snps_in.size() ; n++){
        // cout << "filtering contqoflmj ig " << n << endl;
        auto snps = snps_in[n];
        auto snp_out = vector<Column>();

        vector<Partition> partitions;
        int lastposition = -5;

        //iterate over all SNPs in the contig
        for (auto snp : snps){
            //go through the partitions to see if this suspicious position looks like smt we've seen before
            if (snp.pos - lastposition <= 5){
                continue;
            }
            bool found = false;
            auto position = snp.pos;
            for (auto p = 0 ; p < partitions.size() ; p++){
                //if the partition is too far away, do not bother comparing
                if (std::abs(snp.pos-partitions[p].get_right())>50000){
                    continue;
                }
                distancePartition dis = distance(partitions[p], snp, snp.ref_base);
                auto comparable = dis.n00 + dis.n11 + + dis.n01 + dis.n10;

                //if ((float(dis.n01+dis.n10)/(min(dis.n00,dis.n11)+dis.n01+dis.n10) <= meanDistance*2 || dis.n01+dis.n10 <= 2)  && dis.augmented && comparable > min(10.0, 0.3*numberOfReads)){
                if (dis.n01 < 0.1 * (dis.n00+dis.n01) && dis.n10 < 0.1 * (dis.n11+dis.n10) && comparable >= snp.readIdxs.size()/2){
                
                    found = true;
                    // if (p == 94){
                    //     cout << "quaomugj emja aaa " << position << " " << endl;
                    //     // partitions[p].print();
                    //     // for (auto r : dis.partition_to_augment.content){
                    //     //     cout << r;
                    //     // }
                    //     // cout << endl;
                    //     int idx = 0;
                    //     for (auto r : debug_reads_of_interest){
                    //         while(idx < dis.partition_to_augment.readIdxs.size() && dis.partition_to_augment.readIdxs[idx] < r){
                    //             idx++;
                    //         }
                    //         if (idx < dis.partition_to_augment.readIdxs.size() && dis.partition_to_augment.readIdxs[idx] == r){
                    //             cout << dis.partition_to_augment.content[idx];
                    //         }
                    //         else{
                    //             cout << " ";
                    //         }
                    //     }
                    //     cout << endl;
                    // }
                    partitions[p].augmentPartition(dis.partition_to_augment, position);

                    break;
                }
            }
            if (!found){    // the second condition is here to create partitions only at specific spots of the backbone
                partitions.push_back(Partition(snp, position, snp.ref_base));
            }
            else{
                lastposition = snp.pos;      //two suspect positions next to each other can be artificially correlated through alignement artefacts
            }
        }

        if (partitions.size() == 0){ //there are no position of interest
            continue;
        }
        
        float threshold = min(4, max(2, (int) (mean_error*100)));

        vector<Partition> listOfFinalPartitions;
        for (auto p1 = 0 ; p1 < partitions.size() ; p1++){

            // if (partitions[p1].number() > 5){
            //     cout << "iqdoudofq non informative partition : "  << p1 << " " << snps.size()<< endl;
            //     partitions[p1].print();
            // }
            
            if (partitions[p1].number() > threshold && partitions[p1].isInformative(false, mean_error)){

                bool different = true;
                
                for (auto p2 = 0 ; p2 < listOfFinalPartitions.size() ; p2++){

                    distancePartition dis = distance(listOfFinalPartitions[p2], partitions[p1], 2);

                    if (dis.augmented 
                        && (dis.n00+dis.n11 > 5*(dis.n01+dis.n10) || dis.n10 + dis.n01 > 5*(dis.n00+dis.n11))
                        && dis.n10 < max(2,2*dis.n01) && dis.n01 < max(2,2*dis.n10)){
                        Partition newPart = listOfFinalPartitions[p2];

                        newPart.mergePartition(partitions[p1], dis.phased);

                        // newPart.print();
                        // for (auto r = 0 ; r < newPart.getPartition().size() ; r++){
                        //     if (newPart.getPartition()[r] == 1){// newPart.getReads()[r] == 496){
                        //         cout << newPart.getReads()[r] << " " << newPart.getPartition()[r] << ":" << newPart.getMore()[r] << ":" << newPart.getLess()[r] << " ";
                        //     }
                        // }
                        // cout << endl;
                        
                        // exit(1);

                        //see if confidence is improved by merging the two partitions, meaning differences were shaky
                        if (dis.n01+dis.n10 < 0.1*(dis.n00+dis.n11) || newPart.compute_conf() > listOfFinalPartitions[p2].compute_conf()){

                            listOfFinalPartitions[p2].mergePartition(partitions[p1], dis.phased);
                            different = false;
                            break;
                        }
                    }
                }
                
                if (different){
                    listOfFinalPartitions.push_back(partitions[p1]);
                }
            }
        }

        // print the final partitions
        // cout << "final partitions : " << endl;
        // for (auto p = 0 ; p < listOfFinalPartitions.size() ; p++){
        //     cout << "partition " << p << " : " << endl;
        //     listOfFinalPartitions[p].print();
        // }
        // exit(1);

        //list the interesting positions
        for (auto snp : snps){
            // cout << "suspecct possisssiion : " << endl;
            // print_snp(snp);
            for (auto p = 0 ; p < listOfFinalPartitions.size() ; p++){
                distancePartition dis = distance(listOfFinalPartitions[p], snp, snp.ref_base);
                if (dis.n00 + dis.n01 + dis.n10 + dis.n11 > 0.5*snp.content.size() 
                    && computeChiSquare(dis) > 5){
                                        
                    snps_out[n].push_back(snp);
                    break;
                }
            }
        }
        // cout << "number of intesdresting positions : " << snps_out[n].size() << endl;
    }
}

/**
 * @brief Generates the MSA of all reads against a backbone
 * 
 * @param bbcontig Backbone read
 * @param allOverlaps All the overlaps of the input reads 
 * @param allreads All the input reads
 * @param snps Result of the function: a vector of Column, each column corresponding to one position on the MSA
 * @param backboneReadIndex Numerotation of the backbone read
 * @param readLimits Limits of the reads on the backbone. Used to recompute coverage of the backbone
 * @param newref Sequence of backbone read in the MSA space
 * @param tmpFolder Folder where to store temporary files
 * @param DEBUG If true, print debug information
 * @return The mean distance between the aligned reads and the consensus backbone
 */
float generate_msa(
    long int bbcontig, 
    std::vector <Overlap> &allOverlaps, 
    std::vector <Read> &allreads, 
    std::vector<Column> &snps, 
    robin_hood::unordered_map<int, int> &insertionPos,
    int backboneReadIndex, 
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits, 
    std::string &newref,
    std::string &tmpFolder,
    bool DEBUG){

    string ACGT = "ACGT-";

    //go through the neighbors of the backbone read and align it

    //keep count of the distance between two reads to know the mean distance
    float totalDistance = 0;
    double totalLengthOfAlignment = 1; //1 to avoid dividing by 0

    //to remember on what backbone the backbone read is leaning -> itself
    allreads[bbcontig].new_backbone(std::make_pair(backboneReadIndex, allreads[bbcontig].neighbors_.size()), allreads[bbcontig].neighbors_.size()+1);
    string read_str = allreads[bbcontig].sequence_.str();

    //small loop to compute truePartition DEBUG
    // if (DEBUG){
    //     for (auto n = 0 ; n<allreads[bbcontig].neighbors_.size() ; n++){

    //         long int neighbor = allreads[bbcontig].neighbors_[n];
    //         Overlap overlap = allOverlaps[neighbor];
            
    //         if (overlap.sequence1 == bbcontig){

    //             truePar.push_back(allreads[overlap.sequence2].name[1]);
    //             // cout << "name : " << allreads[overlap.sequence2].name << " " << allreads[overlap.sequence2].name[1] << " " << truePartition[truePartition.size()-1] << endl;
    //         }
    //         else{
    //             truePar.push_back(allreads[overlap.sequence1].name[1]);
    //             // cout << "name : " << allreads[overlap.sequence1].name << " " << allreads[overlap.sequence1].name[1] << " " << truePartition[truePartition.size()-1] << endl;
    //         }    
    //     }
    // }

    // string fileOut = "/home/rfaure/Documents/these/overlap_filtering/species/Escherichia/triploid/answer_"+std::to_string(read)+".tsv";
    // std::ofstream out(fileOut);
    // for (auto c : truePart){out<<c;}out.close();
    // cout << "name : " << allreads[read].name << " " << allreads[read].name[1] << " " << truePartition[truePartition.size()-1] << endl;
    // /* compute only true partition
    //now mark down on which backbone read those reads are leaning
    //in the same loop, inventoriate all the polishing reads
    vector <string> polishingReads;
    vector <std::pair <int,int>> positionOfReads; //position of polishing reads on the consensus
    vector <string> CIGARs;
    for (auto n = 0 ; n < allreads[bbcontig].neighbors_.size() ; n++){
        long int neighbor = allreads[bbcontig].neighbors_[n];
        Overlap overlap = allOverlaps[neighbor];
        
        // cout << "overlap : " << allreads[overlap.sequence1].name << " on " <<allreads[overlap.sequence2].name << ", " << overlap.position_1_1 << " " << overlap.position_1_2 << " " << overlap.position_2_1 << " "
        //     << overlap.position_2_2 << endl;

        if (overlap.CIGAR != ""){
            CIGARs.push_back(overlap.CIGAR);
            if (overlap.sequence1 == bbcontig){
                allreads[overlap.sequence2].new_backbone(make_pair(backboneReadIndex,n), allreads[bbcontig].neighbors_.size()+1);
                if (overlap.strand){
                    polishingReads.push_back(allreads[overlap.sequence2].sequence_.str());
                    positionOfReads.push_back(make_pair(overlap.position_1_1, overlap.position_1_2));
                }
                else{
                    polishingReads.push_back(allreads[overlap.sequence2].sequence_.reverse_complement().str());
                    positionOfReads.push_back(make_pair(overlap.position_1_1, overlap.position_1_2));
                }
            }
            else { //generally go through here if SAM file is used
                //cout << "iiop Neighbor of " << allreads[overlap.sequence2].name << " : " << allreads[overlap.sequence1].name << endl;
                allreads[overlap.sequence1].new_backbone(make_pair(backboneReadIndex,n), allreads[bbcontig].neighbors_.size()+1);
                if (overlap.strand){
                    polishingReads.push_back(allreads[overlap.sequence1].sequence_.str());
                    positionOfReads.push_back(make_pair(overlap.position_2_1, overlap.position_2_2));
                }
                else{
                    polishingReads.push_back(allreads[overlap.sequence1].sequence_.reverse_complement().str());
                    positionOfReads.push_back(make_pair(overlap.position_2_1, overlap.position_2_2));
                }
            }
        }
        else{

            if (overlap.sequence1 == bbcontig){
                // cout << "Neighbor of " << allreads[overlap.sequence1].name << " : " << allreads[overlap.sequence2].name << endl;
                allreads[overlap.sequence2].new_backbone(make_pair(backboneReadIndex,n), allreads[bbcontig].neighbors_.size()+1);
                if (overlap.strand){
                    polishingReads.push_back(allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1, overlap.position_2_2-overlap.position_2_1).str());
                    positionOfReads.push_back(make_pair(overlap.position_1_1, overlap.position_1_2));
                }
                else{
                    polishingReads.push_back(allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1, overlap.position_2_2-overlap.position_2_1).reverse_complement().str());
                    positionOfReads.push_back(make_pair(overlap.position_1_1, overlap.position_1_2));
                }
            }
            else { //generally go through here if PAF file is used
                // if (polishingReads.size() == 2999 || true){
                //     cout << "iiop2 Neighbor of " << allreads[overlap.sequence1].name << " : " << allreads[overlap.sequence2].name << endl;
                //     cout << "positions of overlap: " << overlap.strand << " " << overlap.position_1_1 << " " << overlap.position_1_2 << " " << overlap.position_2_1 << " " << overlap.position_2_2 << endl;
                // }
                allreads[overlap.sequence1].new_backbone(make_pair(backboneReadIndex,n), allreads[bbcontig].neighbors_.size()+1);
                if (overlap.strand){
                    polishingReads.push_back(allreads[overlap.sequence1].sequence_.subseq(overlap.position_1_1, overlap.position_1_2-overlap.position_1_1).str());
                    positionOfReads.push_back(make_pair(overlap.position_2_1, overlap.position_2_2));
                }
                else{
                    polishingReads.push_back(allreads[overlap.sequence1].sequence_.subseq(overlap.position_1_1, overlap.position_1_2-overlap.position_1_1).reverse_complement().str());
                    positionOfReads.push_back(make_pair(overlap.position_2_1, overlap.position_2_2));
                }
            }
        }
    }

    string consensus = read_str; //if the input assembly is already polished

    snps = vector<Column>(consensus.size());
    for (auto c = 0 ; c < consensus.size() ; c++){
        snps[c].pos = c;
    }

    //while all the alignments are computed, build the positions
    vector<int> numberOfInsertionsHere (consensus.size()+1, 0);

    float alignmentTime = 0;
    float MSAtime = 0;

    for (auto n = 0 ; n < polishingReads.size() ; n++){

        if (DEBUG && n%10 == 0){
            cout << "Aligned " << n << " reads out of " << allreads[bbcontig].neighbors_.size() << " on the backbone\r" << std::flush;
        }

        string alignment;
        if (CIGARs.size() != polishingReads.size()){ //compute the exact alignment if it is not already computed
            cout << "ERROR: CIGARs not well computed" << endl;
            exit(1);
        }
        else{
            alignment = convert_cigar(CIGARs[n]);
            // if (n == 29){
            //     cout << "Alignkldjmfq kdidid " << alignment << " " << allreads[allOverlaps[allreads[bbcontig].neighbors_[n]].sequence1].name.substr(0,10) << " " << endl;
            // }
        }
        // cout << alignment << " ";
        //cout << alignment.size() << " " << std::count(alignment.begin(), alignment.end(), 'I') << " " << std::count(alignment.begin(), alignment.end(), 'D') << "\n";

        totalLengthOfAlignment += alignment.size();
        // mappingQuality.push_back(float(-aligner.getAlignmentScore())/alignment.size());

        // if (n == 10) {break;}

        // //a loop going through the CIGAR and modifyning snps
        int indexQuery = positionOfReads[n].first; //query corresponds to consensus
        int indexTarget = 0; //target corresponds to the read
        int numberOfInsertionsThere = 0;
        int maxNumberOfInsertions = 20; //maximum number of insertions allowed in a row, then they are discarded

        int lengthOfWindow = 100; //number of bases to consider for the local divergence
        int localDivergence = 0; //number of mismatches/insertions/deletions in the last 100 bases
        vector<bool> divergences (lengthOfWindow, false); //vector of booleans to keep track of the last 100 bases

        // if (indexQuery < 20){
        //     cout << "beginning of query : " << indexQuery << " " << polishingReads[n].substr(indexQuery,100) << endl;
        // }
        //cout the alignment
        // cout << "alignment : " << alignment << endl;
        // cout << "ccxoaa ";
        // for (auto de = 0 ; de < min(indexQuery,20000) ; de+= 100){
        //     cout << " ";
        // }

        bool solidAlignment = false;
        int size_of_deletion = 0; //to avoid excluding homopolymers in long deletions
        
        unsigned char unaligned_pos = '!'+25*6;
        char previous_previous_previous_char = 'A';
        char previous_previous_char = 'C';
        char previous_char = 'G';
        for (int l = 0; l < alignment.size(); l++) {

            if (indexQuery < consensus.size()){

                if (l >= lengthOfWindow){
                    if (divergences[l%lengthOfWindow]){
                        localDivergence -= 1;
                        divergences[l%lengthOfWindow] = false;
                    }
                }

                if (alignment[l] == '=' || alignment[l] == 'X' || alignment[l] == 'M'){
                    //fill inserted columns with '-' just before that position
                    // numberOfConsecutiveMatches += 1;
                    // for (int ins = numberOfInsertionsThere ; ins < numberOfInsertionsHere[indexQuery] ; ins++){ //in these positions, insert '-' instead of ' '
                    //     snps[insertionPos[10000*indexQuery+ins]].readIdxs.push_back(n);
                    //     snps[insertionPos[10000*indexQuery+ins]].content.push_back('-');
                    // }

                    previous_previous_previous_char = previous_previous_char;
                    previous_previous_char = previous_char;
                    previous_char = polishingReads[n][indexTarget];

                    unsigned char three_mer = '!' + 5*ACGT.find(previous_previous_previous_char) + ACGT.find(previous_previous_char) + 25*ACGT.find(previous_char);
                    snps[indexQuery].readIdxs.push_back(n);
                    snps[indexQuery].content.push_back(three_mer);

                    // if (indexQuery == 393){
                    //     cout << "addinggg " << three_mer << " " << previous_previous_previous_char << " " << previous_previous_char << " " << previous_char << " " << endl;
                    // }

                    // if ((n == 22 || n == 109 || n == 4 || n == 5 || n == 131) && indexQuery < 40){
                    //     cout << polishingReads[n][indexTarget];
                    //     if (indexQuery == 28){
                    //         cout << " " << previous_previous_previous_char << previous_previous_char << previous_char << " ";
                    //     }
                    // }

                    indexQuery++;
                    indexTarget++;
                    numberOfInsertionsThere = 0;

                    if (alignment[l] != 'M' && alignment[l] != '='){
                        localDivergence += 1;
                        divergences[l%lengthOfWindow] = true;
                    }
                    size_of_deletion = 0;
                }
                else if (alignment[l] == 'S' || alignment[l] == 'H'){
                    indexTarget++;
                    localDivergence += 1;
                    divergences[l%lengthOfWindow] = true;
                }
                else if (alignment[l] == 'D'){
                    // //fill inserted columns with '-' just before that position
                    // for (int ins = numberOfInsertionsThere ; ins < numberOfInsertionsHere[indexQuery] ; ins++){ //in these positions, insert '-' instead of ' '
                    //     snps[insertionPos[10000*indexQuery+ins]].readIdxs.push_back(n);
                    //     snps[insertionPos[10000*indexQuery+ins]].content.push_back('-');
                    // }

                    size_of_deletion += 1;

                    previous_previous_previous_char = previous_previous_char;
                    previous_previous_char = previous_char;
                    previous_char = '-';

                    unsigned char three_mer = '!' + 5*ACGT.find(previous_previous_previous_char) + ACGT.find(previous_previous_char) + 25*ACGT.find(previous_char);

                    snps[indexQuery].readIdxs.push_back(n);
                    snps[indexQuery].content.push_back(three_mer);                    
                    indexQuery++;
                    numberOfInsertionsThere = 0;

                    // if (indexQuery == 29){
                    //     cout << "musdhg: " << previous_previous_previous_char << " " << previous_previous_char << " " << previous_char << " " << three_mer << " " << endl;
                    // }
                    // if ((n == 22 || n == 109 || n == 4 || n == 5 || n == 131) && indexQuery < 40){
                    //     cout << "-";
                    //     if (indexQuery == 28){
                    //         cout << " ";
                    //     }
                    // }

                    totalDistance += 1;

                    localDivergence += 1;
                    divergences[l%lengthOfWindow] = true;
                }
                else if (alignment[l] == 'I'){ 
                    // if (numberOfInsertionsHere[indexQuery] <= maxNumberOfInsertions && indexQuery > positionOfReads[n].first) {

                    //     if (numberOfInsertionsThere >= numberOfInsertionsHere[indexQuery]) { //i.e. this is a new column
                    //         insertionPos[10000*indexQuery+numberOfInsertionsHere[indexQuery]] = snps.size();
                    //         numberOfInsertionsHere[indexQuery] += 1;
                            
                    //         Column newInsertedPos;
                    //         newInsertedPos.readIdxs = snps[indexQuery-1].readIdxs;    
                    //         newInsertedPos.content = vector<char>(snps[indexQuery-1].content.size() , '-');
                    //         newInsertedPos.content[newInsertedPos.content.size()-1] = polishingReads[n][indexTarget];
                    //         newInsertedPos.pos = l;
                    //         snps.push_back(newInsertedPos);
                    //     }
                    //     else{
                    //         snps[insertionPos[10000*indexQuery+numberOfInsertionsThere]].readIdxs.push_back(n);
                    //         snps[insertionPos[10000*indexQuery+numberOfInsertionsThere]].content.push_back(polishingReads[n][indexTarget]);
                    //     }
                    //     numberOfInsertionsThere ++;
                    // }

                    previous_previous_previous_char = previous_previous_char;
                    previous_previous_char = previous_char;
                    previous_char = polishingReads[n][indexTarget];

                    indexTarget++;
                    totalDistance += 1;
                    
                    localDivergence += 1;
                    divergences[l%lengthOfWindow] = true;
                }

                size_of_deletion = 0;

            }

            // if (l%100 == 0 && l >= lengthOfWindow && n == 5){
            //     cout << localDivergence << " " << alignment.substr(l-lengthOfWindow, lengthOfWindow) << " " << l << "\n";
            // }
            //cout << l << " " << result.alignmentLength <<  "\n";
        }

        positionOfReads[n].second = indexQuery;
        
        // cout << " et voiqsl " << n << endl;


        // if ((n == 22 || n == 109 || n == 4 || n == 6 || n == 131)){
        //     cout << endl;
        // }

        readLimits[bbcontig] = positionOfReads;
    }

    string newRef = "";
    unsigned char previous_previous_previous_char = 'A';
    unsigned char previous_previous_char = 'C';
    unsigned char previous_char = 'G';
    for(auto i : consensus){
        previous_previous_previous_char = previous_previous_char;
        previous_previous_char = previous_char;
        previous_char = i;
        newRef += (unsigned char) ('!' + 5*ACGT.find(previous_previous_previous_char) + ACGT.find(previous_previous_char) + 25*ACGT.find(previous_char));
    }
    newref = newRef;

    //print snps (just for debugging)
    // cout << "Printing SNPs, in split_read: cicizzx" << endl;
    // int step = 1; //porportions of reads
    // int prop = 1; //proportion of positions
    // int firstRead = 0;
    // int lastRead = 10;
    // int numberOfReads = lastRead-firstRead;
    // int start = 200;
    // int end = 400;
    // vector<string> reads (int(numberOfReads/step));
    // string cons = "";
    // for (unsigned int i = start ; i < end; i+=prop){
    //     for (short n = 0 ; n < numberOfReads ; n+= step){
    //         unsigned char c = ' ';
    //         int ri = 0;
    //         int soughtRead = firstRead+n;
    //         for (auto r : snps[i].readIdxs){
    //             if (r == soughtRead){
    //                 c = snps[i].content[ri];
    //             }
    //             ri ++;
    //         }
    //         reads[n/step] += std::min(c, (unsigned char) 126);
    //     }
    //     // for (short insert = 0 ; insert < min(9999,numberOfInsertionsHere[i]) ; insert++ ){
    //     //     int snpidx = insertionPos[10000*i+insert];
    //     //     for (short n = 0 ; n < numberOfReads*step ; n+= step){
    //     //         char c = ' ';
    //     //         int ri = 0;
    //     //         for (auto r : snps[snpidx].readIdxs){
    //     //             if (r == n){
    //     //                 c = snps[snpidx].content[ri];
    //     //             }
    //     //             ri ++;
    //     //         }
    //     //         reads[n/step] += c;
    //     //     }
    //     // }
    // }
    // cout << "Here are the aligned reads : " << endl;
    // int index = firstRead;
    // for (auto neighbor : reads){
    //     if (neighbor[neighbor.size()-1] != ' '){
    //         cout << neighbor << " " << index << " " << allreads[allOverlaps[allreads[bbcontig].neighbors_[index]].sequence1].name << endl;
    //     }
    //     index+= step;
    // }
    // int n =start;
    // for(unsigned char i : consensus.substr(start, end-start)){
    //     cout << newRef[n];
    //     n+=prop;
    // } cout << endl;
    // cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;
    // exit(1);
    
    
    return totalDistance/totalLengthOfAlignment;

    //*/
}


int main(int argc, char *argv[])
{

    if (argc != 8){
        std::cout << "Usage: ./filter_variants <columns> <error_rate> <threads> <DEBUG> <col_out> <vcf_in> <vcf_out>\n";
        return 1;
    }

    std::string columnsFile = argv[1];
    float meanError = atof(argv[2]);
    int threads = atoi(argv[3]);
    bool DEBUG = bool(atoi(argv[4]));
    string output_file = argv[5];
    string vcf_in = argv[6];
    string vcf_out = argv[7];

    std::unordered_map<string, int> name_of_contigs;
    std::unordered_map<int, std::string> name_of_contigs2;
    std::vector<std::vector<Column>> snps_in;
    std::vector<int> numberOfReads; //useless here
    parse_column_file(columnsFile, snps_in, name_of_contigs, name_of_contigs2, numberOfReads);

    // cout << "name of contigs : " << endl;
    // for (auto c : name_of_contigs){
    //     cout << c.first << " " << c.second << endl;
    // }
    // cout << snps_in.size() << endl;
    // exit(1);

    std::vector<std::vector<Column>> snps_out;
    keep_only_robust_variants(snps_in, snps_out, name_of_contigs, meanError, threads);

    output_new_column_file(columnsFile, snps_out, name_of_contigs, output_file);

    output_new_vcf_file(vcf_in, snps_out, name_of_contigs2, vcf_out);

}
