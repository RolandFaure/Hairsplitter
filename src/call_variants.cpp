#include "call_variants.h"

#include <iostream>
#include <fstream> //for reading files
#include <sstream> //for parsing strings
#include <string>
#include <omp.h> //for efficient parallelization
#include <set>
#include <cmath>


#include "input_output.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ifstream;
using robin_hood::unordered_map;
using std::pair;
using std::make_pair;


string convert_base_to_triplet(int base){

    base -= 33;
    char middle = "ACGT-"[base%5];
    base /= 5;
    char first = "ACGT-"[base%5];
    base /= 5;
    char last = "ACGT-"[base%5];
    
    return string(1,first)+string(1,middle)+string(1,last);
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
    allreads[bbcontig].new_backbone(make_pair(backboneReadIndex, allreads[bbcontig].neighbors_.size()), allreads[bbcontig].neighbors_.size()+1);
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
    vector <pair <int,int>> positionOfReads; //position of polishing reads on the consensus
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

/**
 * @brief Rudimentary variant calling based on the pileup of the reads
 * 
 * @param snps Pilup of the reads
 * @param suspectPostitions Vector of positions where there are suspected variants
 * @param meanError error rate of the reads
 * @param tmpFolder folder to write the temporary files
 */
vector<Column> call_variants(
    std::vector<Column> &snps, 
    std::vector <Read> &allreads,
    std::vector <Overlap> &allOverlaps, 
    long int contig,
    std::string &ref,
    std::vector<size_t> &suspectPostitions, 
    float &meanError,
    std::string &tmpFolder,
    bool DEBUG){

    vector<int> suspectPositions;
    vector<Column> suspiciousColumns;

    int minimumNumberOfReadsToBeConsideredSuspect = 5;
    if (meanError < 0.015){ //HiFi reads, we can be more stringent
        minimumNumberOfReadsToBeConsideredSuspect = 3;
    }

    
    double depthOfCoverage = 0;
    int posoflastsnp = -5;
    for (int position = 0 ; position < ref.size() ; position++){ 

        if (DEBUG && position%100 == 0){
            cout << "Going through the positions, " << position << "/" << ref.size() << "          \r" << std::flush;
        }
        //count how many time each char appears at this position
        unordered_map<unsigned char, int> content;
        // bool here496 = false;
        for (short n = 0 ; n < snps[position].content.size() ; n++){
            unsigned char base = snps[position].content[n];
            if (content.find(base) == content.end()){
                content[base] = 0;
            }
            if (base != ' '){
                content[base] += 1;
                depthOfCoverage += 1;
            }
            // if (snps[position].readIdxs[n] == 2114){
            //     here496 = true;
            // }
        }
        content[(unsigned char) 0 ] = 0; //to make sure that there are at least 3 different chars
        content[(unsigned char) 1 ] = 0; //to make sure that there are at least 3 different chars
        content[(unsigned char) 2 ] = 0; //to make sure that there are at least 3 different chars

        //find the most frequent chars in content
        vector<pair<unsigned char, int>> content_sorted;
        for (auto it = content.begin() ; it != content.end() ; it++){
            content_sorted.push_back(make_pair(it->first, it->second));
        }
        std::sort(content_sorted.begin(), content_sorted.end(), [](const pair<unsigned char, int>& a, const pair<unsigned char, int>& b) {return a.second > b.second;});

        if (content_sorted.size() > 1){
            snps[position].ref_base = content_sorted[0].first;
            snps[position].second_base = content_sorted[1].first;
            snps[position].pos = position;
        }
        // if (position == 654 ){//&& allreads[contig].name == "edge_2"){
        //     cout << "at pos " << position << " the two mcall_variatns.cppost frequent bases are : " << convert_base_to_triplet((int) content_sorted[0].first) << " "
        //          << convert_base_to_triplet( (int) content_sorted[1].first ) << " " << convert_base_to_triplet( (int) content_sorted[2].first ) << " " <<
        //          convert_base_to_triplet((int) content_sorted[3].first) << endl;
        //     cout << "with frequencies : " << content_sorted[0].second << " " << content_sorted[1].second << " " << content_sorted[2].second << " " << content_sorted[3].second<< endl;
        // }
        if (content_sorted[1].second > minimumNumberOfReadsToBeConsideredSuspect //there is a "frequent" base other than the ref base
            && (content_sorted[1].second > content_sorted[2].second * 5 || minimumNumberOfReadsToBeConsideredSuspect == 2) //this other base occurs much more often than the third most frequent base (could it really be chance ?)
            && content_sorted[0].first%5 != content_sorted[1].first%5 //the central base differs
            && ((content_sorted[1].first - '!')%5 != 4 || (content_sorted[1].first/5%5 != content_sorted[0].first%5 && content_sorted[1].first/25%5 != content_sorted[0].first%5) ) //don't call indels that are adjacent to homopolymers, it's the best way to call false positives
            && position - posoflastsnp > 5){ //the snp is not too close to the previous one
            
            posoflastsnp = position;
            suspectPositions.push_back(position);

            char ref_base;
            if (position < ref.size()){
                ref_base = ref[position];
            } else {
                ref_base = '-';
            }
            // cout << "\nAt first snp, dlsq , " << position << " " << (int) content_sorted[0].first << " " << (int) content_sorted[1].first << " " << (int) ref_base << endl;
            // int code = content_sorted[0].first-'!';
            // cout << "first alternative : " << "ACGT-"[(code%25) / 5] <<  "ACGT-"[code%5] << "ACGT-"[code/25] << endl;
            // cout << "second alternative : " << "ACGT-"[((content_sorted[1].first-'!')%25)/5] << "ACGT-"[(content_sorted[1].first-'!')%5] << "ACGT-"[((content_sorted[1].first-'!')/25)] << endl;
            // exit(1);
            //find the most frequent char in the reads (eexcept the ref base)
            char second_frequent_base = content_sorted[0].first;
            if (second_frequent_base == ref_base){
                second_frequent_base = content_sorted[1].first;
            }

            Column snp;
            snp.pos = position;
            snp.ref_base = ref_base;
            snp.second_base = second_frequent_base;
            auto content_tmp = snps[position].content;
            int n = 0;
            for (int read : snps[position].readIdxs){
                snp.readIdxs.push_back(read);
                snp.content.push_back(content_tmp[n]);
                n++;
            }

            suspiciousColumns.push_back(snp);

        }
    }

    allreads[contig].depth = depthOfCoverage/ref.size();
    return suspiciousColumns;
}

/**
 * @brief Keep only robust variants, i.e. variants that occur on several different positions in the contig
 * 
 * @param snps_in Input snps
 * @param snps_out Filtered snps
 * @param num_threads 
 */
void keep_only_robust_variants(
    std::vector<Column> &msa,
    std::vector<Column>  &snps_in, 
    std::vector<Column>  &snps_out, 
    float mean_error,
    std::vector <Partition> &parts){

    snps_out = vector<Column>();

    // cout << "filtering contqoflmj ig " << n << endl;

    vector<Partition> partitions;
    int lastposition = -5;

    //iterate over all SNPs in the contig
    for (auto snp : snps_in){
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
        return;
    }
    
    float threshold = std::min(4, std::max(2, (int) (mean_error*100)));

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
                    && dis.n10 < std::max(2,2*dis.n01) && dis.n01 < std::max(2,2*dis.n10)){
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
    float chisquare_sum = 0;
    int number_of_interesting_positions = 0;
    for (auto snp : snps_in){
        // cout << "suspecct possisssiion : " << snp.pos << endl;
        // print_snp(snp);
        for (auto p = 0 ; p < listOfFinalPartitions.size() ; p++){
            distancePartition dis = distance(listOfFinalPartitions[p], snp, snp.ref_base);

            // cout << "distance between " << snp.pos << endl;
            // // listOfFinalPartitions[p].print();
            // // print_snp(snp);
            // cout << " is " << dis.n00 << " " << dis.n01 << " " << dis.n10 << " " << dis.n11 << " " << " " << 0.5*snp.content.size() << computeChiSquare(dis)<< endl;

            float chisqu = computeChiSquare(dis);
            if (dis.n00 + dis.n01 + dis.n10 + dis.n11 > 0.5*snp.content.size() 
                && chisqu > 15){

                snps_out.push_back(snp);
                chisquare_sum += chisqu;
                number_of_interesting_positions += 1;
                break;
            }
        }
    }

    double mean_chisquare = chisquare_sum/number_of_interesting_positions; //now see if there aren't a few position we missed
    int idxSnps = 0;
    vector<Column> snps_out_tmp = snps_out;
    snps_out = vector<Column>();

    for (auto position = 0 ; position < msa.size() ; position++){
        if (idxSnps < snps_out_tmp.size() && snps_out_tmp[idxSnps].pos == position){ //this was already suspicious
            snps_out.push_back(snps_out_tmp[idxSnps]);
            idxSnps++;
        }
        else{  //now see if this should be rescued
            if (msa[position].ref_base %5 !=  msa[position].second_base%5 //the central base differs
            && ((msa[position].second_base - '!')%5 != 4 || (msa[position].second_base/5%5 != msa[position].ref_base%5 && msa[position].second_base/25%5 != msa[position].ref_base%5) )) //don't call indels that are adjacent to homopolymers, it's the best way to call false positives
            {
                for (auto p = 0 ; p < listOfFinalPartitions.size() ; p++){
                    distancePartition dis = distance(listOfFinalPartitions[p], msa[position], msa[position].ref_base);
                    // if (position == 502) {
                    //     cout << "dfjlkiox " << dis.n10 << " " << dis.n00 << " " << dis.n01 << " " << dis.n11 << " " << computeChiSquare(dis) << " " << mean_chisquare << endl;
                    //     // print_snp(msa[position]);
                    //     exit(1);
                    // }
                    if (computeChiSquare(dis) > std::min(mean_chisquare, 20.0)){
                        snps_out.push_back(msa[position]);
                        // cout << "fkldjsqdsmlj new sus pos : " << position << " " << computeChiSquare(dis) << endl;
                        break;
                    }
                }
            }
        }
    }

    // cout << "number of intesdresting positions : " << snps_out[n].size() << endl;
    parts = listOfFinalPartitions;
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
 * @brief output the variants in a .col and a .vcf file
 * 
 * @param variants map associating the index of the contig to the list of variants on this contig
 * @param allreads vector containing all the reads (including the contigs)
 * @param allOverlaps vector containing all the overlaps
 * @param col_file output file for the variants in .col format
 * @param vcf_file output file for the variants in vcf format
 */
void output_files(std::unordered_map<int, std::vector<Column>> &variants, std::vector<Overlap> &allOverlaps, std::vector<Read> &allreads, std::string col_file, std::string vcf_file){

    std::ofstream out(col_file);
    std::ofstream vcf(vcf_file);
    
    for (auto c : variants){
        long int contig = c.first;
        vector<Column> suspiciousColumns = c.second;

        //output the suspect positions
        out << "CONTIG\t" << allreads[contig].name << "\t" << allreads[contig].sequence_.size() << "\t" << allreads[contig].depth << "\n";
        for (long int n : allreads[contig].neighbors_){
            out << "READ\t" << allreads[ allOverlaps[n].sequence1 ].name 
                << "\t" << allOverlaps[n].position_1_1 << "\t" << allOverlaps[n].position_1_2 
                << "\t" << allOverlaps[n].position_2_1 << "\t" << allOverlaps[n].position_2_2
                << "\t" << allOverlaps[n].strand << "\n";
        }
        int numberOfReads = allreads[contig].neighbors_.size();
        //output the list of reads
        for (auto c : suspiciousColumns){
            out << "SNPS\t" << c.pos << "\t" << (int) c.ref_base << "\t" << (int) c.second_base << "\t";
            string idxs;
            string bases;
            for (auto r = 0 ; r < c.readIdxs.size() ; r++){
                idxs += std::to_string(c.readIdxs[r]) + ",";
                bases += std::to_string((int) c.content[r]) + ",";
            }
            out << idxs << "\t" << bases << "\n";

            //output the vcf file
            vcf << allreads[contig].name << "\t" << c.pos << "\t.\t" << "ACGT-"[(c.ref_base-'!') % 5]<< "\t" << "ACGT-"[(c.second_base-'!') % 5] << "\t.\t.\tDP=" << c.readIdxs.size()
                << "\n";
        }
        out << endl;
        vcf << endl;
    }
}

int main(int argc, char *argv[])
{
    std::vector <Read> allreads; 
    robin_hood::unordered_map<std::string, unsigned long int> indices;
    vector<unsigned long int> backbone_reads;
    std::vector <Overlap> allOverlaps;
    vector <Link> allLinks;

    if (argc < 10){
        std::cout << "Usage: ./call_variants <gfa_file> <reads_file> <sam_file> <num_threads> <tmpDir> <error_rate_out> <DEBUG> <file_out> <vcfFile>\n";
        return 0;
    }
    std::string gfafile = argv[1];
    std::string readsFile = argv[2];
    std::string samFile = argv[3];
    int num_threads = std::stoi(argv[4]);
    std::string tmpFolder = argv[5];
    std::string error_rate_out = argv[6];
    bool DEBUG = bool(std::stoi(argv[7]));
    std::string file_out = argv[8];
    std::string vcfFile = argv[9];
    //erase the output files
    std::ofstream out(file_out);
    out.close();
    //write the header of the vcf file
    std::ofstream vcf(vcfFile);
    vcf << "##fileformat=VCFv4.2\n";
    vcf << "##source=call_variants\n";
    vcf << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
    vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    vcf.close();

    cout << " - Loading all reads from " << readsFile << " in memory\n";
    parse_reads(readsFile, allreads, indices);

    cout << " - Loading all contigs from " << gfafile << " in memory\n";
    parse_assembly(gfafile, allreads, indices, backbone_reads, allLinks);

    cout << " - Loading alignments of the reads on the contigs from " << samFile << "\n";
    if (samFile.substr(samFile.size()-4,4) == ".paf"){
        cout << "ERROR: please provide a .sam file as input for the alignments of the reads on the contigs." << endl;
        exit(EXIT_FAILURE);
        // parse_PAF(alnOnRefFile, allOverlaps, allreads, indices, backbone_reads, false);
    }
    else if (samFile.substr(samFile.size()-4,4) == ".sam"){
        parse_SAM(samFile, allOverlaps, allreads, indices);
    }
    else{
        cout << "ERROR: the file containing the alignments on the assembly should be .sam" << endl;
        exit(EXIT_FAILURE);
    }

    cout << " - Calling variants on each contig\n";

    int index = 0;
    float totalErrorRate = 0;
    int numberOfContigsWHereErrorRateIsComputed = 0;
    std::unordered_map<int, vector<Column>> variants;

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
    #pragma omp for
        for (int contigIndex = 0 ; contigIndex < backbone_reads.size() ; contigIndex++){

            long int contig = backbone_reads[contigIndex];
            if (allreads[contig].name != "edge_522@00"){ //for debugging purposes

                // if (DEBUG){
                //     #pragma omp critical
                //     {
                //         cout << "Looking at contig number " << index << " out of " << backbone_reads.size() << " (" << allreads[contig].name << ")" << ". By thread " 
                //             << omp_get_thread_num() << ", " << allreads[contig].neighbors_.size() << " reads align here.\n";
                //     }
                // }
                
                #pragma omp critical
                {
                    parse_reads_on_contig(readsFile, contig, allOverlaps, allreads);
                }
    
                vector<Column> snps;  //vector containing list of position, with SNPs at each position
                //build a pileup of the reads on the contig

                robin_hood::unordered_map<int, int> insertionPositions;
                string ref3mers;
                std::unordered_map <int, std::vector<std::pair<int,int>>> readLimits;
                
                float meanDistance = generate_msa(contig, allOverlaps, allreads, snps, insertionPositions, 
                    contigIndex, readLimits, ref3mers, tmpFolder, DEBUG);
                
                #pragma omp critical (errorRate)
                {
                    if (meanDistance > 0){ //to avoid contigs with no reads aligned
                        totalErrorRate += meanDistance;
                        numberOfContigsWHereErrorRateIsComputed += 1;
                    }
                }
                
                //call variants

                vector<size_t> suspectPostitions;
                vector<Column> variants_here = call_variants(snps, allreads, allOverlaps, contig, ref3mers, suspectPostitions, meanDistance, tmpFolder, DEBUG);

                //filter the variants
                vector<Column> filteredSnps;
                vector<Partition> partitions;
                keep_only_robust_variants(snps, variants_here, filteredSnps, meanDistance, partitions);

                // //now do a second round variant calling to rescue the positions that the first round missed
                // if (partitions.size() > 0){
                //     vector<Column> variants_here2 = rescue_snps(snps, meanDistance, partitions, suspectPostitions);
                // }

                variants[contig] = filteredSnps;

                //free up memory by deleting the sequence of the reads used there
                string empty = "";
                for (auto n : allreads[contig].neighbors_){
                    if (allOverlaps[n].sequence1 != contig){
                        allreads[allOverlaps[n].sequence1].set_sequence(empty);
                    }
                    else{
                        allreads[allOverlaps[n].sequence2].set_sequence(empty);
                    }
                }
        
            }
            index++;
        }

    }

    //output the errorRate
    std::ofstream errorRateFile;
    errorRateFile.open(error_rate_out);
    errorRateFile << totalErrorRate/numberOfContigsWHereErrorRateIsComputed << endl;
    errorRateFile.close();

    //output the variants
    output_files(variants, allOverlaps, allreads, file_out, vcfFile);

    return 0;

}

