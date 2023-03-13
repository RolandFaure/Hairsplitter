#include "check_overlaps.h"
// #include "spoa/spoa.hpp"
//#include "cluster_graph.h"
// using namespace wfa;

#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iterator>
#include <omp.h> //for efficient parallelization
#include <map>

#include "input_output.h"
#include "tools.h"

using std::string;
using std::cout;
using std::endl;
using std::list;
using std::vector;
using std::min;
using std::max;
using std::begin;
using std::end;
using std::pair;
using std::make_pair;
using std::unordered_map;
using std::set;
using std::max_element;
using std::thread;
using std::ref;
using std::ofstream;
using std::cerr;

using namespace std::chrono;

extern string MINIMAP; //path to the minimap executable
// extern string MINIASM; //path to the miniasm executable
extern string RACON; //path to the racon executable
extern string HAIRSPLITTER; //path to the hairsplitter folder
extern bool DEBUG; //running debug mode ?

//definition of a small struct that will be useful later
struct distPart{
    float distance = 1;
    short phased = 1;
};
bool comp (distPart i, distPart j){
    return i.distance < j.distance;
}

//input : the set of all overlaps and the backbone reads
//output : a partition for all backbone reads. 
//All reads also have updated backbone_seqs, i.e. the list of backbone reads they are leaning on
//readLimits contains the border of all reads aligning on each backbone, used later to recompute the coverage
void checkOverlaps(std::string fileReads, std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, 
        vector<unsigned long int> &backbones_reads, std::unordered_map<unsigned long int, vector< pair<pair<int,int>, pair<vector<int>, unordered_map<int, string>>  > >> &partitions, 
        bool assemble_on_assembly, unordered_map <int, vector<pair<int,int>>> &readLimits,
        bool polish, int num_threads){

    //main loop : for each backbone read, build MSA (Multiple Sequence Alignment) and separate the reads
    int index = 0;

omp_set_num_threads(num_threads);
#pragma omp parallel
    {
#pragma omp for
        for (long int read : backbones_reads){

            cout << "Looking at backbone read number " << index << " out of " << backbones_reads.size() << " (" << allreads[read].name << ")" << ". By thread " << omp_get_thread_num() << ", " << allreads[read].neighbors_.size() << " reads align here." << endl;
            
            if (allreads[read].neighbors_.size() > 10 && allreads[read].name == "edge_3@0@0" ){

                if (DEBUG){
                    #pragma omp critical
                    {
                        cout << "Looking at backbone read number " << index << " out of " << backbones_reads.size() << " (" << allreads[read].name << ")" << ". By thread " << omp_get_thread_num() << "\n";
                    }
                }
                
                compute_partition_on_this_contig(fileReads, read, allreads, allOverlaps, backbones_reads, partitions, assemble_on_assembly,
                    readLimits, polish);
        
            }
            index++;
        }
    }
    cout << "                                                                       \n";
}

/**
 * @brief Computes the partition of contig and adds it in the map "partitions"
 * 
 * @param fileReads file containing all the reads
 * @param contig index of the contig/read we're aligning on 
 * @param allreads 
 * @param allOverlaps 
 * @param backbones_reads 
 * @param partitions map containing the partitions (with all their interval) of all backbone reads
 * @param assemble_on_assembly set to false if assembling without a reference
 * @param readLimits 
 * @param polish set to true if assembly needs to be polished
 * @param thread identifier of the thread running this function
 */
void compute_partition_on_this_contig(
    std::string fileReads, 
    long int contig, 
    std::vector <Read> &allreads,
    std::vector <Overlap> &allOverlaps, 
    std::vector<unsigned long int> &backbones_reads, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions,
    bool assemble_on_assembly,
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits,
    bool polish
    ){

    parse_reads_on_contig(fileReads, contig, allOverlaps, allreads);
    
    vector<Column> snps;  //vector containing list of position, with SNPs at each position
    //first build an MSA
    string truePar; //for debugging

    vector<bool> misalignedReads(allreads[contig].neighbors_.size(), false);
    string consensus;


    // while building the MSA, check if some reads are too far away from the consensus and build another MSA for them
    // cout << "...creating MSA\n";
    vector<string> consensuses; //vector containing all the consensus sequences (we hope jeust one, but we might have to do several MSA if the reads are too far away from each other)
    robin_hood::unordered_map<int, int> insertionPositions;
    string ref3mers;
    float meanDistance = generate_msa(contig, allOverlaps, allreads, snps, insertionPositions, 
        partitions.size(), truePar, assemble_on_assembly, readLimits, misalignedReads, polish, ref3mers);

    //then separate the MSA
    // string ref = allreads[contig].sequence_.str();
    vector<pair<pair<int,int>, vector<int>> > par = separate_reads(ref3mers, snps,
                                                        allreads[contig].neighbors_.size()+1-int(assemble_on_assembly));
    
    // cout << "hereee are the dkdk partitions : " << endl;
    // for (auto interval : par){
    //     cout << interval.first.first << " <-> " << interval.first.second << " : " << endl;
    //     for (auto cluster : interval.second){
    //         cout << cluster << " ";
    //     }
    //     cout << endl;
    // }

    //now compute the consensus for each new contig compute_consensus_in_partitions does not work with 3mers for now
    // compute_consensus_in_partitions(contig, par, allreads, allOverlaps, snps, insertionPositions, partitions);
    for (auto interval : par){
        unordered_map <int, string> consensus_sequences;
        for (auto cluster : interval.second){
            consensus_sequences[cluster] = allreads[contig].sequence_.subseq(interval.first.first, interval.first.second-interval.first.first+1).str();
        }
        partitions[contig].push_back(make_pair(interval.first, make_pair(interval.second, consensus_sequences)));
    }


    // print all elements of the map partitions[contig][0].second.second
    // auto map = partitions[contig][0].second.second;
    // for (auto it = map.begin(); it != map.end(); ++it) {
    //     cout << it->first << " => " << it->second << '\n';
    // }

    //partitions[contig] = par;

    // cout << "dfopuq True partition:" << truePar << endl << endl;


    //free up memory by deleting snps
    snps.clear();
    vector<Column>().swap(snps);

    //free up memory by deleting the sequence of the reads used there
    for (auto n : allreads[contig].neighbors_){
        if (allOverlaps[n].sequence1 != contig){
            allreads[allOverlaps[n].sequence1].free_sequence();
        }
        else{
            allreads[allOverlaps[n].sequence2].free_sequence();
        }
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
 * @param truePar Vector containing the true partitions (debug only)
 * @param assemble_on_assembly Is backbone a read like any other (false), or rather a reference (true) ?
 * @param readLimits Limits of the reads on the backbone. Used to recompute coverage of the backbone
 * @param misalignedReads Reads aligned with <80% identity: the alignment is not good enough to consider phasing there
 * @param polish True if the backbone reads are not polished
 * @return The mean distance between the aligned reads and the consensus backbone
 */
float generate_msa(long int bbcontig, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, 
    std::vector<Column> &snps, robin_hood::unordered_map<int, int> &insertionPos, int backboneReadIndex, string &truePar, bool assemble_on_assembly, 
    unordered_map <int, vector<pair<int,int>>> &readLimits, std::vector<bool> &misalignedReads, bool polish,
    string &newref){

    string ACGT = "ACGT-";

    //go through the neighbors of the backbone read and align it

    //keep count of the distance between two reads to know the mean distance
    float totalDistance = 0;
    double totalLengthOfAlignment = 0;

    //to remember on what backbone the backbone read is leaning -> itself
    allreads[bbcontig].new_backbone(make_pair(backboneReadIndex, allreads[bbcontig].neighbors_.size()), allreads[bbcontig].neighbors_.size()+1);
    string read_str = allreads[bbcontig].sequence_.str();

    //small loop to compute truePartition DEBUG
    if (DEBUG){
        for (auto n = 0 ; n<allreads[bbcontig].neighbors_.size() ; n++){

            long int neighbor = allreads[bbcontig].neighbors_[n];
            Overlap overlap = allOverlaps[neighbor];
            
            if (overlap.sequence1 == bbcontig){

                truePar.push_back(allreads[overlap.sequence2].name[1]);
                // cout << "name : " << allreads[overlap.sequence2].name << " " << allreads[overlap.sequence2].name[1] << " " << truePartition[truePartition.size()-1] << endl;
            }
            else{
                truePar.push_back(allreads[overlap.sequence1].name[1]);
                // cout << "name : " << allreads[overlap.sequence1].name << " " << allreads[overlap.sequence1].name[1] << " " << truePartition[truePartition.size()-1] << endl;
            }    
        }
    }

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

        if (overlap.CIGAR != "" && !polish){
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
    readLimits[bbcontig] = positionOfReads;

    string consensus;
    if (polish){
        string thread_id = std::to_string(omp_get_thread_num());
        consensus = consensus_reads(read_str , polishingReads, thread_id);
        if (DEBUG){
            cout << "Done polishing the contig" << endl;
        }
    }
    else{
        consensus = read_str; //if the input assembly is already polished
    } 

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

        auto t1 = high_resolution_clock::now();

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

        //a loop going through the CIGAR and modifyning snps

        auto t2 = high_resolution_clock::now();

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

                    if (polishingReads[n][indexTarget] != previous_char){ //to avoid all homopolymers errors
                        previous_previous_previous_char = previous_previous_char;
                        previous_previous_char = previous_char;
                        previous_char = polishingReads[n][indexTarget];
                    }

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

                    if ( '-' != previous_char || size_of_deletion > 10){ //to avoid all homopolymers, except for long deletions
                        previous_previous_previous_char = previous_previous_char;
                        previous_previous_char = previous_char;
                        previous_char = '-';
                    }

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

                    if (polishingReads[n][indexTarget] != previous_char && indexQuery > 0){ //to avoid all homopolymers
                        previous_previous_previous_char = previous_previous_char;
                        previous_previous_char = previous_char;
                        previous_char = polishingReads[n][indexTarget];
                    }

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
        
        // cout << " et voiqsl " << n << endl;


        // if ((n == 22 || n == 109 || n == 4 || n == 6 || n == 131)){
        //     cout << endl;
        // }

        auto t3 = high_resolution_clock::now();

        alignmentTime += duration_cast<milliseconds>(t2-t1).count();
        MSAtime += duration_cast<milliseconds>(t3-t2).count();
        
    }

    string newRef = "";
    unsigned char previous_previous_previous_char = 'A';
    unsigned char previous_previous_char = 'C';
    unsigned char previous_char = 'G';
    for(auto i : consensus){
        if (i != previous_char){
            previous_previous_previous_char = previous_previous_char;
            previous_previous_char = previous_char;
            previous_char = i;
        }
        newRef += (unsigned char) ('!' + 5*ACGT.find(previous_previous_previous_char) + ACGT.find(previous_previous_char) + 25*ACGT.find(previous_char));
    }
    newref = newRef;

    // if (DEBUG){
    //     cout << "Building MSA took time... " << alignmentTime << " for WFA and " << MSAtime << " for filling the vector" << endl;
    //     //print the vector of bool misaligned reads
    //     cout << "jkjksz Misaligned reads : ";
    //     int n = 0;
    //     for (auto r : misalignedReads){
    //         cout << r << " ";
    //         if (!r){
    //             cout << "Read " << allreads[allOverlaps[allreads[read].neighbors_[n]].sequence1].name << " is aligned: " <<" \n";
    //         }
    //         n++;
    //     }
    //     cout << endl;
    // }
    
    /*
    //print snps (just for debugging)
    int step = 1; //porportions of reads
    int prop = 1; //proportion of positions
    int firstRead = 0;
    int lastRead = polishingReads.size();
    int numberOfReads = lastRead-firstRead;
    int start = 47728;
    int end = 47828;
    vector<string> reads (int(numberOfReads/step));
    string cons = "";
    for (unsigned int i = start ; i < end; i+=prop){
        for (short n = 0 ; n < numberOfReads ; n+= step){
            unsigned char c = ' ';
            int ri = 0;
            int soughtRead = firstRead+n;
            for (auto r : snps[i].readIdxs){
                if (r == soughtRead){
                    c = snps[i].content[ri];
                }
                ri ++;
            }
            reads[n/step] += min(c, (unsigned char) 126);
        }
        // for (short insert = 0 ; insert < min(9999,numberOfInsertionsHere[i]) ; insert++ ){
        //     int snpidx = insertionPos[10000*i+insert];
        //     for (short n = 0 ; n < numberOfReads*step ; n+= step){
        //         char c = ' ';
        //         int ri = 0;
        //         for (auto r : snps[snpidx].readIdxs){
        //             if (r == n){
        //                 c = snps[snpidx].content[ri];
        //             }
        //             ri ++;
        //         }
        //         reads[n/step] += c;
        //     }
        // }
    }
    cout << "Here are the aligned reads : " << endl;
    int index = firstRead;
    for (auto neighbor : reads){
        if (neighbor[neighbor.size()-1] != ' '){
            cout << neighbor << " " << index << " " << allreads[allOverlaps[allreads[bbcontig].neighbors_[index]].sequence1].name.substr(0,10) << endl;
        }
        index+= step;
    }
    int n =start;
    for(unsigned char i : consensus.substr(start, end-start)){
        cout << newRef[n];
        n+=prop;
    } cout << endl;
    cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;
    exit(1);
    */

    return totalDistance/totalLengthOfAlignment;

    //*/
}

/**
 * @brief Function that uses racon to polish a sequence
 * 
 * @param backbone sequence to be polished
 * @param polishingReads list of reads to polish it
 * @param overhang length of the unpolished ends to be used
 * @param id an id (typically a thread id) to be sure intermediate files do not get mixed up with other threads 
 * @return polished sequence 
 */
string consensus_reads(string const &backbone, vector <string> &polishingReads, string &id){
    
    if (polishingReads.size() == 0){
        return backbone;
    }

    system("mkdir tmp/ 2> trash.txt");
    std::ofstream outseq("tmp/unpolished_"+id+".fasta");
    outseq << ">seq\n" << backbone;
    outseq.close();

    std::ofstream polishseqs("tmp/reads_"+id+".fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        polishseqs << ">read"+std::to_string(read)+"\n" << polishingReads[read] << "\n";
    }
    polishseqs.close();

    // assemble de novo using wtdbg2
    /*
    string comAsm = "/home/rfaure/Documents/software/wtdbg2/wtdbg2 -e 5 -l 1000 -L 3000 -S 1 -R -o tmp/wtdbg2"+id+" -i tmp/reads_"+id+".fasta 2>tmp/trash.txt";
    system(comAsm.c_str());

    string cons_wtdbg2 = "/home/rfaure/Documents/software/wtdbg2/wtpoa-cns -i tmp/wtdbg2"+id+".ctg.lay.gz -fo tmp/dbg.raw"+id+".fa 2>tmp/trash.txt";
    int res_wtdbg2 = system(cons_wtdbg2.c_str());

    if (res_wtdbg2){
        // polish consensus, not necessary if you want to polish the assemblies using other tools
        string comMap = "minimap2 -ax map-pb -r2k tmp/dbg.raw"+id+".fa tmp/reads_"+id+".fasta 2>tmp/trash.txt | samtools sort >tmp/dbg"+id+".bam 2>tmp/trash.txt";
        system(comMap.c_str());

        string comSamtools = "samtools view -F0x900 tmp/dbg"+id+".bam 2>tmp/trash.txt | /home/rfaure/Documents/software/wtdbg2/wtpoa-cns -d tmp/dbg.raw"+id+".fa -i - -fo tmp/dbg.cns"+id+".fa 2>tmp/trash.txt";
        system(comSamtools.c_str());

        string comUnfold = "awk '{if(\">\" == substr($1,1,1)){ printf \"\\n\"; print;} else printf $1;}' tmp/dbg.cns"+id+".fa > tmp/consensus"+id+".fa  2>tmp/trash.txt";
        system(comUnfold.c_str());

        //now read the consensus and return it
        std::ifstream in("tmp/consensus"+id+".fa");
        string consensus;
        string line;
        while (std::getline(in, line)){
            if (line[0] != '>'){
                consensus += line;
            }
        }
        in.close();
        return consensus;
    }

    else{
    */

    string com = " -t 1 -x map-pb tmp/unpolished_"+id+".fasta tmp/reads_"+id+".fasta > tmp/mapped_"+id+".paf 2>tmp/trash.txt";
    string commandMap = MINIMAP + com; 
    system(commandMap.c_str());

    com = " -w 500 -e 1 -t 1 tmp/reads_"+id+".fasta tmp/mapped_"+id+".paf tmp/unpolished_"+id+".fasta > tmp/polished_"+id+".fasta 2>tmp/trash.txt";
    string commandPolish = RACON + com;
    system(commandPolish.c_str());

    //polish using CONSENT
    // string CONSENT = "/home/rfaure/Documents/software/CONSENT/CONSENT-polish";
    // string com = " -S 30 --contigs tmp/unpolished_"+id+".fasta --reads tmp/reads_"+id+".fasta --out tmp/polished_"+id+".paf >tmp/log_CONSENT.txt 2>tmp/trash.txt";
    // string commandPolish = CONSENT + com;
    // system(commandPolish.c_str());

    std::ifstream polishedRead("tmp/polished_"+id+".fasta");
    string line;
    string consensus;
    while(getline(polishedRead, line)){
        if (line[0] != '>'){
            consensus = line; 
        }
    }

    if (consensus == ""){ //that's typically if there are too few reads to consensus
        return backbone;
    }

    //racon tends to drop the ends of the sequence, so attach them back.
    //This is an adaptation in C++ of a Minipolish (Ryan Wick) code snippet 
    auto before_size = min(size_t(300), backbone.size());
    auto after_size = min(size_t(200), consensus.size());

    // Do the alignment for the beginning of the sequence.
    string before_start = backbone.substr(0,before_size);
    string after_start = consensus.substr(0,after_size);

    EdlibAlignResult result = edlibAlign(after_start.c_str(), after_start.size(), before_start.c_str(), before_start.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));


    int start_pos = result.startLocations[0];
    string additional_start_seq = before_start.substr(0, start_pos);

    edlibFreeAlignResult(result);


    // And do the alignment for the end of the sequence.
    string before_end = backbone.substr(backbone.size()-before_size, before_size);
    string after_end = consensus.substr(consensus.size()-after_size , after_size);

    result = edlibAlign(after_end.c_str(), after_end.size(), before_end.c_str(), before_end.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    int end_pos = result.endLocations[0]+1;
    string additional_end_seq = before_end.substr(end_pos , before_end.size()-end_pos);
    edlibFreeAlignResult(result);
    
    // cout << "consensus : " << endl << (additional_start_seq + consensus + additional_end_seq).substr(100,150) << endl;
    // cout << backbone.substr(100,150) << endl;
    
    return additional_start_seq + consensus + additional_end_seq;

}

/**
 * @brief separates the reads of the MSA into groups
 * 
 * @param ref reference sequence against which the reads are aligned (in 3mer space)
 * @param snps MSA
 * @param meanDistance estimation of the error rate
 * @param numberOfReads number of reads of the MSA
 * @return vector<pair<pair<int,int>, vector<int>> >  pairs of coordinates to which are attached a certain partition of reads
 */
vector<pair<pair<int,int>, vector<int>> > separate_reads(string& ref, std::vector<Column> &snps, int numberOfReads){

    vector<Partition> partitions; //list of all partitions of the reads, with the number of times each occurs

    vector <bool> no_mask (numberOfReads, true);
    vector<size_t> suspectPostitions;
    float meanError = 0;
    partitions = get_solid_partitions(ref, snps, no_mask, suspectPostitions, meanError, numberOfReads);

    // cout << "partition dnooww: " << partitions.size() << endl;
    // int n = 0;
    // for (auto p : partitions){
    //     cout << "dqdsoi " << n << " ";
    //     p.print();
    //     n++;
    // }

    // exit(1);

    // vector <Partition> strengthened_partition = second_pass(ref, partitions, snps, suspectPostitions, meanError, numberOfReads);

    // cout << "partition strengththththened: " << strengthened_partition.size() << endl;
    // for (auto p : strengthened_partition){
    //     p.print();
    // }
    // exit(0);

    // there is a first list of patitions, now check if there are some sub-partitions that were missed

    // int numberOfNewPartitions = partitions.size();
    // vector<vector<bool>> all_masks_ever_tested; //to avoid testing the same mask twice
    // while (numberOfNewPartitions > 0){
    //     numberOfNewPartitions = 0;
        
    //     // partitions = {partitions[0], partitions[1]};
    //     vector<vector<bool>> masks = create_masks(partitions, numberOfReads, all_masks_ever_tested);
    //     // cout << "qaeefzqf here are all the masks : " << endl;
    //     // for (auto mask : masks){
    //     //     for (int i = 0 ; i < numberOfReads ; i+=1){
    //     //         cout << mask[i];
    //     //     }
    //     //     cout << endl;
    //     // }
    //     // cout << endl << endl << endl;
    //     all_masks_ever_tested.insert(all_masks_ever_tested.end(), masks.begin(), masks.end());
    //     for (auto mask : masks){
    //         float meanErrorHere = 0;
    //         vector<Partition> newPartitions = get_solid_partitions(ref, snps, mask, suspectPostitions, meanErrorHere, numberOfReads);
    //         // cout << "mask : reads, hereeee they are:" << endl;
    //         // for (int i = 0 ; i < numberOfReads ; i+=1){
    //         //     cout << mask[i];
    //         // }
    //         // cout << endl;
    //         // cout << "new partidssdstions : " << newPartitions.size() << endl;
    //         // for (auto p : newPartitions){
    //         //     p.print();
    //         // }
    //         // exit(1);
    //         numberOfNewPartitions += newPartitions.size();
    //         partitions.insert(partitions.end(), newPartitions.begin(), newPartitions.end());
    //         // break;
    //     }    
    // }

    // //now we have the list of final partitions : there may be several, especially if there are more than two haplotypes

    // if (DEBUG){
    //     cout << "final dddaaz partitions : " << partitions.size() << endl;
    //     for (auto p : partitions){
    //         p.print();
    //     }
    // }

    //now go through windows of width 1000 along the reference and create local partitions

    //list the interesting positions
    vector<size_t> interestingPositions;
    for (auto pos : suspectPostitions){
        // cout << "suspecct possisssiion : " << endl;
        // print_snp(snps[pos], no_mask);
        for (auto p = 0 ; p < partitions.size() ; p++){
            distancePartition dis = distance(partitions[p], snps[pos], ref[pos]);
            if (dis.n00 + dis.n01 + dis.n10 + dis.n11 > 0.5*snps[pos].content.size() 
                && computeChiSquare(dis) > 5 
                && ((interestingPositions.size()> 0 && interestingPositions[interestingPositions.size()-1] != pos) || interestingPositions.size() == 0)){
                // cout << "posssz : ";
                // print_snp(snps[pos], mask_at_this_position);
                // if (pos ==  4224){
                //     cout << "correlationg with " << ref_base << " " << dis.n00 << " " << dis.n01 << " " << dis.n10 << " " << dis.n11 << endl;
                //     partitions[p].print();
                //     exit(1);
                // }
                
                interestingPositions.push_back(pos);
                break;
            }
        }
    }
    cout << "listing similaritiees and differences between reads ddd ..." << endl;
    vector<vector<pair<int,int>>> sims_and_diffs (numberOfReads, vector<pair<int,int>>(numberOfReads, make_pair(0,0)));
    list_similarities_and_differences_between_reads(no_mask, snps, interestingPositions, sims_and_diffs);

    // cout << "interesting posiaaations : " <<  endl;
    // for (auto pos : interestingPositions){
    //     cout << pos << " ";
    // }
    // cout << endl;

    // cout << "similrities and differences computed between reads 0 and other" << endl;
    // for (auto i : sims_and_diffs[0]){
    //     cout << i.first << " " << i.second << endl;
    // }

    vector<pair<pair<int,int>, vector<int>>> threadedReads;
    int suspectPostitionIdx = 0;
    int sizeOfWindow = 2000;
    for (auto chunk = 0 ; chunk < double(ref.size())/sizeOfWindow ; chunk++){
        vector<Partition> localPartitions;

        //only consider reads that span the whole length of the window: create a mask
        vector<bool> mask_at_this_position (numberOfReads, false); 
        if ((chunk+1)*sizeOfWindow < ref.size()){
            int n = 0;
            for (auto r : snps[(chunk+1)*sizeOfWindow].readIdxs){
                if (snps[(chunk+1)*sizeOfWindow].content[n] != ' '){
                    mask_at_this_position[r] = true;
                }
                n++;
            }
        }
        else{
            mask_at_this_position = no_mask;
        }
        int index = 0;
        if (chunk > 0){ //at position exactly 0 there are very few reads, so don't exclude any reads based on that
            for (auto r = 0 ; r < mask_at_this_position.size(); r++){
                while(index < snps[chunk*sizeOfWindow].readIdxs.size() && snps[chunk*sizeOfWindow].readIdxs[index] < r){
                    index++;
                }
                if (index >= snps[chunk*sizeOfWindow].readIdxs.size() || snps[chunk*sizeOfWindow].readIdxs[index] != r){
                    mask_at_this_position[r] = false;
                }
            }
        }

        // for (auto p = chunk*sizeOfWindow ; p < (chunk+1)*sizeOfWindow ; p++){
        //     int indexsnp = 0;
        //     cout << p << " : ";
        //     for (auto read = 0 ; read<numberOfReads ; read++){
        //         if (mask_at_this_position[read]){
        //             while (indexsnp < snps[p].readIdxs.size() && snps[p].readIdxs[indexsnp] < read){
        //                 indexsnp++;
        //             }
        //             if (indexsnp >= snps[p].readIdxs.size() || snps[p].readIdxs[indexsnp] != read){
        //                 cout << " ";
        //             }
        //             else {
        //                 cout << snps[p].content[indexsnp];
        //             }
        //         }
        //     }
        //     cout << endl;
        // }

        for (auto p : partitions){
            Partition new_partition;
            auto m = p.get_mask();
            new_partition.apply_mask(m);
            if (chunk < ref.size()/sizeOfWindow - 1 && chunk > 0){
                new_partition.apply_mask(mask_at_this_position);
            }
            localPartitions.push_back(new_partition);
        }

        // //let's see what partitions are found on the local window
        // while(suspectPostitionIdx < suspectPostitions.size() && suspectPostitions[suspectPostitionIdx] < (chunk+1)*sizeOfWindow){

        //     auto pos = suspectPostitions[suspectPostitionIdx];
        //     // if (pos < 4000){
        //     //     cout << "her edfkju ifnfif 4000 : " << pos << endl;
        //     //     print_snp(snps[pos], mask_at_this_position);
        //     // }
        //     //find the closest partition
        //     for (auto p = 0 ; p < localPartitions.size() ; p++){
        //         //compute the distance between the partition and the suspect position
        //         if (pos <= partitions[p].get_right() && pos >= partitions[p].get_left()){
        //             char ref_base = ref[pos];
        //             // listOfFinalPartitions[p].print();
        //             //print snps[pos];
        //             // for (auto i : snps[pos].content){
        //             //     cout << i;
        //             // }
        //             // cout << endl;
        //             distancePartition dis = distance(partitions[p], snps[pos], ref_base);

        //             if (( min(dis.n00,dis.n11) >= 5 && float(dis.n01+dis.n10)/min(dis.n00 , dis.n11) < 0.2 && dis.phased == 1) 
        //                 || (min(dis.n01,dis.n10) >= 5 && float(dis.n00+dis.n11)/min(dis.n01 , dis.n10) < 0.2 && dis.phased == -1)){
        //             // if (dis.solid11 >= 5 && dis.solid00 >= 5
        //             //     && dis.solid10 < min(5.0, 0.2*dis.solid11) 
        //             //     && dis.solid01 < min(5.0, 0.2*dis.solid00)
        //             //     && dis.n01 < 0.4*(dis.n00+dis.n01)){
        //                 localPartitions[p].augmentPartition(snps[pos], pos);
        //                 // cout << "augmenting ddddpartition " << p << " with position " << pos << endl;
        //                 // if (chunk == 4 && p == 4){
        //                 //     cout << "on windhhow " << chunk << " found partition " << p << endl;
        //                 //     // for (auto r = 0 ; r < localPartitions[p].getReads().size() ; r++){
        //                 //     //     if (localPartitions[p].getPartition()[r] == 1 || localPartitions[p].getPartition()[r] == -1){
        //                 //     //         cout << localPartitions[p].getReads()[r] << "_" << localPartitions[p].getPartition()[r] << "_" << localPartitions[p].getMore()[r] << "_" << localPartitions[p].getLess()[r] << " ; ";
        //                 //     //     }
        //                 //     // }
        //                 //     // cout << endl;
        //                 //     //print snps[pos].readIdx;
        //                 //     // for (auto i : snps[pos].readIdxs){
        //                 //     //     cout << i << " ";
        //                 //     // }
        //                 //     // cout << endl;
        //                 //     auto index = 0;
        //                 //     for (auto i = 0; i < snps[pos].content.size() ; i++){
        //                 //         while(index < localPartitions[p].getReads().size() && localPartitions[p].getReads()[index] < snps[pos].readIdxs[i]){
        //                 //             index++;
        //                 //         }
        //                 //         if (localPartitions[p].getReads()[index] == snps[pos].readIdxs[i] && localPartitions[p].getPartition()[index] != -2){
        //                 //             cout << snps[pos].content[i];
        //                 //         }
        //                 //     }
        //                 //     cout << endl;
        //                 //     cout << "dis.n00 : " << dis.n10 << " " << dis.n01 << " " << dis.n11 << " " << dis.n00 << endl;
        //                 //     localPartitions[p].print();
        //                 //     // partitions[p].print();
        //                 //     // cout << localPartitions[p].getReads()[localPartitions[p].getReads().size()-1] << endl;
        //                 //     // while(true){}
        //                 //     // if (localPartitions[p].number() > 3){
        //                 //     //     exit(0);
        //                 //     // }
        //                 // }
        //             }
        //         }
        //     }
        //     suspectPostitionIdx += 1;
        // }

        // //go through the local partitions and strengthen them
        // for (auto p = 0 ; p < localPartitions.size() ; p++){
        //     localPartitions[p].strengthen_partition(partitions[p]);
            
        // }

        // list and strengthen local partitions that have a numberOfOccurences >= 1
        // if (chunk == 89){
        //     for (auto p = 0 ; p < localPartitions.size() ; p++){
        //         if (localPartitions[p].number() >= 1){
        //             cout << "on windhhow " << chunk << " found partition " << p << endl;
        //             // localPartitions[p].print();
        //             cout << localPartitions[p].getReads()[localPartitions[p].getReads().size()-1] << endl;
        //         }
        //     }
        //     cout << endl << endl;
        //     while(true){}
        // }

        //create non-null-partitions, compased of local partitions that have a number() >= 1
        // vector<Partition> non_null_partitions;
        // for (auto p = 0 ; p < localPartitions.size() ; p++){
        //     if (localPartitions[p].number() >= 1){
        //         non_null_partitions.push_back(localPartitions[p]);
        //     }
        // }

        //select compatible partitions
        // vector<Partition> listOfCompatiblePartitions = select_compatible_partitions(non_null_partitions, numberOfReads, meanError);
        // cout << "mcopqodn sefjdoi" << endl;
        // for (auto p = 0 ; p < listOfCompatiblePartitions.size() ; p++){
        //     cout << "on window " << chunk << " found coimpatdsible partition " << endl;
        //     listOfCompatiblePartitions[p].print();
        // }

        // for (auto p : listOfCompatiblePartitions){
        //     auto content = p.getPartition();
        //     for (auto i : content){
        //         if (i != 0 && i != 1 && i != -1 && i != -2){
        //             cout << "error in pdljqssqddartition " << endl;
        //             for (auto j : content){
        //                 cout << j << ",";
        //             }
        //             cout << endl;
        //             exit(1);
        //         }
        //     }
        // }



        // thread partition on this window
        // vector<int> clustReads = threadHaplotypes_in_interval(non_null_partitions, numberOfReads);
        // for (auto r = 0 ; r < numberOfReads ; r++){
        //     if (mask_at_this_position[r] == false){
        //         clustReads[r] = -2;
        //     }
        // }

        // cout << "threaded reads  edfqc"  << endl;
        // for (auto r = 0 ; r < clustReads.size() ; r++){
        //     if (clustReads[r] != -2){
        //         cout << clustReads[r] << " ";
        //     }
        // }
        // cout << endl;



        // cout << "fqmlsjdm local partitions" << endl;
        // for (auto p = 0 ; p < non_null_partitions.size() ; p++){
        //     cout << "on windssow " << chunk << " found local partition " << p << endl;
        //     non_null_partitions[p].print();
        // }

        vector<vector<int>> adjacency_matrix (mask_at_this_position.size(), vector<int>(mask_at_this_position.size(), 0));
        create_read_graph(mask_at_this_position, snps, chunk, interestingPositions, sizeOfWindow, sims_and_diffs, adjacency_matrix);
        vector<int> clustersStart (adjacency_matrix.size(), 0);
        for (auto r = 0 ; r < adjacency_matrix.size() ; r++){
            clustersStart[r] = r;
        }
        auto strengthened_adjacency_matrix = strengthen_adjacency_matrix(adjacency_matrix);

        vector<vector<int>> allclusters_debug;
        vector<int> clusteredReads1 = chinese_whispers(adjacency_matrix, clustersStart);
        allclusters_debug.push_back(clusteredReads1);
        vector<vector<int>> localClusters = {};

        // cout << "here are all the interesting positions" << endl;
        // for (auto p : interestingPositions){
        //     cout << p << " ";
        // }
        // cout << endl;

        for (auto position : interestingPositions){
            if (position >= chunk*sizeOfWindow && position < chunk*sizeOfWindow + sizeOfWindow){
                // cout << "in dldjk position " << position << " : " << endl;
                // print_snp(snps[position], mask_at_this_position);

                unordered_map<unsigned char, int> charToIndex;
                vector<int> clustersStart2 (adjacency_matrix.size(), 0);
                for (auto r = 0 ; r < adjacency_matrix.size() ; r++){
                    clustersStart2[r] = r;
                }
                for (auto r = 0 ; r < snps[position].content.size() ; r++){
                    if (mask_at_this_position[snps[position].readIdxs[r]] == true){
                        if (charToIndex.find(snps[position].content[r]) == charToIndex.end()){
                            charToIndex[snps[position].content[r]] = snps[position].readIdxs[r];
                        }
                        int indexHere = charToIndex[snps[position].content[r]];
                        clustersStart2[snps[position].readIdxs[r]] = indexHere;
                    }
                }
                
                vector<int> clusteredReads_local = chinese_whispers(strengthened_adjacency_matrix, clustersStart2); 
                localClusters.push_back(clusteredReads_local);

                if (chunk == 3){

                    // outputGraph(strengthened_adjacency_matrix, clusteredReads_local, "graphs/local_graph_"+std::to_string(position)+".gdf");

                    // cout << "clustered reads local : " << position << endl;
                    // for (auto r = 0 ; r < clusteredReads_local.size() ; r++){
                    //     if (mask_at_this_position[r] == true){
                    //         cout << clusteredReads_local[r] << " ";
                    //     }
                    // }
                    // cout << endl;
                    allclusters_debug.push_back(clusteredReads_local);
                }
            }         
        }

        vector<int> new_clusters = merge_clusterings(localClusters, adjacency_matrix);
        vector<int> clusteredReads = new_clusters;

        // find all clusters that have a size < 5 and give them a -1 cluster, so they can be rescued later
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
            if (clusterSizes[clusteredReads[r]] < 5 && clusteredReads[r] != -2){
                clusteredReads[r] = -1;
            }
        }

        // cout << "clustereReads : " << endl;
        // for (auto r : clusteredReads){
        //     if (r > -1){
        //         cout << r <<" ";
        //     }
        // }
        // cout << endl;

        vector<int> mergedHaplotypes = merge_wrongly_split_haplotypes(clusteredReads, snps, chunk, interestingPositions, adjacency_matrix, sizeOfWindow);

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

        vector<int> haplotypes = rescue_reads(mergedHaplotypes, snps, chunk, interestingPositions, sizeOfWindow);

        cout << "haploutypes : " << chunk << endl;
        for (auto h : haplotypes){
            if (h > -1){
                cout << h;
            }
            // else{
            //     cout << "_";
            // }
        }
        cout << endl;

        // if (chunk == 14){
        //     outputGraph(adjacency_matrix, clusteredReads, "tmp/adjacency_matrix_14.gdf");
        // }

        if (chunk == 3){
            allclusters_debug.push_back(haplotypes);
            outputGraph_several_clusterings(adjacency_matrix, allclusters_debug, "graphs/cluster_final.gdf");
            exit(1);
        }

        // cout << "already separated qldfjp : " << endl;


        threadedReads.push_back(make_pair(make_pair(chunk*sizeOfWindow, min((chunk+1)*sizeOfWindow-1, int(ref.size()-1))), haplotypes));

    }
    return threadedReads;
    
}

/**
 * @brief Among a certain subset of reads, find the best partition if there is one
 * 
 * @param ref sequence of the original contig
 * @param snps MSA
 * @param mask vector of bool, indicating by "false" what reads not to take
 * @param suspectPostitions list of the suspect positions already found
 * @param meanError will be changed in this function
 * @param numberOfReads 
 * @return Partition 
 */
vector<Partition> get_solid_partitions(std::string& ref, 
    std::vector<Column> &snps,
    vector<bool> &mask,
    vector <size_t> &suspectPostitions,
    float &meanError,
    int numberOfReads){

    vector<Partition> partitions; //list of all partitions of the reads, with the number of times each occurs

    int numberOfSuspectPostion = 0;
    int numberOfNeighbors = 0;

    //two variables to know how much a position diverges from the consensus on average
    int numberOfExtensions = 0;
    double meanErrorHere = 0;
    vector<size_t> suspectPostitionsHere;
    vector<Column> suspiciousColumns;

    vector<unsigned int> debug_reads_of_interest;

    for (int position = 0 ; position < ref.size() ; position++){ 

        if (DEBUG && position%100 == 0){
            cout << "Going through the positions, " << position << "/" << ref.size() << "          \r" << std::flush;
        }
        //count how many time each char appears at this position
        unordered_map<char, int> content;
        int numberOfReadsHere = 0;
        // bool here496 = false;
        for (short n = 0 ; n < snps[position].content.size() ; n++){
            if (mask[snps[position].readIdxs[n]]){
                char base = snps[position].content[n];
                if (content.find(base) == content.end()){
                    content[base] = 0;
                }
                if (base != ' '){
                    content[base] += 1;
                    numberOfReadsHere += 1;
                }
            }
            // if (snps[position].readIdxs[n] == 2114){
            //     here496 = true;
            // }
        }
        content[(unsigned char) 0 ] = 0; //to make sure that there are at least 3 different chars
        content[(unsigned char) 1 ] = 0; //to make sure that there are at least 3 different chars
        content[(unsigned char) 2 ] = 0; //to make sure that there are at least 3 different chars

        //find the most frequent chars in content
        vector<pair<char, int>> content_sorted;
        for (auto it = content.begin() ; it != content.end() ; it++){
            content_sorted.push_back(make_pair(it->first, it->second));
        }
        sort(content_sorted.begin(), content_sorted.end(), [](const pair<char, int>& a, const pair<char, int>& b) {return a.second > b.second;});

        if (content_sorted[1].second > 5 //there is a "frequent" base other than the ref base
            && content_sorted[1].second > content_sorted[2].second * 5 //this other base occurs much more often than the third most frequent base (could it really be chance ?)
            && content_sorted[0].first%5 != content_sorted[1].first%5 ){ //the central base differs
           
            suspectPostitionsHere.push_back(position);

            // cout << "iouocxccv " << threshold << " " << position << " ;bases : " << content[0] << " " << content[1] << " " << content[2] << " " << content[3] << " " << content[4] << endl;
            // cout << "fmqjldqmjl " << endl;
            // print_snp(snps[position], mask);

            char ref_base;
            if (position < ref.size()){
                ref_base = ref[position];
            } else {
                ref_base = '-';
            }
            //find the most frequent char in the reads (eexcept the ref base)
            char second_frequent_base = content_sorted[0].first;
            if (second_frequent_base == ref_base){
                second_frequent_base = content_sorted[1].first;
            }

            //creat the masked snp, with only reads that are not masked
            Column snp;
            snp.pos = position;
            snp.ref_base = ref_base;
            snp.second_base = second_frequent_base;
            auto content_tmp = snps[position].content;
            int n = 0;
            for (int read : snps[position].readIdxs){
                if (mask[read]){
                    snp.readIdxs.push_back(read);
                    snp.content.push_back(content_tmp[n]);
                }
                n++;
            }

            suspiciousColumns.push_back(snp);

            //go through the partitions to see if this suspicious position looks like smt we've seen before
            bool found = false;
            for (auto p = 0 ; p < partitions.size() ; p++){
                //if the partition is too far away, do not bother comparing
                if (std::abs(snp.pos-partitions[p].get_right())>50000){
                    continue;
                }
                distancePartition dis = distance(partitions[p], snp, ref_base);
                auto comparable = dis.n00 + dis.n11 + + dis.n01 + dis.n10;


                //if ((float(dis.n01+dis.n10)/(min(dis.n00,dis.n11)+dis.n01+dis.n10) <= meanDistance*2 || dis.n01+dis.n10 <= 2)  && dis.augmented && comparable > min(10.0, 0.3*numberOfReads)){
                if (dis.n01 < 0.1 * (dis.n00+dis.n01) && dis.n10 < 0.1 * (dis.n11+dis.n10) && comparable >= snp.readIdxs.size()/2){
                    int pos = -1;
                    if (position < ref.size()){
                        pos = position;
                    }
                
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
                    partitions[p].augmentPartition(dis.partition_to_augment, pos);
                    // if (p == 133){
                    //     cout << position <<  " distance with 133: ";
                    //     distancePartition dis2 = distance(partitions[133], snps[position],  ref_base);
                    //     partitions[133].print();
                    //     for (auto r : snps[position].content){
                    //         if (r <= 126){
                    //             cout << r;
                    //         }
                    //         else{
                    //             cout << (char)(r-80);
                    //         }
                    //     }
                    //     cout << endl;
                    //     cout << dis2.n00 << " " << dis2.n11 << " " << dis2.n01 << " " << dis2.n10 << endl;
                    //     cout << float(dis2.n01+dis2.n10)/(min(dis2.n00,dis2.n11) + dis2.n01 + dis2.n10) << endl;
                    // }

                    meanError += float(dis.n01+dis.n10)/(dis.n00 + dis.n11 + dis.n01 + dis.n10);
                    numberOfExtensions += 1;

                    break;
                }
            }

            if (!found && position < ref.size()){    // the second condition is here to create partitions only at specific spots of the backbone
                partitions.push_back(Partition(snp, position, mask, ref_base)); 
                if (partitions.size() == 95){
                    debug_reads_of_interest = snp.readIdxs;
                }
            }
            else{
                position += 5;      //two suspect positions next to each other can be artificially correlated through alignement artefacts
            }
            
            numberOfSuspectPostion += 1;

        }

    }

    /*
    vector<vector<double>> distances (suspiciousColumns.size(), vector<double>(suspiciousColumns.size(), 1));
    for (auto i = 0 ; i < suspiciousColumns.size() ; i++){
        for (auto j = i+1 ; j < suspiciousColumns.size() ; j++){
            if (suspiciousColumns[j].pos - suspiciousColumns[i].pos < 20000){
                distancePartition dis = distance(suspiciousColumns[i], suspiciousColumns[j]);
                distances[i][j] = float(dis.n01+dis.n10)/(dis.n00 + dis.n01 + dis.n10);
                distances[j][i] = distances[i][j];
            }
        }
    }
    //to build the adjacency matrix take the five closest positions to each position
    vector<vector<int>> adjacency (suspiciousColumns.size(), vector<int>(suspiciousColumns.size(), 0));
    for (auto i = 0 ; i < suspiciousColumns.size() ; i++){
        vector<pair<double, int>> distancesToI;
        for (auto j = 0 ; j < suspiciousColumns.size() ; j++){
            if (i != j){
                distancesToI.push_back(make_pair(distances[i][j], j));
            }
        }
        sort(distancesToI.begin(), distancesToI.end());
        int j = 0;
        while (j < min(5, (int) distancesToI.size())  && distancesToI[j].first < 0.1 ){
            adjacency[i][distancesToI[j].second] = 1;
            j++;
        }
    }

    //ouptut the graph
    vector<int> vectorEmpty (suspiciousColumns.size(), 0);
    vector<int> clusters = chinese_whispers(adjacency);
    outputGraph(adjacency, clusters, "tmp/adjacency_matrix.gdf");
    cout << "qfmldjcjc" << endl;

    unordered_map<int, Partition> clustPartitions;
    for (auto i = 0 ; i < suspiciousColumns.size() ; i++){
        if (clustPartitions.find(clusters[i]) == clustPartitions.end()){
            clustPartitions[clusters[i]] = Partition(suspiciousColumns[i], suspiciousColumns[i].pos, mask, suspiciousColumns[i].ref_base);
        }
        else{
            if (clusters[i] == 961){
                cout << " qcccs ";
                print_snp(suspiciousColumns[i], mask);
            }
            clustPartitions[clusters[i]].augmentPartition(suspiciousColumns[i], suspiciousColumns[i].pos);
        }
    }

    for (auto p : clustPartitions){
        if (p.second.number() > 4){
            cout << "partitionds " << p.first << " " << p.second.number() << endl;
            p.second.print();
        }
    }

    exit(1);
    */


    // meanError /= numberOfExtensions+1;

    // // cout << "dodoozozz" << endl;

    // int np = 133;
    // partitions[np].print();
    // for (auto r = 0 ; r < partitions[np].getPartition().size() ; r++){
    //     if ((partitions[np].getPartition()[r] == -1 || partitions[np].getPartition()[r] == 1) && partitions[np].getMore()[r] >= 3){
    //         cout << "fdiuc " <<  partitions[np].getPartition()[r] << " " << partitions[np].getReads()[r] << " : " << partitions[np].getMore()[r] << " " << partitions[np].getLess()[r] << endl;
    //     }
    // }

    // cout << "qlkfjdjql" << endl;
    // int n = 0;
    // for (auto p : partitions){
    //     if (p.number() > 2){
    //         cout << n << " ";
    //         p.print();
    //     }
    //     n++;
    // }

    // if (DEBUG){
    //     cout << "found " << numberOfSuspectPostion << " suspect positions" << endl;
    // }

    if (partitions.size() == 0){ //there are no position of interest
        return vector<Partition>();
    }

    
    float threshold = 4;

    vector<Partition> listOfFinalPartitions;
    for (auto p1 = 0 ; p1 < partitions.size() ; p1++){

        // if (partitions[p1].number() > 0){
        //     cout << "iqdoudofq non informative partition : "  << p1 << " " << snps.size() << " " << partitions[p1].isInformative(true) << endl;
        //     partitions[p1].print();
        // }
        
        if (partitions[p1].number() > threshold && partitions[p1].isInformative(true)){

            bool different = true;
            
            for (auto p2 = 0 ; p2 < listOfFinalPartitions.size() ; p2++){

                distancePartition dis = distance(listOfFinalPartitions[p2], partitions[p1], 2);
                if (dis.augmented 
                    && (dis.n00+dis.n11 > 5*(dis.n01+dis.n10) || dis.n10 + dis.n01 > 5*(dis.n00+dis.n11))
                    && dis.n10 < 2*dis.n01 && dis.n01 < 2*dis.n10){
                    Partition newPart = listOfFinalPartitions[p2];

                    // cout << endl <<  "mergddding : " << endl;
                    // newPart.print();
                    // partitions[p1].print();

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
                partitions[p1].apply_mask(mask);
                listOfFinalPartitions.push_back(partitions[p1]);
            }
        }
    }

    // cout << "jqjlk num of partitions : " << listOfFinalPartitions.size() << endl;

    //select only compatible partitions. In the worst case scenario, trashed partitions will be recovered when masking on the next iteration
    vector<Partition> compatiblePartitions = select_compatible_partitions(listOfFinalPartitions, numberOfReads, meanError);

    vector<size_t> newSuspectPositions;
    int index_previous = 0;
    int index_here = 0;
    for (size_t pos = 0 ; pos < snps.size() ; pos++){
        bool toAdd = false;
        if (index_previous < suspectPostitions.size() && suspectPostitions[index_previous] == pos){
            toAdd = true;
            index_previous += 1;
        }
        if (index_here < suspectPostitionsHere.size() && suspectPostitionsHere[index_here] == pos){
            toAdd = true;
            index_here += 1;
        }
        if (toAdd && (newSuspectPositions.size () == 0 || pos-newSuspectPositions.back()>4)){
            newSuspectPositions.push_back(pos);
        }
    }
    suspectPostitions = newSuspectPositions;

    if (compatiblePartitions.size() > 0 && compatiblePartitions[0].number() < min(min(15.0, 0.01*suspectPostitions.size()), 0.01*(compatiblePartitions[0].get_right()-compatiblePartitions[0].get_left()))){
        return vector<Partition>();
    }

    return compatiblePartitions;
    
}

/**
 * @brief Given positions and snps, tries to solidify the partitions by adding new positions to them
 * 
 * @param ref 
 * @param partitions 
 * @param snps 
 * @param mask 
 * @param suspectPostitions 
 * @param meanError 
 * @param numberOfReads 
 * @return std::vector<Partition> 
 */
std::vector<Partition> second_pass(
    std::string& ref,
    std::vector<Partition> &partitions, 
    std::vector<Column> &snps,
    std::vector<size_t> &suspectPostitions,
    float &meanError,
    int numberOfReads){

    // for each partition, identify a group of confident '1' (and of confident '-1'). At positions where they all agree with consensus, exclude all reads that do not agree with consensus

    vector<Partition> newPartitions;
    for (auto p : partitions){

        for (auto conf : p.getConfidence()){
            cout << int((1-conf)*10) << " ";
        }
        cout << endl;

        auto readsIdx = p.getReads();
        auto part = p.getPartition();
        auto mores = p.getMore();
        vector <float> conf = p.getConfidence();

        for (int side : {1,-1}){

            vector <int> excludedReads (numberOfReads, 0);
            vector <int> includedReads (numberOfReads, 0);

            vector <bool> confidentReads (numberOfReads, false);

            //first find a group of confident reads that are on 'side'

            for (auto pos = 0 ; pos < ref.size() ; pos += 1000){
                //find the three most confident reads at this position
                vector <pair<int, float>> confidentReadsAtPos;
                int indexInSnp = 0;
                for (int i = 0 ; i < readsIdx.size() ; i++){
                    while (indexInSnp < snps[pos].content.size() && snps[pos].readIdxs[indexInSnp] < readsIdx[i]){
                        indexInSnp++;
                    }
                    if (indexInSnp >= snps[pos].content.size()){
                        break;
                    }
                    if (snps[pos].readIdxs[indexInSnp] == readsIdx[i] && part[i] == side && mores[i] > 5){
                        confidentReadsAtPos.push_back({readsIdx[i], conf[i]});
                    }
                }
                sort(confidentReadsAtPos.begin(), confidentReadsAtPos.end(), [](pair<int, float> a, pair<int, float> b){return a.second > b.second;});
                confidentReadsAtPos.resize(min(3, (int)confidentReadsAtPos.size()));
                for (auto r : confidentReadsAtPos){
                    confidentReads[r.first] = true;
                }
            }
            
            cout << "confiiiiiiident reads : " << confidentReads.size() << endl;
            for (auto r = 0 ; r < confidentReads.size() ; r++){
                if (confidentReads[r]){
                    cout << r << " ";
                }
            }
            cout << endl;

            //now go through snps and take all positions that agree with the confident reads
            for (auto position : suspectPostitions){

                //first find the two most frequent alleles at this position
                unordered_map<char, int> alleleCount;
                int indexsnp = 0;
                for (auto r = 0 ; r < readsIdx.size() ; r++){
                    while(indexsnp < snps[position].readIdxs.size() && snps[position].readIdxs[indexsnp] < readsIdx[r]){
                        indexsnp += 1;
                    }
                    if (indexsnp >= snps[position].readIdxs.size()){
                        break;
                    }
                    if (part[r] != -2){
                        char ch = snps[position].content[indexsnp];
                        if (alleleCount.find(ch) == alleleCount.end()){
                            alleleCount[ch] = 0;
                        }
                        alleleCount[ch] += 1;
                    }
                }
                unsigned char refBase = ref[position];
                unsigned char altBase = ' ';
                int maxCount = 0;
                int secondMaxCount = 0;
                for (auto it = alleleCount.begin() ; it != alleleCount.end() ; it++){
                    if (it->first != refBase && it->second >= maxCount){
                        secondMaxCount = maxCount;
                        maxCount = it->second;
                        altBase = it->first;
                    }
                    else if (it->first != refBase && it->second >= secondMaxCount){
                        secondMaxCount = it->second;
                    }
                }
                if (maxCount < secondMaxCount *5 || maxCount < 5){ //this is not so much a biallelic loci
                    continue;
                }

                int not1 = 0;
                int numberOfConfidentReadsHere = 0;
                for (int r = 0 ; r < snps[position].readIdxs.size() ; r++){

                    //looking for positions where all the confident reads agree with consensus
                    if (confidentReads[snps[position].readIdxs[r]]){
                        if (snps[position].content[r] == altBase){
                            not1 += 1;
                        }

                        numberOfConfidentReadsHere += 1;
                    }
                }

                cout << "Posiddtion: " << position << " not1: " << not1 << " numberOfConfidentReadsHere: " << numberOfConfidentReadsHere << " " << refBase << " "<< altBase << endl;
                int n = 0;
                int cri = 0;
                for (unsigned char i : snps[position].content){

                    if (confidentReads[snps[position].readIdxs[n]]){
                        cout << ".";
                    }
                    if (i > 126){
                        cout << (unsigned char) i - 80 << " ";
                    }
                    else if (i < 33){
                        cout << (unsigned char) i + 33 << " ";
                    }
                    else{
                        cout << i << " ";
                    }
                    n += 1;
                }
                cout << endl;

                //now, if all confident reads agree with consensus, exclude all reads that do not agree with consensus
                if (not1 < min(1.0 ,0.1*numberOfConfidentReadsHere) && numberOfConfidentReadsHere >= 3){
                    for (int r = 0 ; r < snps[position].readIdxs.size() ; r++){
                        if (snps[position].content[r] == altBase){
                            excludedReads[snps[position].readIdxs[r]] += 1;
                        }
                        else if (snps[position].content[r] == refBase){
                            includedReads[snps[position].readIdxs[r]] += 1;
                        }
                    }
                }
            }

            cout << "Recomputed a partition: " << side << endl;
            p.print();
            int allEx = 0;
            int allIn = 0;
            cout << "Excluded reads: " << endl;
            for (int i = 0 ; i < excludedReads.size() ; i++){
                if (includedReads[i] > 9){
                    cout << excludedReads[i] << " ";
                }
                else{
                    cout << "  ";
                }
                if (excludedReads[i] < 5 ){
                    allIn += includedReads[i];
                    allEx += excludedReads[i];
                }
            }
            cout << endl;
            cout << "All excluddded reads: " << allEx << " all included reads: " << allIn << endl;

        }


    }
    cout << "fiouqisodpouqfp" << endl;
    exit(1);
    return {};
 
}

/**
 * @brief Create all the masks, each mask isolating a group of reads that is not separated in partitions 
 * 
 * @param partitions 
 * @param numberOfReads 
 * @param all_maks_ever_tested list of maks already tested, so that the same mask are not computed over and over again
 * @return vector<vector<bool>> A vector of masks
 */
vector<vector<bool>> create_masks(vector<Partition> &partitions, int numberOfReads, vector<vector<bool>> &all_maks_ever_tested){

    //for each partition, create a mask that isolates the 1s and the -1s
    // vector<vector<bool>> masks;
    // for (int npart = 0 ; npart < partitions.size() ; npart++){
    //     for (int side : {1, -1}){
    //         vector<bool> mask (numberOfReads, false);
    //         vector<int> reads = partitions[npart].getReads();
    //         vector<short> part = partitions[npart].getPartition();
    //         for (int i = 0 ; i < reads.size() ; i++){
    //             if (part[i] == side){
    //                 mask[reads[i]] = true;
    //             }
    //         }

    //         // test if this mask has already been tested
    //         bool already_tested = false;
    //         for (auto testedMask : all_maks_ever_tested){
    //             int n11 = 0;
    //             int n00 = 0;
    //             int n01 = 0;
    //             int n10 = 0;
    //             for (auto r = 0 ; r < mask.size() ; r++){
    //                 if (mask[r] && testedMask[r]){
    //                     n11 += 1;
    //                 }
    //                 else if (!mask[r] && !testedMask[r]){
    //                     n00 += 1;
    //                 }
    //                 else if (mask[r] && !testedMask[r]){
    //                     n10 += 1;
    //                 }
    //                 else if (!mask[r] && testedMask[r]){
    //                     n01 += 1;
    //                 }
    //             }
    //             // cout << "diddddddfferent: " << different << " equal: " << equal << endl;
    //             if (n01 <= 5){ //the this mask is an extension of the tested mask
    //                 already_tested = true;
    //             }
    //         }
    //         if (!already_tested){
    //             masks.push_back(mask);
    //         }

    //     }
    // }
    // return masks;

    vector<vector<bool>> masks;

    vector<int> readsRepartition (numberOfReads, 0);
    for (int npart = 0 ; npart < partitions.size() ; npart++){
        vector<int> reads = partitions[npart].getReads();
        vector<short> part = partitions[npart].getPartition();
        for (int i = 0 ; i < reads.size() ; i++){
            if (part[i] == 1 || part[i] == -1){
                readsRepartition[reads[i]] += (part[i]+3)/2 * pow(3, npart) ;
            }
        }
    }

    //check what integer is present more than five times in readsRepartition
    std::unordered_map<int, int> count;
    for (int i = 0 ; i < readsRepartition.size() ; i++){
        count[readsRepartition[i]] += 1;
    }

    // cout << "here isee count: " << endl;
    // for (auto cluster : count){
    //     cout << cluster.first << " " << cluster.second << endl;
    // }
    
    for (auto cluster : count){
        if (cluster.second > 5 && cluster.first != 0){ //0 means present on no partitions
            vector<bool> mask(numberOfReads, false);
            for (int i = 0 ; i < readsRepartition.size() ; i++){
                if (readsRepartition[i] == cluster.first){
                    mask[i] = true;
                }
            }

            bool already_tested = false;
            for (auto testedMask : all_maks_ever_tested){
                int n11 = 0;
                int n00 = 0;
                int n01 = 0;
                int n10 = 0;
                for (auto r = 0 ; r < mask.size() ; r++){
                    if (mask[r] && testedMask[r]){
                        n11 += 1;
                    }
                    else if (!mask[r] && !testedMask[r]){
                        n00 += 1;
                    }
                    else if (mask[r] && !testedMask[r]){
                        n10 += 1;
                    }
                    else if (!mask[r] && testedMask[r]){
                        n01 += 1;
                    }
                }
                // cout << "diddddddfferent: " << different << " equal: " << equal << endl;
                if (n01 <= 5){ //the this mask is an extension of the tested mask
                    already_tested = true;
                }
            }
            if (!already_tested){
                masks.push_back(mask);
            }

        }
    }

    return masks;

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
 * @brief Gives a distance between two columns
 * 
 * @param col1 
 * @param col2 
 * @return distancePartition 
 */
distancePartition distance(Column& col1, Column& col2){
    distancePartition res;
    res.n00 = 0;
    res.n01 = 0;
    res.n10 = 0;
    res.n11 = 0;
    res.phased = 1;

    int n1 = 0;
    int n2 = 0;
    auto it1 = col1.readIdxs.begin();
    auto it2 = col2.readIdxs.begin();
    while (it1 != col1.readIdxs.end() && it2 != col2.readIdxs.end()){

        if (*it1 == *it2){
            ++it1;
            ++it2;

            if (col2.content[n2] == col2.ref_base){

                if (col1.content[n1] == col1.ref_base){
                    res.n11 += 1;
                }
                else if (col1.content[n1] == col1.second_base){
                    res.n01 += 1;
                }
            }
            else if (col2.content[n2] == col2.second_base){

                if (col1.content[n1] == col1.ref_base){
                    res.n10 += 1;
                }
                else if (col1.content[n1] == col1.second_base){
                    res.n00 += 1;
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

    return res;
}

//input : a distancePartition 
//output : the chisquare (one degree of freedom)
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

//input : two partitions and thresholds for comparing the partitions
//output : true if the two partitions are the same given the thresholds. In that case, merge partitions into par1
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
        else if (more1[r1] > 4 && more2[r2] > 4){

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
 * @brief Takes all partition and returns a subset of compatible partitions
 * 
 * @param partitions 
 * @param numberOfReads 
 * @param errorRate 
 * @return vector<Partition>  subset of compatible partitions
 */
vector<Partition> select_compatible_partitions(vector<Partition> &partitions, int numberOfReads, float errorRate){
    
    vector<vector<int>> allIdxs;
    for (auto i : partitions){
        allIdxs.push_back(i.getReads());
    }
    vector<vector<short>> allPartitions;
    for (auto i : partitions){
        allPartitions.push_back(i.getPartition());
    }
    vector<vector<float>> allConfidences;
    for (auto i : partitions){
        allConfidences.push_back(i.getConfidence());
    }

    //as a first step, we'll try to throw away the partitions that do not look very sure of themselves
    /*to to that, we'll look at each read, see if it's well clustered in one partition. 
    If yes, all partitions where it's badly clustered are suspicious
    Indeed, we expect the different errors to compensate each other
    */
    

    vector<bool> readsClassified (numberOfReads, true);
    for (int par = 0 ; par < partitions.size() ; par++){
        int c = 0;
        for (auto read : allIdxs[par]){
            if (allPartitions[par][c] != 0 && allConfidences[par][c] >= 0.7){
                readsClassified[read] = true;
            }
            c++;
        }
    }

    //we know what reads are not well classified anywhere.
    //now go through the partitions to discard those that are not compatible with each other. Take only the best one

    //as a first step, we'll compute a score for each partition evaluating on how confident it is and order the partitions by that score   
    for (int par = 0 ; par < partitions.size() ; par++){
        //compute a score evaluating the certainty of the partition
        partitions[par].compute_conf();
    }

    struct {
        bool operator()(Partition a, Partition b) const { return a.get_conf() > b.get_conf(); }
    } customLess;
    std::sort(partitions.begin(), partitions.end(), customLess);

    // if (DEBUG){
    //         cout << "Here are all the partitions, sorted : " << endl;
    //     for (auto p : partitions){
    //         p.print();
    //     }
    // }

    //draw a list of compatible partitions
    
    vector<bool> compatibles(partitions.size(), false);
    for (int p = 0 ; p < partitions.size() ; p++){
        bool compatible = true;
        for (int p2 = 0 ; p2 < p ; p2++){
            if (compatibles[p2]){ //for all more sturdy validated partitions, check if this one is compatible
                int compCode = compatible_partitions(partitions[p], partitions[p2]);
                if (compCode == 0){
                    compatible = false;
                    break;
                }
                else if (compCode == 2){ //means that they're the same partition
                    compatible = false;
                    partitions[p2].mergePartition(partitions[p]);
                    break;
                }
            }
        }
        if (compatible){
            compatibles[p] = true;
        }
    }
    
    vector<Partition> compatiblePartitions;
    for (auto p = 0 ; p < partitions.size(); p++){
        if (compatibles[p]){
            compatiblePartitions.push_back(partitions[p]);
            // if (DEBUG){partitions[p].print();}
        }
    }

    return compatiblePartitions;
}

/**
 * @brief Point out the partitions that are very sure of themselves and that should be kept
 * 
 * @param partitions list of candidate partitions
 * @param trimmedListOfFinalPartitionBool true if a partition is kept, false otherwise. May be initialized with some pre-selected partitions
 * @param numberOfReads 
 * @param errorRate 
 * @return vector<Partition> All the confident partitions
 */
vector<Partition> select_confident_partitions(vector<Partition> &partitions, std::vector<bool> trimmedListOfFinalPartitionBool, int numberOfReads, float errorRate, int numberOfSuspectPositions){

    vector<vector<int>> allIdxs;
    for (auto i : partitions){
        allIdxs.push_back(i.getReads());
    }
    vector<vector<short>> allPartitions;
    for (auto i : partitions){
        allPartitions.push_back(i.getPartition());
    }
    vector<vector<float>> allConfidences;
    for (auto i : partitions){
        allConfidences.push_back(i.getConfidence());
    }

    //first see which reads are well classified in at least one partition
    vector<bool> readsClassified (numberOfReads, true);
    for (int par = 0 ; par < partitions.size() ; par++){
        int c = 0;
        for (auto read : allIdxs[par]){
            if (allPartitions[par][c] != 0 && allConfidences[par][c] >= 0.7){
                readsClassified[read] = true;
            }
            c++;
        }
    }

    for (int par = 0 ; par < partitions.size() ; par++){

        //because of the ref, the two haplotypes are not strictly equivalent : try to counterbalance that
        int numberOfPartitionnedRead = 0;
        float means [2] = {0,0}; //it's the mean confidence we have in each haplotype
        float n1 = 0.001; //not 0 just to be sure never to divide by 0
        float n0 = 0.001;
        int c = 0;
        for (auto read : allIdxs[par]){
            if (allPartitions[par][c] != 0){
                numberOfPartitionnedRead++;
                if (allPartitions[par][c] == 1){
                    means[1] += allConfidences[par][c];
                    n1++;
                }
                else if (allPartitions[par][c] == -1){
                    means[0] += allConfidences[par][c];
                    n0++;
                }
            }
            c++;
        }

        //compute the confidence level we must have to be sure (different for both haplotypes because of the bias)
        double errors [2] = {min((1-means[0]/n0)/(2-means[0]/n0-means[1]/n1) * errorRate*4 + errorRate, float(0.3))
                            , min((1-means[1]/n1)/(2-means[0]/n0-means[1]/n1) * errorRate*4 + errorRate, float(0.3)) }; //tolerate on average errorRate*3
        // cout << "the centers are : " << means[0]/n0 << " and " << means[1]/n1 << endl;

        int numberOfUnsureReads0 = 0;
        int numberOfUnsureReads1 = 0;
        c = 0;
        for (auto read : allIdxs[par]){
            if (allPartitions[par][c] == 1 && allConfidences[par][c] < 1-errors[1] && readsClassified[read]){ //oops
                numberOfUnsureReads1++;
            }
            else if (allPartitions[par][c] == -1 && allConfidences[par][c] < 1-errors[0] && readsClassified[read]){ //oops
                numberOfUnsureReads0++;
            }
            c++;
        }
        int numberOfUnsureReads = numberOfUnsureReads0+numberOfUnsureReads1;

        if (DEBUG){
            cout << "filtering partition, tolerating " << errors[0] << "," << errors[1] << endl;
            partitions[par].print();
            cout << "There are " << numberOfUnsureReads << " unsure reads, (" << numberOfUnsureReads0 << "+" << numberOfUnsureReads1 << ") out of " 
                << numberOfPartitionnedRead << " partitionned reads, and in terms of size, the het rate is " <<
            partitions[par].number() << " for a length of " << (partitions[par].get_right()-partitions[par].get_left());
        }

        if (float(numberOfUnsureReads)/numberOfPartitionnedRead < 0.15 //then the partition is sure enough of itself 
            && partitions[par].number() > min(20.0, 0.05*numberOfSuspectPositions)){ //and big enough //max(10.0, min(50.0,float(numberOfUnsureReads+1)/numberOfPartitionnedRead/0.15*0.01*(partitions[par].get_right()-partitions[par].get_left()))
            trimmedListOfFinalPartitionBool[par] = true;
            numberOfSuspectPositions -= partitions[par].number();
        }
        else if (float(numberOfUnsureReads)/numberOfPartitionnedRead > 0.2 && (par>0||partitions[par].number() < 30)){ //if it's not the best partition it may be a weird version of the best partition
            trimmedListOfFinalPartitionBool[par] = false;
        }
        else if (float(numberOfUnsureReads)/numberOfPartitionnedRead > 0.3){ //a partition could have been marked as true if it augmented enough when recomputed, but more than 30% is too much
            trimmedListOfFinalPartitionBool[par] = false;
        }
        if(DEBUG){cout << ". Overall, do I take it : " << trimmedListOfFinalPartitionBool[par] << endl;}

    }

    vector<Partition> trimmedListOfFinalPartition;
    for (auto p = 0 ; p < partitions.size() ; p++){
        if (trimmedListOfFinalPartitionBool[p]){
            trimmedListOfFinalPartition.push_back(partitions[p]);
            // cout << "remaining partition : " << endl;
            // listOfFinalPartitions[p].print();
        }
    }

    return trimmedListOfFinalPartition;
}

//input : two partitions
//output : are those two partitions compatible ? : 0 No, 1 Yes, 2 they're the same partition !
int compatible_partitions(Partition &p1 , Partition &p2){

    if (p1.get_right() <= p2.get_left() || p2.get_right() <= p1.get_left()){ //the two partitions are not defined on the same stretch
        return true;
    }

    //compatibility is defined as : either the 0s or the 1s of the extension all fall squarely within one already defined haplotype

    auto idxs = p1.getReads();
    auto content = p1.getPartition();
    auto confidence1 = p1.getConfidence();

    auto idxs2 = p2.getReads();
    auto content2 = p2.getPartition();
    auto confidence2 = p2.getConfidence();

    int n = 0;
    float numberOf1s = 0;
    float numberOf0s = 0;

    int one_in_1 = 0;
    int one_in_0 = 0;
    int zero_in_1 = 0;
    int zero_in_0 = 0;

    int unsureCalls = 0;
    int confidentCalls = 0;

    int n2 = 0;
    for (auto idx1 : idxs){

        while(n2 < idxs2.size() && idxs2[n2] < idx1){
            n2++;
        }

        if (idx1 == idxs2[n2]){

            if (content[n]==1){
                if (content2[n2] == 1){
                    one_in_1 += 1;
                    numberOf1s++;
                }
                else if (content2[n2] == -1){
                    one_in_0 += 1;
                    numberOf1s++;
                }                
            }
            else if (content[n]==-1){
                if (content2[n2] == 1){
                    zero_in_1 += 1;
                    numberOf0s++;
                }
                else if (content2[n2] == -1){
                    zero_in_0 += 1;
                    numberOf0s++;
                } 
            }
        }
        n++;
    }

    //now the idea is to see if the two partitions are compatible
    //additionally, ensure that the most specific parts are made up of ones
    int compatible = 0;

    if (zero_in_0/numberOf0s > 0.8 && zero_in_1/(numberOf0s+numberOf1s) < 0.1 ){
        // p1.flipPartition();
        compatible += 1;
    }
    else if(zero_in_1/numberOf0s > 0.8 && zero_in_0/(numberOf0s+numberOf1s) < 0.1){
        // p1.flipPartition();
        // p2.flipPartition();
        compatible += 1;
    }

    if (one_in_0/numberOf1s > 0.8 && one_in_1/(numberOf0s+numberOf1s) < 0.1){
        compatible += 1;
    }
    else if (one_in_1/numberOf1s > 0.8 && one_in_0/(numberOf0s+numberOf1s) < 0.1){
        // p2.flipPartition();
        compatible += 1;
    }

    // if(DEBUG){cout << "Compatibility : " << numberOf0s << " " << zero_in_0 << " " << zero_in_1 << " ; " << numberOf1s << " " << one_in_0 << " " << one_in_1 << " ; " << compatible << endl;}
    return compatible;
}

/**
 * @brief Concatenates binary partitions to form a single multi-part partition
 * 
 * @param listOfFinalPartitions the binary partitions
 * @param numberOfReads 
 * @param columnThere This is a column of SNPs matrix, to know what reads are present in the interval
 * @return vector<int> the final partition
 */
vector<int> threadHaplotypes_in_interval(vector<Partition> &listOfFinalPartitions, int numberOfReads){

    if (listOfFinalPartitions.size() == 0){
        //return a vector composed of 0s for reads aligning and -1 for reads not aligning (so that output_GAF does not believe that all the reads are aligning)
        vector<int> unpartition (numberOfReads, 0);
        return unpartition;
    }

    vector<vector<short>> allPartitions;
    for (auto i : listOfFinalPartitions){
        allPartitions.push_back(i.getPartition());
    }
    vector<vector<float>> allConfidences;
    for (auto i : listOfFinalPartitions){
        allConfidences.push_back(i.getConfidence());
    }

    vector<vector<int>> allMores;
    for (auto i : listOfFinalPartitions){
        allMores.push_back(i.getMore());
    }

    int numberOfAssignedReads=0;
    //inventroriate the reads that are present on all partitions
    vector<int> numberOfTimeOfEachRead (numberOfReads, 0);
    vector<bool> notAllMasked (numberOfReads, false);
    for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
        int r = 0;
        for (auto read : listOfFinalPartitions[binary].getReads()){
            if (allPartitions[binary][r] != 0){
                numberOfTimeOfEachRead[read] += 1;
            }
            if (allPartitions[binary][r] == 1 || allPartitions[binary][r] == -1){
                notAllMasked[read] = true;
            }
            r += 1;
        }
    }
    //now we have the number of times each read is present in the partitions, create a vector of reads that are present in all partitions
    vector<int> readsPresentInAllPartitions;
    for (auto i = 0 ; i < numberOfReads ; i++){
        if (numberOfTimeOfEachRead[i] == listOfFinalPartitions.size()){
            readsPresentInAllPartitions.push_back(i);
        }
    }
    //print readsPresentInAllPartitions

    vector<long double> clustersD (numberOfReads, -1); 
    vector<long double> clustersAll (numberOfReads, -1); 
    vector<int> presentReads (numberOfReads, listOfFinalPartitions.size()); //this will count to 0 the number of partitions in which the read is present

    for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
        int c = 0;

        //for (auto read : readsPresentInAllPartitions){
        for (auto read : listOfFinalPartitions[binary].getReads()){

            // if (c == listOfFinalPartitions[binary].getReads().size()-1){
            //     cout << "last aiep read: " << numberOfReads << " " << binary << " " << read << " " << numberOfTimeOfEachRead[read] << " " << listOfFinalPartitions.size() << endl;
            // }
            if (numberOfTimeOfEachRead[read] == listOfFinalPartitions.size()){

                presentReads[read] -= 1;

                auto camp = allPartitions[binary][c];
                if (camp == -2){ //mark all the masked reads as belogning to an arbitrary part
                    camp = 1;
                }
                if (clustersAll[read] == -1 && camp != 0){ //means that the read is present in one partition
                    clustersD[read] = 0;
                    clustersAll[read] = 0;
                }

                //only keep a read if it covers at least 50% of the suspect positions
                if ( allConfidences[binary][c] < 0.8 && allMores[binary][c] < 0.5*listOfFinalPartitions[binary].number() /*|| camp == 0 /*|| allMores[binary][c] < 5*/){  //0.7 to be pretty confident about the reads we separate (more than 70% on all partitions)
                    clustersD[read] = -10; //definitely putting that below -1
                }
                else if (camp != 0 && clustersD[read] != -10){
                    clustersD[read] += pow(2, binary)*int(0.5+0.5*camp); //*0 if haplotype -1, *1 if haplotype 1
                }

                if (camp!=0){
                    clustersAll[read] += pow(2, binary)*int(0.5+0.5*camp); //*0 if haplotype -1, *1 if haplotype 1
                }
            }
            c += 1;
        }

    }

    // cout << "clustszeD: " << endl;
    // int nr = 0;
    // for (auto i : clustersD){
    //     if (i == 8191){
    //         cout << " ";
    //     }
    //     else{
    //         cout << i << ",";
    //     }
    // }
    // cout << endl;

    // cout << "clustersall :" << endl;
    // for (auto i : clustersAll){
    //     if (i >= 0){
    //         cout << i << " ";
    //     }
    // }
    // cout << endl;

    set<int> count; //a set of all existing groups
    // vector<int> frequenceOfPart (pow(2, listOfFinalPartitions.size())); 
    unordered_map<int, int> frequenceOfPart;
    unordered_map<long double, int> doubleToIdx;
    int idx = 0;
    vector <int> clusters (numberOfReads, -1);
    int n = 0;
    for (auto id : clustersD){
        if (id >= 0 && notAllMasked[n]){
            if (doubleToIdx.find(id) == doubleToIdx.end()){
                frequenceOfPart[idx] = 0;
                doubleToIdx[id] = idx;
                count.emplace(idx);
                idx += 1;
            }
            frequenceOfPart[doubleToIdx[id]] += 1;
            numberOfAssignedReads ++;
            clusters[n] = doubleToIdx[id];
        }
        else {
            clusters[n] = -1;
        }
        n += 1;
    }

    
    /*
        As a first approximation:
        to make sure we don't over-estimate the number of clusters, we'll make the assumption that all haplotypes have haplotype-specific mutation => the number of final cluster cannot be higher than the number of partitions
        Then:
        also count the number of solid clusters (i.e. frequent)
    */
    // vector<float> proportionOf0;
    // for (auto p : listOfFinalPartitions){
    //     proportionOf0.push_back(p.numberOf0()/numberOfAssignedReads);  
    // }
    // vector<float> listOfLikelihoods;
    // vector<float> listOfProbas;
    // size_t max_number_of_cluster = 0; //count the number of clusters that look solid
    // for (auto group : count){
    //     int id = group;
    //     double proba = 1;
    //     for (int binary = listOfFinalPartitions.size()-1 ; binary > -1 ; binary--){
    //         if (id%2 == 1){
    //             proba *= 1-proportionOf0[binary];
    //         }
    //         else if (id%2==0){
    //             proba *= proportionOf0[binary];
    //         }
    //         id = int(id / 2);
    //     }
    //     if (DEBUG){
    //         cout << "group " << group << " is expected " << numberOfAssignedReads*proba << " times and arrives " << frequenceOfPart[group] << " times " << endl;
    //     }
    //     listOfLikelihoods.push_back(float(frequenceOfPart[group]-numberOfAssignedReads*proba) / frequenceOfPart[group]);
    //     listOfProbas.push_back(proba);
    //     if (frequenceOfPart[group] > 4 && numberOfAssignedReads*proba < frequenceOfPart[group]*3){
    //         max_number_of_cluster++;
    //     }
    //     // cout << "group " << group << " likelihood: " << float(frequenceOfPart[group] -numberOfAssignedReads*proba) / frequenceOfPart[group] << endl;
    // }

    // float minLikelihood = -1000;
    // if (listOfLikelihoods.size() > listOfFinalPartitions.size() && listOfFinalPartitions.size()>1){
    //     int extra = 0; //when there is only 1 partitions there is actually 2 parts
    //     if (listOfFinalPartitions.size()<=2){extra=1;}
    //     vector <float> minElements (listOfLikelihoods.size());
    //     std::partial_sort_copy(listOfLikelihoods.begin(),  listOfLikelihoods.end(), minElements.begin(), minElements.end());
    //     if (minElements.size() >= listOfFinalPartitions.size()+extra){
    //         auto number_of_cluster = max(max_number_of_cluster, listOfFinalPartitions.size()+extra);
    //         cout << "qipdiqp number of clusters " << number_of_cluster << endl;
    //         minLikelihood = minElements[minElements.size()-number_of_cluster];
    //     }
    //     else{
    //         minLikelihood = minElements[0];
    //     }
    // }

    //establish a list of possible clusters
    std::set <int> listOfGroups;
    int indexOfCount = 0;

    int numberOfAssignedReads2 = numberOfAssignedReads;
    float newtotalproba = 1;
    for (auto group : count){
        if (/*listOfLikelihoods[indexOfCount] >= minLikelihood &&*/ frequenceOfPart[group] >= 5){ //>5 because we want at least 5 reads per haplotype
            listOfGroups.emplace(group);
            numberOfAssignedReads2 -= frequenceOfPart[group];
            // newtotalproba -= listOfProbas[indexOfCount];
        }
        indexOfCount++;
    }

    //switch all the clusters present less than 5 times to -1
    for (int i = 0 ; i < clusters.size() ; i++){
        if (clusters[i] >= 0 && listOfGroups.find(clusters[i]) == listOfGroups.end()){
            clusters[i] = -1;
        }
    }
    //now maybe rescue cluster that were not 100% sure but in hindsight seem good
    // indexOfCount = 0;
    // for (auto group : count){
    //     float newproba = listOfProbas[indexOfCount]/newtotalproba;
    //     cout << "bbcxv trying to see if should rescue " << group << " is now expected " << numberOfAssignedReads2*newproba << " times and arrives " << frequenceOfPart[group] << " times " << endl;
    //     if (listOfGroups.find(group) == listOfGroups.end() 
    //         && (numberOfAssignedReads2*newproba < 0.67*frequenceOfPart[group] || (numberOfAssignedReads2-frequenceOfPart[group]) < 0.5*(numberOfAssignedReads2-numberOfAssignedReads2*newproba)) 
    //         && newproba < 0.9 && frequenceOfPart[group] > 5){ //>5 because we want at least 5 reads per haplotype
    //         listOfGroups.emplace(group);
    //         cout << "rescuing " << group << endl;
    //     }
    //     indexOfCount++;
    // }

    // cout << "posssippoble clusters : "; for (auto g : listOfGroups){cout << g << " ";} cout << endl;
    // cout << "Ryes : ";
    // for (auto i : clusters){
    //     cout << i << ",";
    // }
    // cout << endl;

    //if there is only one cluster, return that
    if (listOfGroups.size() <= 1){
        return vector<int> (clusters.size(), 1);
    }

    //now most reads should be assigned to a cluster. Rescue those that have not been assigned or assigned to a rare cluster

    // for (auto read = 0 ; read < clusters.size() ; read++){

    //     if (listOfGroups.find(clusters[read]) == listOfGroups.end() && presentReads[read] == 0 && notAllMasked[read] && clusters[read] == -1){ //this means the read needs to be rescued 

    //         //cout << "rescuing "; for(auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++) {cout << allPartitions[binary][read];} cout << endl;
    //         //iterate through the list of groups and choose the one that can be explained by the less mistakes
    //         int bestGroup = -1;
    //         double bestGroupScore = 0;

    //         //compute a score for assigning the read to each existing group, then assign the read to the best group
    //         for (auto group : listOfGroups){

    //             double score = 0;
    //             int groupDecomposition = group;
    //             long long int readDecomposition = clustersAll[read];

    //             for (int binary = listOfFinalPartitions.size() -1 ; binary > -1  ; binary--){
    //                 int expectedAssignation = groupDecomposition%2;
    //                 int actualAssignation = readDecomposition%2;

    //                 if (actualAssignation != expectedAssignation && allPartitions[binary][read] != 0){
    //                     score -= allConfidences[binary][read]-0.5;
    //                 }
    //                 else if (actualAssignation == expectedAssignation){
    //                     score += allConfidences[binary][read]-0.5;
    //                 }

    //                 groupDecomposition /= 2;
    //                 readDecomposition /= 2;
    //             }
    //             if (score > bestGroupScore){
    //                 bestGroup = group;
    //                 bestGroupScore = score;
    //             }
    //         }
    //         clusters[read] = bestGroup;

    //         // cout << "Rescuing ";
    //         // for (int binary = 0 ; binary < listOfFinalPartitions.size()  ; binary++){
    //         //     cout << listOfFinalPartitions[binary].getMore()[read]<< "/" << listOfFinalPartitions[binary].getLess()[read] << " ";
    //         // }
    //         // cout << " as : " << bestGroup << endl;
            
    //     }
    // }

    // cout << "Res of threazd haplotypes: ";
    // for (auto i : clusters){
    //     if (i > -1){
    //         cout << i << " ";
    //     }
    // }
    // cout << endl;

    return clusters;

}


/**
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
    int chunk, 
    std::vector<size_t> &suspectPostitions,
    std::vector<std::vector<int>> &adjacencyMatrix,
    int sizeOfWindow){

    set<int> listOfGroups;
    int max_cluster = 0;
    unordered_map<int, int> indexOfGroups;

    int index = 0;
    for (auto read = 0 ; read < clusteredReads.size() ; read++){
        if (clusteredReads[read] > -1){
            listOfGroups.emplace(clusteredReads[read]);
            if (indexOfGroups.find(clusteredReads[read]) == indexOfGroups.end()){
                indexOfGroups[clusteredReads[read]] = index;
                index++;
            }
        }
    }
    vector<vector<bool>> imcompatibilities (listOfGroups.size(), vector<bool> (listOfGroups.size(), false));

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

    for (auto position : suspectPostitions){
        if (position >= chunk*sizeOfWindow && position < (chunk+1)*sizeOfWindow){

            //what bases occur in which cluster ?
            unordered_map<int,unsigned char> cluster_to_majority_base;
            unordered_map<int,unordered_map<unsigned char,int>> bases_in_each_cluster;
            unordered_map<int,int> nb_bases_in_each_cluster;
            unordered_map<unsigned char,int> bases_in_total;
            for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
                int read = snps[position].readIdxs[r];
                unsigned char base = snps[position].content[r];
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
                    if (cluster_to_majority_base[group1] != ' ' && cluster_to_majority_base[group2] != ' '){
                        if (cluster_to_majority_base[group1] != cluster_to_majority_base[group2]){
                            imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]] = true;
                            imcompatibilities[indexOfGroups[group2]][indexOfGroups[group1]] = true;

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
    //         if (imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]]){
    //             cout << "incompatiebility btw " << group1 << " and " << group2 << endl;
    //         }
    //     }
    //     cout << endl;
    // }

    // list all the distances between different clusters
    std::map < pair<int, int>, double > number_of_links_between_the_two_clusters;
    unordered_map <int,int> number_of_links_in_each_cluster;
    for (auto read1 = 0 ; read1 < adjacencyMatrix.size() ; read1++){
        for (auto read2 = 0 ; read2 < adjacencyMatrix.size() ; read2++){
            if (adjacencyMatrix[read1][read2]){
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


    //now that all the incompatibilities are clear, we can merge the clusters
    unordered_map <int, int> old_group_to_new_group;
    for (auto group : listOfGroups){
        old_group_to_new_group[group] = group;
    }
    old_group_to_new_group[-1] = -1;
    old_group_to_new_group[-2] = -2;

    for (auto pair_of_clusters : sorted_links){
        if (pair_of_clusters.second > 0.05){
            int cluster1 = pair_of_clusters.first.first;
            int cluster2 = pair_of_clusters.first.second;

            if (old_group_to_new_group[cluster1] == old_group_to_new_group[cluster2]){
                continue;
            }
            //check if there are no incompatibilities between the two clusters
            bool incompatibility = false;
            for (auto group1 : listOfGroups){
                if (old_group_to_new_group[group1] == old_group_to_new_group[cluster1]){
                    for (auto group2 : listOfGroups){
                        if (old_group_to_new_group[group2] == old_group_to_new_group[cluster2]){
                            if (imcompatibilities[indexOfGroups[group1]][indexOfGroups[group2]]){
                                // cout << "incompatiebility btw " << group1 << " and " << group2 << endl;
                                // cout << "merging " << cluster1 << " and " << cluster2 << " would create an incompatibility" << endl;
                                incompatibility = true;
                            }
                        }
                    }
                }
            }
            if (!incompatibility){
                // cout << "merging " << cluster1 << " and " << cluster2 << endl;
                for (auto group2 : listOfGroups){
                    if (old_group_to_new_group[group2] == old_group_to_new_group[cluster2]){
                        old_group_to_new_group[group2] = old_group_to_new_group[cluster1];
                    }
                }
            }
        }
    }

    //re-nomber the clusters
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


/**
 * @brief Reclassifies all the read to rescue some reads
 * 
 * @param clusteredReads First separation of the reads
 * @param snps Vector of columns containing the snps
 * @param chunk To look only at the interesting snps
 * @param suspectPostitions To look only at the interesting snps
 * @param sizeOfWindow To look only at the interesting snps
 * @return std::vector<int> 
 */
std::vector<int> rescue_reads(
    std::vector<int> &clusteredReads, 
    std::vector<Column> &snps, 
    int chunk, 
    std::vector<size_t> &suspectPostitions,
    int sizeOfWindow){

    set<int> listOfGroups;
    unordered_map<int, int> sizeOfGroups;
    int max_cluster = 0;
    for (auto read = 0 ; read < clusteredReads.size() ; read++){
        if (clusteredReads[read] > -1){
            listOfGroups.emplace(clusteredReads[read]);
            if (sizeOfGroups.find(clusteredReads[read]) == sizeOfGroups.end()){
                sizeOfGroups[clusteredReads[read]] = 0;
            }
            sizeOfGroups[clusteredReads[read]]++;
        }
    }

    if (listOfGroups.size() <= 1){
        return clusteredReads;
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

    unordered_map<int, int> numberOfSnpsOfGroups;

    for (auto position : suspectPostitions){
        if (position >= chunk*sizeOfWindow && position < (chunk+1)*sizeOfWindow){

            //what bases occur in which cluster ?
            unordered_map<int,unsigned char> cluster_to_majority_base;
            unordered_map<int,unordered_map<unsigned char,int>> bases_in_each_cluster;
            unordered_map<unsigned char,int> bases_in_total;
            for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
                int read = snps[position].readIdxs[r];
                unsigned char base = snps[position].content[r];
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
                }
            }
            
            //count the majority base of each cluster at this position
            set <unsigned char> maxbases;
            for (auto cluster : bases_in_each_cluster){
                int max = min(3.0, 0.5*sizeOfGroups[cluster.first]);
                char maxBase = ' ';
                for (auto base : cluster.second){
                    if (base.second > max){
                        max = base.second;
                        maxBase = base.first;
                    }
                }
                cluster_to_majority_base[cluster.first] = maxBase;
                if (maxBase != ' '){
                    if (max > 5){
                        maxbases.emplace(maxBase);
                    }
                    numberOfSnpsOfGroups[cluster.first]++;
                }
            }
            
            if (maxbases.size() <= 1){
                continue;
            }

            //fill with_which_clusters_goes_this_read
            for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
                int read = snps[position].readIdxs[r];
                unsigned char base = snps[position].content[r];
                int cluster = clusteredReads[read];

                for (int group : listOfGroups){
                    if (cluster_to_majority_base[group] != base && bases_in_total[base] >= 5 && cluster_to_majority_base[group] != ' ' && base != ' '){
                        with_which_clusters_goes_this_read[read][group]++;
                    }
                }
                if (bases_in_total[base] >= 5){
                    on_how_many_snps_is_this_read_defined[read]++;
                }
            }

            // int nr = 0;
            // cout << "eexxww " << position << " " << endl;
            // for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
            //     int read = snps[position].readIdxs[r];
            //     while (nr < read){
            //         if (clusteredReads[nr] >= 0){
            //             cout << " ";
            //         }
            //         nr += 1;
            //     }
            //     if (clusteredReads[nr] >= 0){
            //         unsigned char c = snps[position].content[r];
            //         if (c > 126){
            //             cout << (unsigned char) (c - 80);
            //         }
            //         else{
            //             cout << (unsigned char) c;
            //         }
            //     }
            //     nr += 1;
            // }
            // cout << endl;
            // for (auto r = 0 ; r < clusteredReads.size() ; r++){
            //     if (clusteredReads[r] >= 0){
            //         cout << clusteredReads[r] << " , ";
            //     }
            // }
            // cout << endl;

            // //print cluster_to_majority_base
            // for (auto c : cluster_to_majority_base){
            //     if (c.second < 126){
            //         cout << c.first << " : " << (unsigned char) c.second << "  -  ";
            //     }
            //     else{
            //         cout << c.first << " : " << (unsigned char) (c.second - 80) << "  -  ";
            //     }
            // }
            // cout << endl;

            // if (chunk == 0 && cluster_to_majority_base.find(10) != cluster_to_majority_base.end() && cluster_to_majority_base.find(9) != cluster_to_majority_base.end()){
            //     // cout << "posddition " << position << " : " << (char) cluster_to_majority_base[3] << " " << (char) cluster_to_majority_base[2] << endl;
            //     if (cluster_to_majority_base[10] != cluster_to_majority_base[9]){

            //         cout << "posddisqdodlmqjmtion " << position << " : ";

            //         cout << "bases_in_each_cluster : " << endl;
            //         for (auto c : bases_in_each_cluster){
            //             cout << c.first << " : ";
            //             for (auto b : c.second){
            //                 cout << (char) b.first << " " << b.second << " ";
            //             }
            //             cout << endl;
            //         }

            //         // cout << "skljiuzau clusters: " << endl;
            //         // int nc = 0;
            //         // for (auto c : clusteredReads){
            //         //     if (c == 9 || c == 10){
            //         //         cout << nc << " " <<  c << endl;
            //         //     }
            //         //     nc++;
            //         // } 
            //         // cout << endl;

            //         cout << "eeddzxxww ";
            //         for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
            //             int nr2 = snps[position].readIdxs[r];
            //             if (clusteredReads[nr2] == 9){
            //                 cout << nr2 << " . ";
            //             }
            //             if (clusteredReads[nr2] == 10){
            //                 cout << nr2 << " , ";
            //             }
            //             if (clusteredReads[nr2] == 9 || clusteredReads[nr2] == 10){
            //                 unsigned char c = snps[position].content[r];
            //                 if (c > 126){
            //                     cout << (unsigned char) (c - 80);
            //                 }
            //                 else{
            //                     cout << (unsigned char) c;
            //                 }
            //                 cout << endl;
            //             }
            //             else{
            //                 cout << " ";
            //             }
            //         }
            //         cout << endl;
            //         exit(1);
            //     }
            // }

        }
    }

    //see if some groups are just not well defined
    int max_number_of_snps = 0;
    for (auto g : numberOfSnpsOfGroups){
        if (g.second > max_number_of_snps){
            max_number_of_snps = g.second;
        }
    }
    set<int> to_erase;
    for (auto g : numberOfSnpsOfGroups){
        if (g.second < max_number_of_snps * 0.5){
            to_erase.emplace(g.first);
        }
    }
    for (auto r = 0 ; r < clusteredReads.size() ; r++){
        if (to_erase.find(clusteredReads[r]) != to_erase.end()){
            clusteredReads[r] = -1;
        }
    }
    for (unordered_map<int, int> g : with_which_clusters_goes_this_read){
        for (auto gg : to_erase){
            g.erase(gg);
        }
    }
    for (auto g : to_erase){
        numberOfSnpsOfGroups.erase(g);
    }


    //count how many times a group of reads is exluded from its own cluster (because of this errors)
    unordered_map <int , vector<double>> how_many_times_a_group_is_excluded_from_its_own_cluster;
    unordered_map <int , vector<double>> how_many_times_a_group_is_excluded_from_different_cluster;
    int read = 0;
    int total_number_of_reads = 0;
    bool add = false;
    for (auto read_exclusions : with_which_clusters_goes_this_read){
        int group_number_of_exclusions = 0;
        for (auto group : listOfGroups){
            if (clusteredReads[read] == group){
                if (how_many_times_a_group_is_excluded_from_its_own_cluster.find(group) == how_many_times_a_group_is_excluded_from_its_own_cluster.end()){
                    how_many_times_a_group_is_excluded_from_its_own_cluster[group] = vector<double>();
                }
                group_number_of_exclusions = read_exclusions[group];
            }
            else if (numberOfSnpsOfGroups[group] > 0 && listOfGroups.find(clusteredReads[read]) != listOfGroups.end()){
                if (how_many_times_a_group_is_excluded_from_different_cluster.find(group) == how_many_times_a_group_is_excluded_from_different_cluster.end()){
                    how_many_times_a_group_is_excluded_from_different_cluster[group] = vector<double>();
                }
                how_many_times_a_group_is_excluded_from_different_cluster[group].push_back(double(read_exclusions[group])/numberOfSnpsOfGroups[group]);
            }
        }
        if (listOfGroups.find(clusteredReads[read]) != listOfGroups.end()){
            how_many_times_a_group_is_excluded_from_its_own_cluster[clusteredReads[read]].push_back(double(group_number_of_exclusions)/numberOfSnpsOfGroups[clusteredReads[read]]);
        }
        read++;
    }

    unordered_map <int, pair<double, double>> confidence_level_clusters; // thresholds below which 90% of the cluster is and above which 90% of non-cluster are
    unordered_map <int, pair<double, double>> mean_confidence_level_clusters; // mean confidence level if in cluster or out of cluster
    for (auto group : listOfGroups){
        if (how_many_times_a_group_is_excluded_from_its_own_cluster.find(group) != how_many_times_a_group_is_excluded_from_its_own_cluster.end()){
            int sum = 0;
            int total = 0;
            //find the first decile of how many times a group is excluded from its own cluster
            sort(how_many_times_a_group_is_excluded_from_its_own_cluster[group].begin(), how_many_times_a_group_is_excluded_from_its_own_cluster[group].end(), std::greater<double>());
            int decile = max(int(how_many_times_a_group_is_excluded_from_its_own_cluster[group].size()/10),1);
            double threshold1 = 0;
            double mean1 = 0;
            if (how_many_times_a_group_is_excluded_from_its_own_cluster[group].size() > 0){
                threshold1 = how_many_times_a_group_is_excluded_from_its_own_cluster[group][decile];
                mean1 = how_many_times_a_group_is_excluded_from_its_own_cluster[group][int(how_many_times_a_group_is_excluded_from_its_own_cluster[group].size()/2)];
            }

            //find the last decile of how many times a group is excluded from a different cluster
            sort(how_many_times_a_group_is_excluded_from_different_cluster[group].begin(), how_many_times_a_group_is_excluded_from_different_cluster[group].end(), std::less<double>());
            int decile2 = max(int(how_many_times_a_group_is_excluded_from_different_cluster[group].size()/20),1);
            double threshold2 = 0;
            double mean2 = 0;
            if (how_many_times_a_group_is_excluded_from_different_cluster[group].size() > 0){
                threshold2 = how_many_times_a_group_is_excluded_from_different_cluster[group][decile2];
                mean2 = how_many_times_a_group_is_excluded_from_different_cluster[group][int(how_many_times_a_group_is_excluded_from_different_cluster[group].size()/2)];
            }

            // cout << "gzgrroup " << group << " is excluded " << how_many_times_a_group_is_excluded_from_its_own_cluster[group][decile] << " times while accepted " <<
            //     how_many_times_a_group_is_excluded_from_different_cluster[group][decile2] << endl;
            
            //confidence level of the cluster
            confidence_level_clusters[group] = make_pair(threshold1, threshold2);
            mean_confidence_level_clusters[group] = make_pair(mean1, mean2);
        }
    }

    //go through the reads and classify them
    vector <int> haplotypes = clusteredReads;
    for (auto read = 0 ; read < clusteredReads.size() ; read++){
        if (clusteredReads[read] == -1){
            unordered_map<int, bool> in_this_cluster;
            for (auto group : listOfGroups){
                in_this_cluster[group] = false;
            }
            int clusterLessExcluded = 0;
            double minExclusion = 1;
            double secondMinExclusion = 1;
            for (auto g : with_which_clusters_goes_this_read[read]){ //g.first: group ; g.second : number of times this group is excluded from this read

                double exclusionLevel = double(g.second)/numberOfSnpsOfGroups[g.first];

                // if (read % 1 == 0){
                //     cout << "read " << read << ", defined on " << numberOfSnpsOfGroups[g.first] << " is excluded " << exclusionLevel << " " << g.second << " times from cluster " << g.first 
                //         << " wiht confidence levells " << confidence_level_clusters[g.first].first
                //         << " " << confidence_level_clusters[g.first].second << " " << mean_confidence_level_clusters[g.first].second << endl;
                // }
                if (exclusionLevel <= confidence_level_clusters[g.first].first && exclusionLevel <= confidence_level_clusters[g.first].second){
                    in_this_cluster[g.first] = true;
                }
                if (exclusionLevel < minExclusion ){
                    secondMinExclusion = minExclusion;
                    minExclusion = exclusionLevel;
                    clusterLessExcluded = g.first;
                }
                else if (exclusionLevel < secondMinExclusion){
                    secondMinExclusion = exclusionLevel;
                }
            }
            //if there is exactly one cluster that is not excluded, then this read is in this cluster
            int n = 0;
            int cluster = -1;
            bool foundTwiceOrMore = false;
            for (auto g : in_this_cluster){
                if (g.second && cluster == -1){
                    cluster = g.first;
                }
                else if (g.second && cluster != -1){
                    foundTwiceOrMore = true;
                }
                n++;
            }

            if (!foundTwiceOrMore && cluster != -1){
                haplotypes[read] = cluster;
            }
            else if (minExclusion*1.5 < secondMinExclusion){ //means that it is less exluded from one cluster than the mean
                haplotypes[read] = clusterLessExcluded;
            }
            else {
                // haplotypes[read] = clusteredReads[read];
            }
            // cout << "jlmjpoic " << haplotypes[read]<< endl;


        }
    }

    // cout << "with_which_clusters_goes_this_read : " << endl;
    // int n = 0;
    // for (auto r : with_which_clusters_goes_this_read){
    //     if (clusteredReads[n] != -2){
    //         cout << n << " : ";
    //         for (auto g : r ){
    //             cout << g.first << " " << double(g.second)/on_how_many_snps_is_this_read_defined[n] << ", ";
    //         }
    //         cout << endl;
    //     }
    //     n++;
    // }


    return haplotypes;
}

//input : all the already threaded haplotypes and a new partition
//output : a bool telling if the new partition is compatible with already threaded haplotypes. If yes, alreadyThreadedHaplotype modified
bool extend_with_partition_if_compatible(vector<int> &alreadyThreadedHaplotypes, Partition &extension, int partitionIndex,
    unordered_map <int, pair<int,int>>& clusterLimits){

    //compatibility is defined as : either the 0s or the 1s of the extension all fall squarely within one already defined haplotype
    std::unordered_map <int, int> repartitionOf0s;
    std::unordered_map <int, int> repartitionOf1s;

    bool compatible = false;

    auto idxs = extension.getReads();
    auto content = extension.getPartition();

    int n = 0;
    int numberOf1s = 0;
    int numberOf0s = 0;
    vector<int> extension_vector (alreadyThreadedHaplotypes.size(), 0);
    for (auto idx : idxs){
        if (alreadyThreadedHaplotypes[idx] != -1)
        {
            if (content[n]==1){
                if (repartitionOf1s.find(alreadyThreadedHaplotypes[idx]) == repartitionOf1s.end()){
                    repartitionOf1s[alreadyThreadedHaplotypes[idx]] = 1;
                }
                else {
                    repartitionOf1s[alreadyThreadedHaplotypes[idx]] += 1;
                }
                numberOf1s++;
            }
            else if (content[n]==-1){
                if (repartitionOf0s.find(alreadyThreadedHaplotypes[idx]) == repartitionOf0s.end()){
                    repartitionOf0s[alreadyThreadedHaplotypes[idx]] = 1;
                }
                else {
                    repartitionOf0s[alreadyThreadedHaplotypes[idx]] += 1;
                }
                numberOf0s++;
            }
        }
        extension_vector[idx] = content[n];
        n++;
    }

    //find the best haplotype for 1s
    float max1 = 0;
    int maxClust1 = -1;
    for (auto pair : repartitionOf1s){
        if (pair.first != -1){
            if (pair.second > max1){
                maxClust1 = pair.first;
                max1 = pair.second;
            }
        }
    }

    //find the best haplotype for 0s
    float max0 = 0;
    int maxClust0 = -1;
    for (auto pair : repartitionOf0s){
        if (pair.first != -1){
            if (pair.second > max0){
                maxClust0 = pair.first;
                max0 = pair.second;
            }
        }
    }

    //first see if one haplotype falls squarely in a non-haplotyped zone
    if (max1 == 0 || max0 == 0){
        compatible = true;
        for (int r = 0 ; r < alreadyThreadedHaplotypes.size() ; r++){
            if (extension_vector[r] == 1){
                alreadyThreadedHaplotypes[r] = partitionIndex*2+1;
            }
            else if (extension_vector[r] == -1){
                alreadyThreadedHaplotypes[r] = partitionIndex*2;
            } 
        }
        clusterLimits[partitionIndex*2] = make_pair(extension.get_left(), extension.get_right());
        clusterLimits[partitionIndex*2+1] = make_pair(extension.get_left(), extension.get_right());
    }

    //see if all the 1s fall squarely within one already threaded haplotype
    if (numberOf1s > 0 && max1/numberOf1s > 0.9){ //yes !
        compatible = true;
        if (repartitionOf0s.find(maxClust1) == repartitionOf0s.end() || repartitionOf0s[maxClust1] < 0.1*max1){ //the 0s and the one donn't share maxClust
            for (int r = 0 ; r < alreadyThreadedHaplotypes.size() ; r++){
                if (alreadyThreadedHaplotypes[r] == maxClust1 && extension_vector[r] == -1){ //whuu, was it really well clustered ?
                    alreadyThreadedHaplotypes[r] = -1;
                }
                else if (alreadyThreadedHaplotypes[r] == -1 && extension_vector[r] == 1){ // let's extend maxClust
                    alreadyThreadedHaplotypes[r] = maxClust1;
                }
            }

            //extend the limits of maxClust1
            if (extension.get_left() < clusterLimits[maxClust1].first){
                clusterLimits[maxClust1].first = extension.get_left();
            }
            if (extension.get_right() > clusterLimits[maxClust1].second){
                clusterLimits[maxClust1].second = extension.get_right();
            }
        }
        else{ //then there are 0s on the same already threaded cluster as where the 1s are now : the already threaded cluster was too big
            for (int r = 0 ; r < alreadyThreadedHaplotypes.size() ; r++){
                if (alreadyThreadedHaplotypes[r] == maxClust1 && extension_vector[r] != 1){
                    alreadyThreadedHaplotypes[r] = -1;
                }
                else if (extension_vector[r] == 1){  
                    alreadyThreadedHaplotypes[r] = partitionIndex*2+1;
                }
            }
            clusterLimits.erase(maxClust1);
            clusterLimits[partitionIndex*2+1] = make_pair(extension.get_left(), extension.get_right());
        }
    }

    //see if all the 0s fall squarely within one already threaded haplotype
    if (numberOf0s > 0 && max0/numberOf0s > 0.9 && clusterLimits.find(maxClust0) != clusterLimits.end()){ //yes !
        compatible = true;
        if (repartitionOf1s.find(maxClust0) == repartitionOf1s.end() || repartitionOf1s[maxClust0] < 0.1*max0){ //the 0s and the one donn't share maxClust
            for (int r = 0 ; r < alreadyThreadedHaplotypes.size() ; r++){
                if (alreadyThreadedHaplotypes[r] == maxClust0 && extension_vector[r] == 1){ //whuu, was it really well clustered ?
                    alreadyThreadedHaplotypes[r] = -1;
                }
                else if (alreadyThreadedHaplotypes[r] == -1 && extension_vector[r] == -1){ // let's extend maxClust
                    alreadyThreadedHaplotypes[r] = maxClust0;
                }
            }

            //extend the limits of maxClust0
            if (extension.get_left() < clusterLimits[maxClust0].first){
                clusterLimits[maxClust0].first = extension.get_left();
            }
            if (extension.get_right() > clusterLimits[maxClust0].second){
                clusterLimits[maxClust0].second = extension.get_right();
            }
        }
        else{ //then there are 1s on the same already threaded cluster as where the 0s are now : the already threaded cluster was too big
            for (int r = 0 ; r < alreadyThreadedHaplotypes.size() ; r++){
                if (alreadyThreadedHaplotypes[r] == maxClust0 && extension_vector[r] != -1){
                    alreadyThreadedHaplotypes[r] = -1;
                }
                else if (extension_vector[r] == -1){  //whuu, was it really well clustered ?
                    alreadyThreadedHaplotypes[r] = partitionIndex*2;
                }
            }
            clusterLimits.erase(maxClust0);
            clusterLimits[partitionIndex*2] = make_pair(extension.get_left(), extension.get_right());
        }
    }

    // cout << "compatible ? : " << compatible << " ; " << numberOf0s << " " << numberOf1s << " " << max0 << " " << max1 << endl; 

    return compatible;

}

/**
 * @brief Create a graph linking very similar reads
 * 
 * @param mask 
 * @param snps 
 * @param chunk 
 * @param suspectPostitions 
 * @param sizeOfWindow 
 * @param adjacency_matrix Result
 */
void create_read_graph(
    vector <bool> &mask,
    std::vector<Column> &snps, 
    int chunk, 
    std::vector<size_t> &suspectPostitions,
    int sizeOfWindow,
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs,
    std::vector< std::vector<int>> &adjacency_matrix){

    set<int> listOfGroups;
    int max_cluster = 0;
    unordered_map<int, int> indexOfGroups;

    for (int read1 = 0 ; read1 < mask.size() ; read1 ++){
        if (mask[read1]){
            vector <float> distance_with_other_reads (mask.size(), 0);
            int max_compat = 5; //to remove reads that match on few positions, see how much compatibility you can find at most
            for (auto r = 0 ; r < distance_with_other_reads.size() ; r++){
                if (mask[r] && r != read1 && sims_and_diffs[read1][r].first > 0){
                    distance_with_other_reads[r] = 1 - float(sims_and_diffs[read1][r].second) / float(sims_and_diffs[read1][r].first+sims_and_diffs[read1][r].second);
                    if (sims_and_diffs[read1][r].first > max_compat){
                        max_compat = sims_and_diffs[read1][r].first;
                    }
                }
            }

            // if (read1 == 16){
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
            float distance_threshold = 1;
            for (auto neighbor : smallest){
                if (neighbor.second > 0.5 && nb_of_neighbors < 5 && mask[neighbor.first]){
                    nb_of_neighbors++;
                    distance_threshold = neighbor.second;
                    adjacency_matrix[read1][neighbor.first] = 1;
                    adjacency_matrix[neighbor.first][read1] = 1;
                }
                if (sims_and_diffs[read1][neighbor.first].second > 5 && (1-neighbor.second) > 2*(1-distance_threshold)){
                    // already_separated[read1][neighbor.first] = true;
                    // already_separated[neighbor.first][read1] = true;
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
 * @brief Compute all the diverging and converging snp between the reads on suspect positions
 * 
 * @param mask 
 * @param snps 
 * @param suspectPostitions 
 * @param sims_and_diffs 
 */
void list_similarities_and_differences_between_reads(
    vector <bool> &mask,
    std::vector<Column> &snps, 
    std::vector<size_t> &suspectPostitions,
    vector<vector<pair<int,int>>> &sims_and_diffs){

    set<int> debug_interesting_reads = {-2, -1};
    for (auto r = 0 ; r < mask.size() ; r++){
        if (mask[r]){
            debug_interesting_reads.insert(r);
        }
    }
    std::unordered_map<int, string> debug_strings;
    for (auto r : debug_interesting_reads){
        debug_strings[r] = "";
    }

    for (auto position : suspectPostitions){

        // if (position <200000){
        //     cout<< "ldjflqmj" << endl;
        //     print_snp(snps[position], mask);
        // }

        //what bases occur in which cluster ?
        unordered_map<unsigned char,int> bases_in_total;
        for (auto r = 0 ; r < snps[position].readIdxs.size() ; r++){
            unsigned char base = snps[position].content[r];
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

        debug_strings[-2] += second_most_frequent_base;
        debug_strings[-1] += most_frequent_base;
        

        int idx1 = 0;
        for (int read1 = 0 ; read1 < mask.size() ; read1 ++){
            if (mask[read1]){
                while (idx1 < snps[position].readIdxs.size() && snps[position].readIdxs[idx1] < read1){
                    idx1++;
                }
                if (snps[position].readIdxs[idx1] > read1){
                    if (debug_interesting_reads.find(read1) != debug_interesting_reads.end()){
                        debug_strings[read1] += " ";
                    }
                    continue;
                }
                else{
                    if (debug_interesting_reads.find(read1) != debug_interesting_reads.end()){
                        if (bases_in_total[snps[position].content[idx1]] >= second_most_frequent_base_count){
                            debug_strings[read1] += snps[position].content[idx1];
                        }
                        else{
                            debug_strings[read1] += " ";
                        }
                    }
                }
                unsigned char base1 = snps[position].content[idx1];

                int idx2 = idx1+1;
                for (int read2 = read1+1 ; read2 < mask.size() ; read2 ++){
                    if (mask[read2]){
                        while (idx2 < snps[position].readIdxs.size() && snps[position].readIdxs[idx2] < read2){
                            idx2++;
                        }
                        if (snps[position].readIdxs[idx2] > read2){
                            continue;
                        }
                        unsigned char base2 = snps[position].content[idx2];

                        // if (bases_in_total[base1] > second_most_frequent_base_count && bases_in_total[base2] >= second_most_frequent_base_count &&
                        //     read1 == 2 && read2 == 83){
                        //     cout << "posodddiicie ";
                        //     print_snp(snps[position], mask);
                        //     cout << "baddese1 " << base1 << " " << bases_in_total[base1] << " base2 " << base2 << " " << bases_in_total[base2] << endl;
                        // }

                        // cout << "qijdioqddsp " << bases_in_total[base1] << " " << second_most_frequent_base_count << " " << bases_in_total[base2] << " " << base1  << " "<< base2 << endl;

                        if (bases_in_total[base1] == second_most_frequent_base_count && bases_in_total[base2] >= second_most_frequent_base_count && base1 != base2){
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
    }

    // cout << "eieieii here are the reads on the interesting positions" << endl;
    // for (auto r : debug_interesting_reads){
    //     cout << debug_strings[r] << " " << r << endl;
    // }
    // cout << "sims and diffs between reads 41, 43, 133" << endl;
    // cout << "41 43 " << sims_and_diffs[41][43].first << " " << sims_and_diffs[41][43].second << endl;
    // cout << "41 133 " << sims_and_diffs[41][133].first << " " << sims_and_diffs[41][133].second << endl;
    // cout << "43 133 " << sims_and_diffs[43][133].first << " " << sims_and_diffs[43][133].second << endl;
}

std::vector<int> merge_clusterings(std::vector<std::vector<int>> &localClusters,
    std::vector< std::vector<int>> &adjacency_matrix){

    if (localClusters.size() == 0){
        return vector<int>(adjacency_matrix.size(), 0);
    }


    vector<double> clusters_aggregated(localClusters[0].size(), 0);
    for (auto i = 0 ; i < localClusters.size() ; i++){
        for (auto j = 0 ; j < localClusters[i].size() ; j++){
            clusters_aggregated[j] += localClusters[i][j]*(std::hash<int>()(i) % 1000);
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

    // cout << "cluster_aggregated_ints : " << endl;
    // for (auto i = 0 ; i < cluster_aggregated_ints.size() ; i++){
    //     cout << cluster_aggregated_ints[i] << " ";
    // }
    // cout << endl;

    auto re_clustered = chinese_whispers(adjacency_matrix, cluster_aggregated_ints);

    // cout << "re_clustered : " << endl;
    // for (auto i = 0 ; i < re_clustered.size() ; i++){
    //     cout << re_clustered[i] << " ";
    // }
    // cout << endl;

    return re_clustered;
}



