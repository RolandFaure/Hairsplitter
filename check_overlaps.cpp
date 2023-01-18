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
            
            if (allreads[read].neighbors_.size() > 10 && allreads[read].name != "consensus@0@00"){

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
    float meanDistance = generate_msa(contig, allOverlaps, allreads, snps, insertionPositions, partitions.size(), truePar, assemble_on_assembly, readLimits, misalignedReads, polish);

    //then separate the MSA
    string ref = allreads[contig].sequence_.str();
    vector<pair<pair<int,int>, vector<int>> > par = separate_reads(ref, snps, meanDistance,
                                                        allreads[contig].neighbors_.size()+1-int(assemble_on_assembly));
    
        
    // auto par = truePar.getPartition();//DEBUG
    // for (auto i = 0 ; i < par.size() ; i++) {par[i]++; }
    // cout << "Proposed partition : " << endl;
    // for (auto i = 0 ; i < par.size() ; i++){cout << par[i];}cout << endl;
    // cout << endl;

    //now compute the consensus for each new contig
    compute_consensus_in_partitions(contig, par, allreads, allOverlaps, snps, insertionPositions, partitions);

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
    unordered_map <int, vector<pair<int,int>>> &readLimits, std::vector<bool> &misalignedReads, bool polish){

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
    for (auto n = 0 ; n< allreads[bbcontig].neighbors_.size() ; n++){
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


    //create a vector of alternative reference sequences (chosen from the polishing reads)
    vector<string> allTheReferences = {allreads[bbcontig].sequence_.str()};
    wfa::WFAlignerGapAffine aligner(4,6,2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);

    for (auto n = 0 ; n < polishingReads.size() ; n++){

        if (DEBUG && n%10 == 0){
            cout << "Aligned " << n << " reads out of " << allreads[bbcontig].neighbors_.size() << " on the backbone\r" << std::flush;
        }

        auto t1 = high_resolution_clock::now();

        string alignment;
        if (CIGARs.size() != polishingReads.size()){ //compute the exact alignment if it is not already computed
            //print the name of the neighbor
            string cons = consensus.substr(positionOfReads[n].first, positionOfReads[n].second-positionOfReads[n].first).c_str();
            string re = polishingReads[n].c_str();

            // wfa::WFAlignerGapAffine aligner(4,6,2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
            // int marginLeftEnd = std::min(int(cons.size()/4), 1000);
            // int marginRightEnd = std::min(int(re.size()/4), 1000);
            // aligner.alignEndsFree(cons,marginLeftEnd,marginLeftEnd, re, marginRightEnd, marginRightEnd); //the ints are the length of gaps that you tolerate for free at the ends of query and sequence
            aligner.alignEnd2End(cons, re);

            alignment = aligner.getAlignmentCigar();

            // if (allreads[allOverlaps[allreads[bbcontig].neighbors_[n]].sequence1].name.substr(0,10) == "@3_ade377e"){
            //     cout << "yuiiee Aligned using substr: " << positionOfReads[n].first << " " << positionOfReads[n].second << endl; 
            //     cout << "fljmqqzzzo Alignment: " << alignment << endl;
            // }
        }
        else{
            alignment = convert_cigar(CIGARs[n]);
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

        int localDivergence = 0; //number of mismatches/insertions/deletions in the last 100 bases
        vector<bool> divergences (100, false); //true if there is a mismatch/insertion/deletion in that position

        // if (indexQuery < 20){
        //     cout << "beginning of query : " << indexQuery << " " << polishingReads[n].substr(indexQuery,100) << endl;
        // }
        //cout the alignment
        // cout << "alignment : " << alignment << endl;

        for (int l = 0; l < alignment.size(); l++) {

            if (indexQuery < consensus.size()){

                // if (l >= 100){
                //     if (divergences[l%100]){
                //         localDivergence -= 1;
                //         divergences[l%100] = false;
                //     }
                // }

                if (alignment[l] == '=' || alignment[l] == 'X' || alignment[l] == 'M'){
                    //fill inserted columns with '-' just before that position
                    // numberOfConsecutiveMatches += 1;
                    for (int ins = numberOfInsertionsThere ; ins < numberOfInsertionsHere[indexQuery] ; ins++){ //in these positions, insert '-' instead of '?'
                        snps[insertionPos[10000*indexQuery+ins]].readIdxs.push_back(n);
                        snps[insertionPos[10000*indexQuery+ins]].content.push_back('-');
                    }

                    if (indexQuery > 0 && polishingReads[n][indexTarget-1] != consensus[indexQuery-1]){
                        totalDistance += 1;
                    }

                    snps[indexQuery].readIdxs.push_back(n);
                    snps[indexQuery].content.push_back(polishingReads[n][indexTarget]);

                    indexQuery++;
                    indexTarget++;
                    numberOfInsertionsThere = 0;

                    if (alignment[l] != 'M' && alignment[l] != '='){
                        divergences[l%100] = true;
                        localDivergence += 1;
                    }
                }
                else if (alignment[l] == 'S'){ //soft-clipped bases
                    indexTarget++;
                    divergences[l%100] = true;
                    localDivergence += 1;
                }
                else if (alignment[l] == 'D'){
                    //fill inserted columns with '-' just before that position
                    for (int ins = numberOfInsertionsThere ; ins < numberOfInsertionsHere[indexQuery] ; ins++){ //in these positions, insert '-' instead of '?'
                        snps[insertionPos[10000*indexQuery+ins]].readIdxs.push_back(n);
                        snps[insertionPos[10000*indexQuery+ins]].content.push_back('-');
                    }

                    snps[indexQuery].readIdxs.push_back(n);
                    snps[indexQuery].content.push_back('-');                    
                    indexQuery++;
                    numberOfInsertionsThere = 0;

                    totalDistance += 1;

                    divergences[l%100] = true;
                    localDivergence += 1;
                }
                else if (alignment[l] == 'I'){ //hardest one
                    if (numberOfInsertionsHere[indexQuery] <= maxNumberOfInsertions && indexQuery > positionOfReads[n].first) {

                        if (numberOfInsertionsThere >= numberOfInsertionsHere[indexQuery]) { //i.e. this is a new column
                            insertionPos[10000*indexQuery+numberOfInsertionsHere[indexQuery]] = snps.size();
                            numberOfInsertionsHere[indexQuery] += 1;
                            
                            Column newInsertedPos;
                            newInsertedPos.readIdxs = snps[indexQuery-1].readIdxs;    
                            newInsertedPos.content = vector<char>(snps[indexQuery-1].content.size() , '-');
                            newInsertedPos.content[newInsertedPos.content.size()-1] = polishingReads[n][indexTarget];
                            newInsertedPos.pos = l;
                            snps.push_back(newInsertedPos);
                        }
                        else{
                            snps[insertionPos[10000*indexQuery+numberOfInsertionsThere]].readIdxs.push_back(n);
                            snps[insertionPos[10000*indexQuery+numberOfInsertionsThere]].content.push_back(polishingReads[n][indexTarget]);
                        }
                        numberOfInsertionsThere ++;
                    }
                    indexTarget++;
                    totalDistance += 1;
                    
                    divergences[l%100] = true;
                    localDivergence += 1;
                }

                if (localDivergence > 30){
                    misalignedReads[n] = true;
                    // if (allreads[allOverlaps[allreads[read].neighbors_[n]].sequence1].name[1] == '3'){
                    //     cout << "Read " << allreads[allOverlaps[allreads[read].neighbors_[n]].sequence1].name << " is misaligned: " << alignment <<" \n";
                    // }
                }

            }
            //cout << l << " " << result.alignmentLength <<  "\n";
        } 

        auto t3 = high_resolution_clock::now();

        alignmentTime += duration_cast<milliseconds>(t2-t1).count();
        MSAtime += duration_cast<milliseconds>(t3-t2).count();
        
    }

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
    
    //print snps (just for debugging)
    /*
    int step = 1; //porportions of reads
    int prop = 1; //proportion of positions
    int firstRead = 0;
    int lastRead = polishingReads.size();
    int numberOfReads = lastRead-firstRead;
    int start = 99900;
    int end = 100000;
    vector<string> reads (int(numberOfReads/step));
    string cons = "";
    for (unsigned int i = start ; i < end; i+=prop){
        for (short n = 0 ; n < numberOfReads ; n+= step){
            char c = '?';
            int ri = 0;
            int soughtRead = firstRead+n;
            for (auto r : snps[i].readIdxs){
                if (r == soughtRead){
                    c = snps[i].content[ri];
                }
                ri ++;
            }
            reads[n/step] += c;
        }
        for (short insert = 0 ; insert < min(9999,numberOfInsertionsHere[i]) ; insert++ ){
            int snpidx = insertionPos[10000*i+insert];
            for (short n = 0 ; n < numberOfReads*step ; n+= step){
                char c = '?';
                int ri = 0;
                for (auto r : snps[snpidx].readIdxs){
                    if (r == n){
                        c = snps[snpidx].content[ri];
                    }
                    ri ++;
                }
                reads[n/step] += c;
            }
        }
    }
    cout << "Here are the aligned reads : " << endl;
    int index = firstRead;
    for (auto neighbor : reads){
        if (neighbor[0] != '?'){
            cout << neighbor << " " << index << " " << allreads[allOverlaps[allreads[bbcontig].neighbors_[index]].sequence1].name.substr(0,10) << endl;
        }
        index+= step;
    }
    int n = 0;
    for(auto i : consensus.substr(start, end-start)){
        if (n%prop == 0){
            cout << i;
        }
        n+=1;
    } cout << endl;
    cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;
    */
    return totalDistance/totalLengthOfAlignment;

    //*/
}

/**
 * @brief Function that uses racon to polish a sequence
 * 
 * @param backbone sequence to be polished
 * @param polishingReads list of reads to polish it
 * @param id an id (typically a thread id) to be sure intermediate files do not get mixed up with other threads 
 * @return polished sequence 
 */
string consensus_reads(string &backbone, vector <string> &polishingReads, string &id){
    
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
    auto before_size = min(size_t(500), backbone.size());
    auto after_size = min(size_t(250), consensus.size());

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

    result = edlibAlign(after_start.c_str(), after_start.size(), before_start.c_str(), before_start.size(),
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
 * @param ref reference sequence against which the reads are aligned
 * @param snps MSA
 * @param meanDistance estimation of the error rate
 * @param numberOfReads number of reads of the MSA
 * @return vector<pair<pair<int,int>, vector<int>> >  pairs of coordinates to which are attached a certain partition of reads
 */
vector<pair<pair<int,int>, vector<int>> > separate_reads(string& ref, std::vector<Column> &snps, float meanDistance, int numberOfReads){

    /*
    The null model is described as uniform error rate -> binomial error distribution on one position
    meanDistance ~= 2*sequencingErrorRate (e) at most
    suspicious positions have p-value < 0.05, i.e. more errors than n*e + 3*sqrt(n*e*(1-e))
    */

    robin_hood::unordered_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['-'] = 4;
    bases2content['?'] = 5;

    vector<Partition> partitions; //list of all partitions of the reads, with the number of times each occurs

    int numberOfSuspectPostion = 0;
    int numberOfNeighbors = 0;

    //two variables to know how much a position diverges from the consensus on average
    float meanError = 0;
    int numberOfExtensions = 0;

    vector<vector<distPart>> distanceBetweenPartitions;
    vector<size_t> suspectPostitions;
    vector <int> interestingParts; //DEBUG
    int localErrors = 0;

    for (int position = 0 ; position < ref.size() ; position++){ 

        if (DEBUG && position%100 == 0){
            cout << "Going through the positions, " << position << "/" << ref.size() << "          \r" << std::flush;
        }
        //first look at the position to see if it is suspect
        int content [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for -
        int numberOfReadsHere = 0;
        for (short n = 0 ; n < snps[position].content.size() ; n++){
                char base = snps[position].content[n];
                if (base != '?' && bases2content.contains(base)){
                    content[bases2content[base]] += 1;
                    numberOfReadsHere += 1;
                }
        }

        //float threshold = min(1 + numberOfReadsHere*meanDistance/2 + 3*sqrt(numberOfReadsHere*meanDistance/2*(1-meanDistance/2)), float(0.1*numberOfReadsHere)); //more than 10% of different bases is always suspicious
        //sort content and store the result in content_sorted
        int content_sorted[5] = {0,0,0,0,0};
        std::partial_sort_copy(content, content+5, content_sorted, content_sorted+5, std::greater<int>());

        float probabilityTrash = float(content_sorted[2]/numberOfReadsHere); //probability of having a trash base at this position

        if (content[4] < content_sorted[1] && content_sorted[1] > content_sorted[2]+3*sqrt(content_sorted[2]*(1-probabilityTrash)) /*3*content_sorted[2]*/ && content_sorted[1] > 4){ //this position is suspect (for now positions with too many gaps are not considered)
           
            // cout << "iouocxccv " << threshold << " " << position << " ;bases : " << content[0] << " " << content[1] << " " << content[2] << " " << content[3] << " " << content[4] << endl;
            char ref_base;
            if (position < ref.size()){
                ref_base = ref[position];
            } else {
                ref_base = '-';
            }

            suspectPostitions.push_back(position);
            //go through the partitions to see if this suspicious position looks like smt we've seen before
            vector<distPart> distances (distanceBetweenPartitions.size());
            bool found = false;
            for (auto p = 0 ; p < partitions.size() ; p++){
                //if the partition is too far away, do not bother comparing
                if (std::abs(snps[position].pos-partitions[p].get_right())>10000){
                    continue;
                }

                distancePartition dis = distance(partitions[p], snps[position], ref_base);
                auto comparable = min(dis.n00,dis.n11) + dis.n01 + dis.n10;
                int trashBases = content_sorted[2] + content_sorted[3] + content_sorted[4]; //bases that are probably errors at this position

                float tolerableNumberOfErrors = float(trashBases)/2;
                if (partitions[p].number() < 10){
                    tolerableNumberOfErrors = trashBases;
                }

                //if ((float(dis.n01+dis.n10)/(min(dis.n00,dis.n11)+dis.n01+dis.n10) <= meanDistance*2 || dis.n01+dis.n10 <= 2)  && dis.augmented && comparable > min(10.0, 0.3*numberOfReads)){
                if(dis.n01 <= localErrors && dis.n10 <= tolerableNumberOfErrors && dis.n00 >= max(5, dis.n01+dis.n10) && dis.augmented && comparable > min(10.0, 0.3*numberOfReads)){  
                    int pos = -1;
                    if (position < ref.size()){
                        pos = position;
                    }

                    partitions[p].augmentPartition(dis.partition_to_augment, pos);
                    found = true;

                    meanError += float(dis.n01+dis.n10)/(dis.n00 + dis.n11 + dis.n01 + dis.n10);
                    numberOfExtensions += 1;

                    break;
                }
                if (comparable > 5){
                    //distances[p].distance = max(float(0), (10-dis.chisquare)/10) ;
                    distances[p].distance = float(dis.n01+dis.n10)/(dis.n00 + dis.n11 + dis.n01 + dis.n10);
                    distances[p].phased = dis.phased;
                }
                else{
                    distances[p].distance = 1;
                }
            }

            if (!found && position < ref.size()){    // the second condition is here to create partitions only at specific spots of the backbone
                distanceBetweenPartitions.push_back(distances);
                partitions.push_back(Partition(snps[position], position, ref_base)); 
            }
            else{
                position += 5;
            }
            
            numberOfSuspectPostion += 1;

            //two suspect positions next to each other can be artificially correlated through alignement artefacts

        }
        else{
            localErrors = content_sorted[2] + content_sorted[3] + content_sorted[4];
        }
    }

    meanError /= numberOfExtensions;

    if (DEBUG){
        cout << "found " << numberOfSuspectPostion << " suspect positions" << endl;
    }

    if (partitions.size() == 0){ //there are no position of interest
        vector<pair<pair<int,int>, vector<int>>> e;
        return e;
    }

    
    float threshold =  max(4.0, min(0.01*numberOfSuspectPostion, 0.001*snps.size()));

    vector<Partition> listOfFinalPartitions;
    for (auto p1 = 0 ; p1 < partitions.size() ; p1++){

        // if (partitions[p1].number() > 5){
        //     cout << "iqdoudofq non informative partition : "  << threshold << " " << numberOfSuspectPostion << " " << snps.size()<< endl;
        //     partitions[p1].print();
        // }
        
        if (partitions[p1].number() > threshold && partitions[p1].isInformative(meanDistance/2, true)){


            bool different = true;
            
            for (auto p2 = 0 ; p2 < listOfFinalPartitions.size() ; p2++){

                distancePartition dis = distance(listOfFinalPartitions[p2], partitions[p1], 2);
                if (dis.augmented){
                    Partition newPart = listOfFinalPartitions[p2];
                    newPart.mergePartition(partitions[p1], dis.phased);

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

    //now we have the list of final partitions : there may be several, especially if there are more than two haplotypes

    //now go through windows of width 1000 along the reference and create local partitions
    vector<pair<pair<int,int>, vector<int>>> threadedReads;
    int suspectPostitionIdx = 0;
    int sizeOfWindow = 1000;
    for (auto chunk = 0 ; chunk < ref.size()/sizeOfWindow ; chunk++){
        vector<Partition> localPartitions(listOfFinalPartitions.size());

        //let's see what partitions are found on the local window
        while(suspectPostitionIdx < suspectPostitions.size() && suspectPostitions[suspectPostitionIdx] < (chunk+1)*sizeOfWindow){

            auto pos = suspectPostitions[suspectPostitionIdx];
            //find the closest partition
            int closestPartition = -1;
            float closestDistance = 1;
            for (auto p = 0 ; p < localPartitions.size() ; p++){
                //compute the distance between the partition and the suspect position
                char ref_base = ref[pos];
                // listOfFinalPartitions[p].print();
                //print snps[pos];
                // for (auto i : snps[pos].content){
                //     cout << i;
                // }
                // cout << endl;
                distancePartition dis = distance(listOfFinalPartitions[p], snps[pos], ref_base);
                float distance = min(float(dis.n01+dis.n10)/(dis.n00 + dis.n11 + dis.n01 + dis.n10), float(dis.n00+dis.n11)/(dis.n00 + dis.n11 + dis.n01 + dis.n10));
                if (distance < closestDistance){
                    closestDistance = distance;
                    closestPartition = p;
                }
            }
            //now we have the closest partition, let's see if we can integrate this position in it
            if (closestDistance <= max(meanError*2, float(0.1)) && closestPartition != -1){
                //we can integrate the position in the partition
                localPartitions[closestPartition].augmentPartition(snps[pos], pos);
            }
            suspectPostitionIdx += 1;
        }

        //go through the local partitions and strengthen them
        for (auto p = 0 ; p < localPartitions.size() ; p++){
            localPartitions[p].strengthen_partition(listOfFinalPartitions[p]);
        }

        //list and strengthen local partitions that have a numberOfOccurences >= 1
        for (auto p = 0 ; p < localPartitions.size() ; p++){
            if (localPartitions[p].number() >= 1){
                cout << "on windhhow " << chunk << " found partition " << endl;
                localPartitions[p].print();
            }
        }

        //create non-null-partitions, compased of local partitions that have a number() >= 1
        vector<Partition> non_null_partitions;
        for (auto p = 0 ; p < localPartitions.size() ; p++){
            if (localPartitions[p].number() >= 1){
                non_null_partitions.push_back(localPartitions[p]);
            }
        }

        //select compatible partitions
        vector<Partition> listOfCompatiblePartitions = select_compatible_partitions(non_null_partitions, numberOfReads, meanError);
        for (auto p = 0 ; p < listOfCompatiblePartitions.size() ; p++){
            cout << "on window " << chunk << " found coimpatdsible partition " << endl;
            listOfCompatiblePartitions[p].print();
        }

        //thread partition on this window
        vector<int> clusteredRead = threadHaplotypes_in_interval(listOfCompatiblePartitions, numberOfReads);

        threadedReads.push_back(make_pair(make_pair(chunk*sizeOfWindow, (chunk+1)*sizeOfWindow-1), clusteredRead));


        // if (localPartitions.size() != 0 && chunk == 20){
        //     cout << "on last window, clustered reads:" << endl;
        //     for (auto i : clusteredRead){
        //         if (i > -1){
        //             cout << i;
        //         }
        //         else{
        //             cout << " ";
        //         }
        //     }
        //     cout << endl;

        //     for (auto p : listOfCompatiblePartitions){
        //         if (p.number() >= 1){
        //             cout << "pausing" << endl;
        //             while(true){}
        //         }
        //     }
        // }

    }
    return threadedReads;
    
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

    vector <unsigned int> idxs2 = par2.readIdxs;
    vector <char> part2 = par2.content;

    robin_hood::unordered_flat_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['-'] = 4;
    bases2content['?'] = 5;
    
    int content2 [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for *, 5 for '?'
    float numberOfBases = 0;

    auto n2 = 0;
    auto n1 = 0;
    for (auto r : idxs2){
        while (n1 < idxs1.size() && idxs1[n1] < idxs2[n2]){
            n1++;
        }
        if (n1 > idxs1.size()){
            break;
        }
        if ( idxs1[n1] == idxs2[n2]){
            numberOfBases+=1;
            content2[bases2content[part2[n2]]] += 1;
        }
        n2++;
    }

    if (numberOfBases == 0){ //not comparable
        res.n00 = 0;
        res.n01 = 0;
        res.n10 = 0;
        res.n11 = 0;
        res.augmented = false;
        return res;
    }

    //determine first and second most frequent bases in par2
    
    char mostFrequent = ref_base;
    auto maxFrequence = content2[bases2content[ref_base]];
    //now find the second most frequent base
    char secondFrequent = 'O';
    int maxFrequence2 = -1;
    for (auto i = 0 ; i < 5 ; i++){
        if (ref_base != "ACGT-"[i]) {
            if (content2[i] > maxFrequence2){
                secondFrequent = "ACGT-"[i];
                maxFrequence2 = content2[i];
            }
        }
    }

    // cout << "poiuy Two most frequent : " << mostFrequent2 << "," << secondFrequent2 << "\n";

    float score = 0; //the scores when directing mostFrequent on either mostfrequent2 or secondFrequent2
    float bestScore = 0;
    //remember all types of matches for the chi square test
    int matches00 = 0;
    int matches01 = 0;
    int matches10 = 0;
    int matches11 = 0;
    Column newPartition;

    newPartition.readIdxs = {};
    for (auto ci = 0 ; ci < par2.content.size() ; ci++){
        char c = par2.content[ci];
        if (c == mostFrequent){
            newPartition.readIdxs.push_back(par2.readIdxs[ci]);
            newPartition.content.push_back('A');

        }
        else if (c == secondFrequent){
            newPartition.readIdxs.push_back(par2.readIdxs[ci]);
            newPartition.content.push_back('a');
        }
        else{
            //the other bases are not considered
        }
    }


    n1 = 0;
    auto it1 = idxs1.begin();
    n2 = 0;
    auto it2 = par2.readIdxs.begin();
    while (it1 != idxs1.end() && it2 != par2.readIdxs.end()){

        float conf = 2*confs1[n2]-1;

        if (*it1 == *it2){
            ++it1;
            ++it2;

            if (par2.content[n2] == mostFrequent){

                if (part1[n1] == 1){
                    score += conf;
                    bestScore += conf;
                    matches11 += 1;
                }
                else if (part1[n1] == -1){
                    score -= conf;
                    bestScore += conf;
                    matches01 += 1;
                }
            }
            else if (par2.content[n2] == secondFrequent){

                if (part1[n1] == 1){
                    score -= conf;
                    bestScore += conf;
                    matches10 += 1;
                }
                else if (part1[n1] == -1){
                    score += conf;
                    bestScore += conf;
                    matches00 += 1;
                }
            }
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
    res.phased = 1;
    res.secondBase = secondFrequent;
    res.partition_to_augment = newPartition;
    //cout << "Computing..." << maxScore << " " << par1.size()-res.nonComparable << " " << res.nmismatch << endl;

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
        else{

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
                }
            }

            r1++;
            r2++;
        }
    }


    distancePartition res;
    res.augmented = true;

    //check if there are too many unenxplainable positions
    if ((ndivergentPositions[0] >= threshold_p && ndivergentPositions[1] >= threshold_p) || numberOfComparableBases == 0){
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
 * @brief Given the original partition, forces each read to choose a camp given all suspect positions : some reads will change camp...
 * 
 * @param originalPartition Partition that we want to clean
 * @param snps 
 * @param suspectPositions Reads will choose based on all suspect positions
 * @return Partition 
 */
void clean_partition(long int backbone, Partition &originalPartition, vector <Read> &allreads, vector <Overlap> &allOverlaps){

    //first compute all the contigs as they would look like on one part or another
    //inventoriate all the reads
    if (DEBUG){
        cout << "cleaning partition..." << endl;
    }
    vector <string> polish_reads_0;
    vector <string> polish_reads_1;

    auto part = originalPartition.getPartition();
        
    std::ofstream polishseqs("tmp/polish.fasta", std::ofstream::trunc);

    int n = 0;
    for (auto ov : originalPartition.getReads()){
        Overlap overlap = allOverlaps[allreads[backbone].neighbors_[ov]];
        string read = allreads[overlap.sequence1].sequence_.str();
        
        if (part[n] == -1){
            string r = read.substr(overlap.position_1_1, overlap.position_1_2-overlap.position_1_1);
            polish_reads_0.push_back(r);
            polishseqs << ">"+std::to_string(n)+"\n" << r << "\n";
        }
        else if (part[n] == 1){
            string r = read.substr(overlap.position_1_1, overlap.position_1_2-overlap.position_1_1);
            polish_reads_1.push_back(r);
            polishseqs << ">"+std::to_string(n)+"\n" << r << "\n";
        }
        n++;
    }
    polishseqs.close();

    //use the reads to polish the backbone
    string toPolish = allreads[backbone].sequence_.str();
    string thread_id = std::to_string(omp_get_thread_num());
    string consensus0 = consensus_reads(toPolish, polish_reads_0, thread_id);
    string consensus1 = consensus_reads(toPolish, polish_reads_1, thread_id);

    //now compare all reads to the two new polished consensus and re-assign them to their best match

    std::ofstream consseqs("tmp/consensusSeqs.fasta");
    consseqs << ">0\n" << consensus0 << "\n>1\n" << consensus1 << "\n";
    consseqs.close();

    string com = " -x map-ont -N 0 -t1 tmp/consensusSeqs.fasta tmp/polish.fasta > tmp/choose_consensus.paf 2> tmp/trash.txt";
    string command = MINIMAP+com;
    system(command.c_str());

    //now read the output
    auto correctedPartition = originalPartition.getPartition();
    std::ifstream in("tmp/choose_consensus.paf");
    string line;
    while(getline(in, line)){
        string field;
        std::istringstream line2(line);
        int fieldnb = 0;
        int readnb = -1;
        while (getline(line2, field, '\t'))
        {
            if (fieldnb == 0){ //read
                readnb = atoi(field.c_str());
            }
            else if (fieldnb == 5){ //the part on which the read mapped
                if (field == "0"){
                    correctedPartition[readnb] = -1;
                }
                else if (field == "1"){
                    correctedPartition[readnb] = 1;
                }
            }
            fieldnb++;
        }
    }
    in.close();

    if (DEBUG){
        cout << "correcting this partition : " << endl;
        originalPartition.print();
        originalPartition.new_corrected_partition(correctedPartition, originalPartition.getReads(), originalPartition.getMore(), originalPartition.getLess());
        originalPartition.print();

        cout << "Done cleaning ! " << endl;
    }
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

    // if (DEBUG){cout << "compatible partitions : " << endl;}
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
        p1.flipPartition();
        compatible += 1;
    }
    else if(zero_in_1/numberOf0s > 0.8 && zero_in_0/(numberOf0s+numberOf1s) < 0.1){
        p1.flipPartition();
        p2.flipPartition();
        compatible += 1;
    }

    if (one_in_0/numberOf1s > 0.8 && one_in_1/(numberOf0s+numberOf1s) < 0.1){
        compatible += 1;
    }
    else if (one_in_1/numberOf1s > 0.8 && one_in_0/(numberOf0s+numberOf1s) < 0.1){
        p2.flipPartition();
        compatible += 1;
    }

    if(DEBUG){cout << "Compatibility : " << numberOf0s << " " << zero_in_0 << " " << zero_in_1 << " ; " << numberOf1s << " " << one_in_0 << " " << one_in_1 << " ; " << compatible << endl;}
    return compatible;
}

/**
 * @brief Concatenates binary partitions to form a single multi-part partition
 * 
 * @param listOfFinalPartitions the binary partitions
 * @param numberOfReads 
 * @return vector<int> the final partition
 */
vector<int> threadHaplotypes_in_interval(vector<Partition> &listOfFinalPartitions, int numberOfReads){

    if (listOfFinalPartitions.size() == 0){
        return vector<int>(numberOfReads, 0);
    }

    vector<vector<short>> allPartitions;
    for (auto i : listOfFinalPartitions){
        allPartitions.push_back(i.getPartition());
    }
    vector<vector<float>> allConfidences;
    for (auto i : listOfFinalPartitions){
        allConfidences.push_back(i.getConfidence());
    }

    int numberOfAssignedReads=0;
    //inventroriate the reads that are present on all partitions
    vector<int> numberOfTimeOfEachRead (numberOfReads, 0);
    for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
        for (auto read : listOfFinalPartitions[binary].getReads()){
            numberOfTimeOfEachRead[read] += 1;
        }
    }
    //now we have the number of times each read is present in the partitions, create a vector of reads that are present in all partitions
    vector<int> readsPresentInAllPartitions;
    for (auto i = 0 ; i < numberOfReads ; i++){
        if (numberOfTimeOfEachRead[i] == listOfFinalPartitions.size()){
            readsPresentInAllPartitions.push_back(i);
        }
    }

    vector<int> clusters (numberOfReads, -1); //this will contain only high-confidence reads
    vector<int> clustersAll (numberOfReads, -1); //this will contain all reads
    vector<int> presentReads (numberOfReads, listOfFinalPartitions.size()); //this will count to 0 the number of partitions in which the read is present

    for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
        int c = 0;

        //for (auto read : readsPresentInAllPartitions){
        for (auto read : listOfFinalPartitions[binary].getReads()){
            if (numberOfTimeOfEachRead[read] == listOfFinalPartitions.size()){

                presentReads[read] -= 1;

                auto camp = allPartitions[binary][c];
                if (clusters[read] == -1 && camp != 0){ //means that the read is present in one partition
                    clusters[read] = 0;
                    clustersAll[read] = 0;
                }

                if (allConfidences[binary][c] < 0.7 || camp == 0){  //0.7 to be pretty confident about the reads we separate (more than 70% on all partitions)
                    clusters[read] = -pow(2, listOfFinalPartitions.size())-2; //definitely putting that below -1
                }
                else if (camp != 0){
                    clusters[read] += pow(2, binary)*int(0.5+0.5*camp); //*0 if haplotype -1, *1 if haplotype 1
                }
                if (camp!=0){
                    clustersAll[read] += pow(2, binary)*int(0.5+0.5*camp); //*0 if haplotype -1, *1 if haplotype 1
                }
            }
            c += 1;
        }

    }

    set<int> count; //a set of all existing groups
    vector<int> frequenceOfPart (pow(2, listOfFinalPartitions.size())); 
    for (auto id : clusters){
        if (id >= 0){
            frequenceOfPart[id] += 1;
            count.emplace(id);
            numberOfAssignedReads ++;
        }
    }
    
    /*
        As a first approximation:
        to make sure we don't over-estimate the number of clusters, we'll make the assumption that all haplotypes have haplotype-specific mutation => the number of final cluster cannot be higher than the number of partitions
        Then:
        also count the number of solid clusters (i.e. frequent)
    */
    vector<float> proportionOf1;
    for (auto p : listOfFinalPartitions){
        proportionOf1.push_back(p.proportionOf1());  
    }
    vector<float> listOfLikelihoods;
    vector<float> listOfProbas;
    size_t max_number_of_cluster = 0; //count the number of clusters that look solid
    for (auto group : count){
        int id = group;
        double proba = 1;
        for (int binary = listOfFinalPartitions.size()-1 ; binary > -1 ; binary--){
            if (id%2 == 1){
                proba *= proportionOf1[binary];
            }
            else if (id%2==0){
                proba *= 1-proportionOf1[binary];
            }
            id = int(id / 2);
        }
        // if (DEBUG){
        //     cout << "group " << group << " is expected " << numberOfAssignedReads*proba << " times and arrives " << frequenceOfPart[group] << " times " << endl;
        // }
        listOfLikelihoods.push_back(float(frequenceOfPart[group]-numberOfAssignedReads*proba) / frequenceOfPart[group]);
        listOfProbas.push_back(proba);
        if (frequenceOfPart[group] > 10 && numberOfAssignedReads*proba < frequenceOfPart[group]*3){
            max_number_of_cluster++;
        }
        // cout << "group " << group << " likelihood: " << float(frequenceOfPart[group] -numberOfAssignedReads*proba) / frequenceOfPart[group] << endl;
    }

    float minLikelihood = -1000;
    if (listOfLikelihoods.size() > listOfFinalPartitions.size() && listOfFinalPartitions.size()>1){
        int extra = 0; //when there is only 1 partitions there is actually 2 parts
        if (listOfFinalPartitions.size()<=2){extra=1;}
        vector <float> minElements (listOfLikelihoods.size());
        std::partial_sort_copy(listOfLikelihoods.begin(),  listOfLikelihoods.end(), minElements.begin(), minElements.end());
        if (minElements.size() >= listOfFinalPartitions.size()+extra){
            auto number_of_cluster = max(max_number_of_cluster, listOfFinalPartitions.size()+extra);
            minLikelihood = minElements[minElements.size()-number_of_cluster];
        }
        else{
            minLikelihood = minElements[0];
        }
    }

    
    //establish a list of possible clusters
    std::set <int> listOfGroups;
    int indexOfCount = 0;

    int numberOfAssignedReads2 = numberOfAssignedReads;
    float newtotalproba = 1;
    for (auto group : count){
        if (listOfLikelihoods[indexOfCount] >= minLikelihood && frequenceOfPart[group] > 5){ //>5 because we want at least 5 reads per haplotype
            listOfGroups.emplace(group);
            numberOfAssignedReads2 -= frequenceOfPart[group];
            newtotalproba -= listOfProbas[indexOfCount];
        }
        indexOfCount++;
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

    // cout << "possible clusters : "; for (auto g : listOfGroups){cout << g << " ";} cout << endl;
    // cout << "Res : ";
    // for (auto i : res){
    //     cout << i << ",";
    // }
    // cout << endl;

    //if there is only one cluster, return that
    if (listOfGroups.size() <= 1){
        return vector<int> (clusters.size(), 1);
    }

    for (auto read = 0 ; read < clusters.size() ; read++){

        if (clusters[read] < 0 || presentReads[read] > 0){
            clusters[read] = -1;
        }
    }


    //now most reads should be assigned to a cluster. Rescue those that have not been assigned or assigned to a rare cluster

    for (auto read = 0 ; read < clusters.size() ; read++){

        if (clusters[read] < 0){
            clusters[read] = -1;
        }

        if (listOfGroups.find(clusters[read]) == listOfGroups.end() && presentReads[read] == 0){ //this means the read needs to be rescued 

            //cout << "rescuing "; for(auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++) {cout << allPartitions[binary][read];} cout << endl;
            //iterate through the list of groups and choose the one that can be explained by the less mistakes
            int bestGroup = -1;
            double bestGroupScore = 0;

            //compute a score for assigning the read to each existing group, then assign the read to the best group
            for (auto group : listOfGroups){

                double score = 0;
                int groupDecomposition = group;
                int readDecomposition = clustersAll[read];

                for (int binary = listOfFinalPartitions.size() -1 ; binary > -1  ; binary--){
                    int expectedAssignation = groupDecomposition%2;
                    int actualAssignation = readDecomposition%2;

                    if (actualAssignation != expectedAssignation && allPartitions[binary][read] != 0){
                        score -= allConfidences[binary][read]-0.5;
                    }
                    else if (actualAssignation == expectedAssignation){
                        score += allConfidences[binary][read]-0.5;
                    }

                    groupDecomposition /= 2;
                    readDecomposition /= 2;
                }
                if (score > bestGroupScore){
                    bestGroup = group;
                    bestGroupScore = score;
                }
            }
            clusters[read] = bestGroup;

            // cout << "Rescuing ";
            // for (int binary = 0 ; binary < listOfFinalPartitions.size()  ; binary++){
            //     cout << listOfFinalPartitions[binary].getMore()[read]<< "/" << listOfFinalPartitions[binary].getLess()[read] << " ";
            // }
            // cout << " as : " << bestGroup << endl;
            
        }
    }

    // cout << "Res of thread haplotypes: ";
    // for (auto i : res){
    //     cout << i << ",";
    // }
    // cout << endl;

    return clusters;

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

//input : the list of threaded clusters and snps
//output : reassign ALL reads to the cluster where they fit best
vector<int> rescue_reads(vector<int> &threadedClusters, vector<Column> &snps, vector<size_t> &suspectPostitions){

    if(DEBUG){cout << "Rescuing reads\r" << endl;}

    //list all the clusters in a map
    std::unordered_map <int, int> clusterIdx;
    int idx = 0;
    for (auto clust : threadedClusters){
        if (clust != -1){
            if (clusterIdx.find(clust) == clusterIdx.end()) {
                clusterIdx[clust] = idx;
                idx++;
            }
        } 
    }

    if (clusterIdx.size() == 0){
        return vector<int> (threadedClusters.size(), 0);
    }

    robin_hood::unordered_flat_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['-'] = 4;
    bases2content['?'] = 5;

    //create a vector counting for each read what cluster fits best
    vector<vector<int>> bestClusters (threadedClusters.size(), vector<int> (clusterIdx.size(), 0));

    //now iterate through the suspect positions
    for (auto position : suspectPostitions){

        //look at what base is normal for each cluster
        vector<vector<int>> basesForEachCluster(clusterIdx.size(), vector<int> (5, 0));
        int c = 0;
        for (auto read : snps[position].readIdxs){
            if (threadedClusters[read] != -1 && snps[position].content[c] != '?'){
                basesForEachCluster[clusterIdx[threadedClusters[read]]][bases2content[snps[position].content[c]]] += 1;
            }
            c++;
        }

        vector<char> clusterBase (clusterIdx.size(), 0);
        bool sure = true; //bool marking if all cluster agree within themselves
        for (auto c = 0 ; c < clusterBase.size() ; c++){
            char bestBase = '?';
            int bestBaseNb = 0;
            int totalBaseNb = 0;
            for (auto b = 0 ; b < 5 ; b++){
                int thisBaseNb = basesForEachCluster[c][b]; 
                totalBaseNb += thisBaseNb;
                if (thisBaseNb > bestBaseNb){
                    bestBaseNb = thisBaseNb;
                    bestBase = "ACGT-"[b];
                }
            }
            clusterBase[c] = bestBase;
            if (float(bestBaseNb)/totalBaseNb < 0.8){
                sure = false;
                break; //this position is not worth looking at
            }
        }

        //now each cluster has its base, let's update bestClusters
        if (sure){
            int c =0;
            for (auto read :snps[position].readIdxs){
                for (auto clust = 0 ; clust < clusterIdx.size() ; clust++){
                    if (snps[position].content[c] == clusterBase[clust]){
                        bestClusters[read][clust] += 1;
                    }
                    else{
                        bestClusters[read][clust] -= 1;
                    }
                }
                c++;
            }
        }
    }

    //now for each read we know at how many positions it agrees with each cluster : find best cluster
    vector<int> newClusters (threadedClusters.size(), -1);
    for (auto r = 0 ; r < newClusters.size() ; r++){

        // cout << "for read " << r << ", here is the bestCluster :"; for(auto i : bestClusters[r]){cout << i << ",";} cout << endl; 
        auto maxIterator = std::max_element(bestClusters[r].begin(), bestClusters[r].end());
        if (*maxIterator > 0){
            newClusters[r] = std::distance(bestClusters[r].begin() , maxIterator );
        }
        // else {
        //     cout << "wow, this read " << r << " has 0 positions, sad sad sad " << r << endl;
        // }
    }

    return newClusters;
}

/**
 * @brief Computes the consensus for each part of each partition of the contig using snps
 * 
 * @param contig 
 * @param partition Partition only on this contig
 * @param allreads 
 * @param allOverlaps 
 * @param snps 
 * @param partitions All the partitions, to be updated
 */
void compute_consensus_in_partitions(long int contig, vector<pair<pair<int,int>, vector<int>> > &partition,  vector <Read> &allreads, vector <Overlap> &allOverlaps, 
    vector<Column> &snps, robin_hood::unordered_map<int, int> &insertionPositions,
    std::unordered_map<unsigned long int, vector< pair<pair<int,int>, pair<vector<int>, std::unordered_map<int, std::string>>  > >> &partitions){


    //define map bases2content
    robin_hood::unordered_flat_map<char, short> bases2content = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'-', 4}, {'?', 5}};
    //get the sequence of contig
    string contigSequence = allreads[contig].sequence_.str();

    int n = 0;
    for (auto interval : partition){

        std::unordered_map <int, string> consensus_sequences; //consensus sequence of group int on this interval

        //be careful around homopolymers, there may be several equivalent alignments
        char currentBase = contigSequence[interval.first.first];
        std::unordered_map<int, int> numberOfDeletionsInHomopolymer;
        for (auto cluster: interval.second){
            numberOfDeletionsInHomopolymer[cluster] = 0;
        }

        //create a map associating each cluster to the number of reads of this cluster
        std::unordered_map<int, int> numberOfReadsInCluster;
        for (auto cluster: interval.second){
            numberOfReadsInCluster[cluster] = 0;
        }
        for (auto cluster: interval.second){
            numberOfReadsInCluster[cluster] += 1;
        }

        //iterate through the snps in this interval
        for (auto position = interval.first.first ; position <= interval.first.second ; position++){

            std::unordered_map <int, vector<int>> basesForEachCluster; //count the number of bases for each cluster on this interval
            for (auto cluster: interval.second){
                basesForEachCluster[cluster] = vector<int> (5, 0);
            }

            //look at what base is normal for each group
            for (int r = 0 ; r < snps[position].readIdxs.size() ; r++){
                int read = snps[position].readIdxs[r];
                
                if (snps[position].content[r] != '?'){
                    basesForEachCluster[interval.second[read]][bases2content[snps[position].content[r]]] += 1;
                }
            }

            if (contigSequence[position] == currentBase){ //we are in a homopolymer
                //if we are in a homopolymer, we need to count the number of deletions in the homopolymer (so that it can maybe compensate an insertion)
                for (auto cluster : numberOfDeletionsInHomopolymer){
                    numberOfDeletionsInHomopolymer[cluster.first] += basesForEachCluster[cluster.first][bases2content['-']];
                    //cluster size is the number of reads for which interval.second[read] == cluster.first
                    int clusterSize = numberOfReadsInCluster[cluster.first];
                    //if the number of deletions is too high, this base is a deletion (even though the deletions may be spread over several positions of the homopolymer)
                    if (numberOfDeletionsInHomopolymer[cluster.first] > clusterSize/2){
                        basesForEachCluster[cluster.first][bases2content['-']] += clusterSize; //making sure a deletion is chosen
                        numberOfDeletionsInHomopolymer[cluster.first] -= clusterSize;
                    }
                }
            }
            else{
                //set the number of deletions in homopolymer to the deletions at this position
                for (auto cluster : numberOfDeletionsInHomopolymer){
                    numberOfDeletionsInHomopolymer[cluster.first] = basesForEachCluster[cluster.first][bases2content['-']];
                }
            }

            //for all clusters, also add the base of contigSequence to basesForEachCluster to be sure there are no empty positions
            for (auto cluster : basesForEachCluster){
                basesForEachCluster[cluster.first][bases2content[contigSequence[position]]] += 3;
            }

            //find the index of highest base in basesForEachCluster
            for (auto cluster : basesForEachCluster){
                auto maxIterator = std::max_element(cluster.second.begin(), cluster.second.end());
                auto d = std::distance(cluster.second.begin() , maxIterator );
                if (d != 4){ //do not add a base if it is a gap
                    consensus_sequences[cluster.first] += "ACGT"[std::distance(cluster.second.begin() , maxIterator )];
                    // cout << "adding " << "ACGT"[std::distance(cluster.second.begin() , maxIterator )] << " to cluster " << cluster.first << " " << position << endl;
                }
            }

            //also iterate through the insertions 
            int i = 0;
            while (insertionPositions.find(10000*position+i) != insertionPositions.end() && i<10000){
                int pos = insertionPositions[10000*position+i];
                for (auto cluster: interval.second){
                    basesForEachCluster[cluster] = vector<int> (5, 0);
                }
                //look at what base is normal for each group
                for (int r = 0 ; r < snps[pos].readIdxs.size() ; r++){
                    int read = snps[pos].readIdxs[r];
                    
                    if (snps[pos].content[r] != '?'){
                        basesForEachCluster[interval.second[read]][bases2content[snps[pos].content[r]]] += 1;
                    }
                }
                //find the index of highest base in basesForEachCluster
                for (auto cluster : basesForEachCluster){
                    auto maxIterator = std::max_element(cluster.second.begin(), cluster.second.end());
                    auto d = std::distance(cluster.second.begin() , maxIterator );
                    if (d != 4){ //do not add a base if it is a gap
                        consensus_sequences[cluster.first] += "ACGT"[std::distance(cluster.second.begin() , maxIterator )];
                        // cout << "inserfqdqseting pheraps : " << "ACGT"[std::distance(cluster.second.begin() , maxIterator )] << " to cluster " << cluster.first << " " << position << endl;
                    }
                }
                i+=1;
            }
        }
        //print consensus_sequences
        // for (auto cluster : consensus_sequences){
        //     cout << "consensus sequence: " << cluster.first << " : " << cluster.second << endl;
        // }
        //add the consensus sequences to the partitions
        partitions[contig].push_back(make_pair(interval.first, make_pair(interval.second, consensus_sequences)));
        n += 1;

        // if (n == partition.size()){
        //     cout << "here are the partitiaeraezons for contig " << allreads[contig].name << " : " << endl;
        //     //print interval.second
        //     for (auto a : interval.second){
        //         cout << a << " ";
        //     }
        //     cout << endl;
        //     cout << "done with contig " << allreads[contig].name << " " << consensus_sequences[2] << endl;
        //     while(true){};
        // }

    }
}







