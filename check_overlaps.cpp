#include "check_overlaps.h"
// #include "spoa/spoa.hpp"
#include "cluster_graph.h"
//#include "WFA2-lib/bindings/cpp/WFAligner.hpp"
//using namespace wfa;

#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <fstream>
#include <chrono>

#include "robin_hood.h"
#include "input_output.h"

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

using namespace std::chrono;

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
void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, 
        vector<unsigned long int> &backbones_reads, unordered_map <unsigned long int ,vector<pair<pair<int,int>, vector<int>>>> &partitions, 
        bool assemble_on_assembly, unordered_map <int, std::pair<int,int>> &clusterLimits, unordered_map <int, vector<pair<int,int>>> &readLimits){

    //main loop : for each backbone read, build MSA (Multiple Sequence Alignment) and separate the reads
    int index = 0;
    for (unsigned long int read : backbones_reads){
        
        if (allreads[read].neighbors_.size() > 20 && (true || allreads[read].get_links_left().size()>0 || allreads[read].get_links_right().size()>0) 
             && allreads[read].name == "edge_159"){

            cout << "Looking at backbone read number " << index << " out of " << backbones_reads.size() << " (" << allreads[read].name << ")" << endl;

            vector<Column> snps;  //vector containing list of position, with SNPs at each position
            //first build an MSA
            cout << "Generating MSA" << endl;
            string truePar; //for debugging
            float meanDistance = generate_msa(read, allOverlaps, allreads, snps, partitions.size(), truePar, assemble_on_assembly, readLimits);

            //then separate the MSA
            cout << "Separating reads" << endl;
            vector<pair<pair<int,int>, vector<int>> > par = separate_reads(read, allOverlaps, allreads, snps, meanDistance, 
                                                                allreads[read].neighbors_.size()+1-int(assemble_on_assembly), clusterLimits);
            
            cout << "True partition : " << endl;
            cout << truePar << endl;
            
            // auto par = truePar.getPartition();//DEBUG
            // for (auto i = 0 ; i < par.size() ; i++) {par[i]++; }
            // cout << "Proposed partition : " << endl;
            // for (auto i = 0 ; i < par.size() ; i++){cout << par[i];}cout << endl;
            // cout << endl;

            partitions[read] = par;
    
        }
        index++;
    }
}

//input: a read with all its neighbor
//outputs : the precise alignment of all reads against input read in the form of matrix snps, return the mean editDistance/lengthOfAlginment
float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, 
    std::vector<Column> &snps, int backboneReadIndex, string &truePar, bool assemble_on_assembly, 
    unordered_map <int, vector<pair<int,int>>> &readLimits){

    // cout << "neighbors of read " << allreads[read].name << " : " << allreads[read].sequence_.str() << endl;
    //go through the neighbors of the backbone read and align it

    //keep count of the distance between two reads to know the mean distance
    float totalDistance = 0;
    double totalLengthOfAlignment = 0;
    
    //first enter the sequence of the read in snps, if the backbone read is a read, not if it's a contig

    // if (!assemble_on_assembly){
    //     string readSequence = allreads[read].sequence_.str();
    //     for (auto i = 0 ; i < readSequence.size() ; i++){
    //         snps[i][snps[i].size()-1] = readSequence[i];
    //     }
    // }

    //to remember on what backbone the backbone read is leaning -> itself
    allreads[read].new_backbone(make_pair(backboneReadIndex, allreads[read].neighbors_.size()), allreads[read].neighbors_.size()+1);
    string read_str = allreads[read].sequence_.str();

    //small loop to compute truePartition DEBUG
    for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){

        long int neighbor = allreads[read].neighbors_[n];
        Overlap overlap = allOverlaps[neighbor];
        
        if (overlap.sequence1 == read){

            truePar.push_back(allreads[overlap.sequence2].name[1]);
            // cout << "name : " << allreads[overlap.sequence2].name << " " << allreads[overlap.sequence2].name[1] << " " << truePartition[truePartition.size()-1] << endl;
        }
        else{
            truePar.push_back(allreads[overlap.sequence1].name[1]);
            // cout << "name : " << allreads[overlap.sequence1].name << " " << allreads[overlap.sequence1].name[1] << " " << truePartition[truePartition.size()-1] << endl;
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
    for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){
        long int neighbor = allreads[read].neighbors_[n];
        Overlap overlap = allOverlaps[neighbor];
        // cout << "overlap : " << overlap.position_1_1 << " " << overlap.position_1_2 << " " << overlap.position_2_1 << " "
        //     << overlap.position_2_2 << endl;
        if (overlap.sequence1 == read){
            // cout << "Neighbor of " << allreads[overlap.sequence1].name << " : " << allreads[overlap.sequence2].name << endl;
            allreads[overlap.sequence2].new_backbone(make_pair(backboneReadIndex,n), allreads[read].neighbors_.size()+1);
            if (overlap.strand){
                polishingReads.push_back(allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1, overlap.position_2_2-overlap.position_2_1).str());
                positionOfReads.push_back(make_pair(overlap.position_1_1, overlap.position_1_2));
            }
            else{
                polishingReads.push_back(allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1, overlap.position_2_2-overlap.position_2_1).reverse_complement().str());
                positionOfReads.push_back(make_pair(overlap.position_1_1, overlap.position_1_2));
            }
        }
        else {
            // cout << "Neighbor of " << allreads[overlap.sequence2].name << " : " << allreads[overlap.sequence1].name << endl;
            allreads[overlap.sequence1].new_backbone(make_pair(backboneReadIndex,n), allreads[read].neighbors_.size()+1);
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
    readLimits[read] = positionOfReads;

    // string consensus = consensus_reads(read_str , polishingReads);
    string consensus = read_str; //DEBUG
    cout << "Done building a consensus of the backbone" << endl;

    //snps = vector<vector<char>>(consensus.size(), vector<char>(polishingReads.size(), '?'));
    snps = vector<Column>(consensus.size());

    //while all the alignments are computed, build the positions
    robin_hood::unordered_map<int, int> insertionPos;
    vector<int> numberOfInsertionsHere (consensus.size()+1, 0);

    float alignmentTime = 0;
    float MSAtime = 0;

    vector<float> mappingQuality; //DEBUG

    for (auto n = 0 ; n < polishingReads.size() ; n++){

        cout << "Aligned " << n << " reads out of " << allreads[read].neighbors_.size() << " on the backbone\r";

        auto t1 = high_resolution_clock::now();

        EdlibAlignResult result = edlibAlign(polishingReads[n].c_str(), polishingReads[n].size(),
                                    consensus.substr(positionOfReads[n].first, positionOfReads[n].second-positionOfReads[n].first).c_str(),
                                    positionOfReads[n].second-positionOfReads[n].first,
                                    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

        string alignment;
        char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
        for (int l = 0; l < result.alignmentLength; l++) {
            alignment += moveCodeToChar[result.alignment[l]];
        }

        // if (n > 1500){
        //     cout << "Aligning read of length " << polishingReads[n].size() << " : " << polishingReads[n].substr(0,100) << endl;
        //     cout << alignment.substr(0,100) << " " << result.startLocations[0]+positionOfReads[n].first << endl;
        // }

        // cout << "Alignment : " << endl;
        // cout << alignment << endl;
        // cout << allOverlaps[allreads[read].neighbors_[n]].CIGAR << endl;

        auto t2 = high_resolution_clock::now();
        
        // string ref = consensus.substr(positionOfReads[n].first, positionOfReads[n].second-positionOfReads[n].first);
        // aligner.alignEnd2End(polishingReads[n],ref);

        auto t2_5 = high_resolution_clock::now();

        totalLengthOfAlignment += result.alignmentLength;
        totalDistance += result.editDistance;
        mappingQuality.push_back(float(result.editDistance)/result.alignmentLength);
        // cout << "Alignment distance : " << float(result.editDistance)/result.alignmentLength << endl;

        // if (n == 10) {break;}

        //a loop going through the CIGAR and modifyning snps
        int indexQuery = result.startLocations[0]+positionOfReads[n].first; //query corresponds to consensus
        int indexTarget = 0; //target corresponds to the read
        int numberOfInsertionsThere = 0;

        //cout << "beginning of query : " << indexQuery << " " << consensus.size() << " " << snps.size() << " " << result.alignmentLength << endl;
        
        // if (n == 169){
        // for (int i = 0; i < result.alignmentLength; i++){cout << moveCodeToChar[result.alignment[i]];} cout << endl;}
        
        for (int l = 0; l < alignment.size(); l++) {

            if (indexQuery < consensus.size()){

                if (alignment[l] == '=' || alignment[l] == 'X' || alignment[l] == 'M'){
                    //fill inserted columns with '-' just before that position
                    for (int ins = numberOfInsertionsThere ; ins < numberOfInsertionsHere[indexQuery] ; ins++){ //in these positions, insert '-' instead of '?'
                        snps[insertionPos[10000*indexQuery+ins]].readIdxs.push_back(n);
                        snps[insertionPos[10000*indexQuery+ins]].content.push_back('-');
                    }

                    snps[indexQuery].readIdxs.push_back(n);
                    snps[indexQuery].content.push_back(polishingReads[n][indexTarget]);

                    indexQuery++;
                    indexTarget++;
                    numberOfInsertionsThere = 0;
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
                }
                else if (alignment[l] == 'I'){ //hardest one
                    if (numberOfInsertionsHere[indexQuery] <= 9999 && indexQuery > positionOfReads[n].first) {

                        if (numberOfInsertionsThere >= numberOfInsertionsHere[indexQuery]) { //i.e. this is a new column
                            insertionPos[10000*indexQuery+numberOfInsertionsHere[indexQuery]] = snps.size();
                            numberOfInsertionsHere[indexQuery] += 1;
                            
                            Column newInsertedPos;
                            newInsertedPos.readIdxs = snps[indexQuery-1].readIdxs;    
                            newInsertedPos.content = vector<char>(snps[indexQuery-1].content.size() , '-');
                            newInsertedPos.content[newInsertedPos.content.size()-1] = polishingReads[n][indexTarget];
                            snps.push_back(newInsertedPos);
                        }
                        else{
                            snps[insertionPos[10000*indexQuery+numberOfInsertionsThere]].readIdxs.push_back(n);
                            snps[insertionPos[10000*indexQuery+numberOfInsertionsThere]].content.push_back(polishingReads[n][indexTarget]);
                        }
                        numberOfInsertionsThere ++;
                    }
                    indexTarget++;
                }
            }
            //cout << l << " " << result.alignmentLength <<  "\n";
        }

        auto t3 = high_resolution_clock::now();

        alignmentTime += duration_cast<milliseconds>(t2-t1).count();
        MSAtime += duration_cast<milliseconds>(t3-t2).count();
        // cout << "index target : " << indexTarget<< endl;
        // cout << "index query : " << indexQuery << endl;
        // cout << "consensus length : " << consensus.size() << endl;
        
        edlibFreeAlignResult(result);
    }

    cout << "Building MSA took time... " << alignmentTime << " for edlib and " << MSAtime << " for filling the vector" << endl;
    
    //print snps (just for debugging)
    // int step = 1;
    // int prop = 10; //1 for every base
    // int firstRead = 1450;
    // int lastRead = polishingReads.size();
    // int numberOfReads = lastRead-firstRead;
    // int start = 011;
    // int end = 1060;
    // vector<string> reads (numberOfReads);
    // string cons = "";
    // for (unsigned short i = start ; i < end; i+=prop){
        
    //     for (short n = 0 ; n < numberOfReads*step ; n+= step){
    //         char c = '?';
    //         int ri = 0;
    //         int soughtRead = firstRead+n;
    //         for (auto r : snps[i].readIdxs){
    //             if (r == soughtRead){
    //                 c = snps[i].content[ri];
    //             }
    //             ri ++;
    //         }
    //         reads[n/step] += c;
    //     }
    //     // for (short insert = 0 ; insert < min(9999,numberOfInsertionsHere[i]) ; insert++ ){
    //     //     int snpidx = insertionPos[10000*i+insert];
    //     //     for (short n = 0 ; n < numberOfReads*step ; n+= step){
    //     //         char c = '?';
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
    //     if (neighbor[0] != '?'){
    //         cout << neighbor << " " << index  << endl;
    //     }
    //     index+= step;
    // }
    // int n = 0;
    // for(auto i : consensus.substr(start, end-start)){
    //     if (n%prop == 0){
    //         cout << i;
    //     }
    //     n+=1;
    // } cout << endl;

    cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;

    return totalDistance/totalLengthOfAlignment;

    //*/
}

//input : a backbone read with a list of all reads aligning on it
//output : a polished backbone read
//function using Racon to consensus all reads
string consensus_reads(string &backbone, vector <string> &polishingReads){
    
    system("mkdir tmp/ 2> trash.txt");
    std::ofstream outseq("tmp/unpolished.fasta");
    outseq << ">seq\n" << backbone;
    outseq.close();

    std::ofstream polishseqs("tmp/reads.fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        polishseqs << ">read"+std::to_string(read)+"\n" << polishingReads[read] << "\n";
    }
    polishseqs.close();

    system("minimap2 -t 1 -x map-ont tmp/unpolished.fasta tmp/reads.fasta > tmp/mapped.paf 2>tmp/trash.txt");
    system("racon -e 1 -t 1 tmp/reads.fasta tmp/mapped.paf tmp/unpolished.fasta > tmp/polished.fasta 2>tmp/trash.txt");

    std::ifstream polishedRead("tmp/polished.fasta");
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

//input : a list of reads
//output : an assembled contig
//function using miniasm+minimap to build an assembly
string local_assembly(vector <string> &reads){

    std::ofstream polishseqs("tmp/reads.fasta");
    for (int read=0 ; read < reads.size() ; read++){
        polishseqs << ">read"+std::to_string(read)+"\n" << reads[read] << "\n";
    }
    polishseqs.close();

    system("minimap2 -x ava-ont -t1 tmp/reads.fasta tmp/reads.fasta > tmp/local_overlaps.paf 2> tmp/trash.txt");
    string command = "miniasm -1 -2 -f tmp/reads.fasta tmp/local_overlaps.paf 2> tmp/trash.txt | awk \'{if ($1 == \"S\") print $3;}\' > tmp/local_assembly.txt 2> tmp/trash.txt";
    cout << command << endl;
    system(command.c_str());

    std::ifstream newcontig("tmp/local_assembly.txt");
    string line;
    string contig;
    while(getline(newcontig, line)){
        contig = line; 
    }
    return contig;
}

//input : a set of reads aligned to read in matrix snps
//output : reads separated by their region of origin
vector<pair<pair<int,int>, vector<int>> > separate_reads(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, 
        std::vector<Column> &snps, float meanDistance, int numberOfReads, unordered_map <int, std::pair<int,int>> &clusterLimits){

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

    vector<vector<distPart>> distanceBetweenPartitions;
    vector<size_t> suspectPostitions;
    vector <int> interestingParts; //DEBUG

    for (int position = 0 ; position < snps.size() ; position++){ 


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

        float threshold = 1 + numberOfReadsHere*meanDistance/2 + 3*sqrt(numberOfReadsHere*meanDistance/2*(1-meanDistance/2));
        //threshold = 3; //DEBUG
        if (content[4] < 0.25*float(numberOfReadsHere) //not too many '-', because '-' are less specific than snps
            && *std::max_element(content, content+4) < numberOfReadsHere-content[4]-threshold){ //this position is suspect
            // cout << threshold << " " << position << " ;bases : " << content[0] << " " << content[1] << " " << content[2] << " " << content[3] << " " << content[4] << endl;
            suspectPostitions.push_back(position);
            //go through the partitions to see if this suspicious position looks like smt we've seen before
            vector<distPart> distances (distanceBetweenPartitions.size());
            bool found = false;
            for (auto p = 0 ; p < partitions.size() ; p++){

                distancePartition dis = distance(partitions[p], snps[position]);
                auto comparable = min(dis.n00,dis.n11) + dis.n01 + dis.n10;
                // if (comparable > 10){
                //     cout << "comparable : " << float(dis.n01+dis.n10)/(dis.n00+dis.n11+dis.n01+dis.n10) << " " << dis.augmented<< " mean distance : " << meanDistance << endl;
                //     partitions[p].print();
                //     Partition(snps[position]).print();
                    
                // }

                if (float(dis.n01+dis.n10)/(dis.n00+dis.n11+dis.n01+dis.n10) <= meanDistance && dis.augmented && comparable > min(10.0, 0.3*numberOfReads)){
                    
                    int pos = -1;
                    if (position < allreads[read].size()){
                        pos = position;
                    }

                    // cout << "augmenting " << p << " now : " << partitions[p].number() << " " << position << " "<< endl;
                    // if (p == 0){
                    //     interestingParts.push_back(position);
                    // }

                    partitions[p].augmentPartition(dis.partition_to_augment, pos);
                    found = true;
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

            if (!found && position < allreads[read].size()){    // the second condition is here to create partitions only at specific spots of the backbone
                distanceBetweenPartitions.push_back(distances);
                partitions.push_back(Partition(snps[position], position)); 
            }
            
            numberOfSuspectPostion += 1;

            //two suspect positions next to each other can be artificially correlated through alignement artefacts
            position += 5;
        }
    }

    cout << "found " << numberOfSuspectPostion << " suspect positions" << endl;


    //for debugging only
    //outputMatrix(snps, suspectPostitions, std::to_string(read));


    if (partitions.size() == 0){ //there are no position of interest
        vector<pair<pair<int,int>, vector<int>>> e;
        return e;
    }

    // cout << "Outputting the graph" << endl;

    /* looking at the graph
    //square the distanceBetweenPartions matrix (which is a triangle matrix until now)

    for (int i = 0 ; i < distanceBetweenPartitions.size() ; i++){
        distPart diag;
        distanceBetweenPartitions[i].push_back(diag); //that is the diagonal
        while(distanceBetweenPartitions[i].size() < distanceBetweenPartitions.size()){
            distanceBetweenPartitions[i].push_back(diag);
        }
    }
    for (int i = 0 ; i < distanceBetweenPartitions.size() ; i++){
       for (int j = i+1 ; j < distanceBetweenPartitions.size() ; j++){
            distanceBetweenPartitions[i][j].distance = distanceBetweenPartitions[j][i].distance;
        }
    }

    // build the adjacency matrix
    vector<vector <float>> adj (distanceBetweenPartitions.size(), vector<float> (distanceBetweenPartitions.size(), 0));
    //each node keeps only the links to very close elements
    int numberOfNeighborsKept = 5;
    for (auto i = 0 ; i < distanceBetweenPartitions.size() ; i++){
        if (distanceBetweenPartitions[i].size() > numberOfNeighborsKept){
            vector <distPart> minElements (numberOfNeighborsKept);
            std::partial_sort_copy(distanceBetweenPartitions[i].begin(), distanceBetweenPartitions[i].end(), minElements.begin(), minElements.end(), comp);
            float maxDiff = minElements[numberOfNeighborsKept-1].distance;

            //cout << "maxdiff of partition " << i << " : " << maxDiff << endl;
            if (i == 1){
                // cout << "1 : " << endl;
                // for(auto j:adj[i]){cout<<j << " ";}
                // cout << endl;
            }
            
            int numberOf1s = 0;
            for (int j = 0 ; j < distanceBetweenPartitions[i].size() ; j++){
                if (distanceBetweenPartitions[i][j].distance > maxDiff || distanceBetweenPartitions[i][j].distance == 1 || numberOf1s >= numberOfNeighborsKept){
                    
                }
                else {
                    if (distanceBetweenPartitions[i][j].distance < meanDistance){
                        // cout << "distance : " << distanceBetweenPartitions[i][j].distance << endl;
                        adj[i][j] = 1;
                        adj[i][j] = 1;
                        numberOf1s++;
                    }
                }
            }
        }
    }

    */

   /* to look at the graph

    vector<int> clusters (partitions.size());

    cluster_graph_chinese_whispers(adj, clusters);
    outputGraph(adj, clusters, "graph.gdf");

    //filter out too small clusters
    vector <int> sizeOfCluster;
    for (auto p = 0 ; p < partitions.size() ; p++){
        while (clusters[p] >= sizeOfCluster.size() ){ 
            sizeOfCluster.push_back(0);
        }
        sizeOfCluster[clusters[p]] += 1;
    }
    vector<Partition> listOfFinalPartitionsdebug(sizeOfCluster.size(), Partition(partitions[0].size()));

    //now merge each cluster in one partition
    for (auto p = 0 ; p < partitions.size() ; p++){

        if (sizeOfCluster[clusters[p]] > 3){
            // if (sizeOfCluster[clusters[p]] == 41){
            //     partitions[p].print();
            // }
            listOfFinalPartitionsdebug[clusters[p]].mergePartition(partitions[p]);
        }
    }

    //filter out empty partitions
    vector<Partition> listOfFinalPartitions2;
    for (auto p : listOfFinalPartitionsdebug){
        if (p.number()> 0 && p.isInformative(meanDistance/2, true)) {
            listOfFinalPartitions2.push_back(p);
        }
    }
    listOfFinalPartitionsdebug = listOfFinalPartitions2;

    for (auto p : listOfFinalPartitionsdebug){
        cout << "final : " << endl;
        p.print();
    }

    */
    
    float threshold =  max(4.0, min(0.01*numberOfSuspectPostion, 0.001*snps.size()));

    vector<Partition> listOfFinalPartitions;
    for (auto p1 = 0 ; p1 < partitions.size() ; p1++){

        // if (partitions[p1].number() > 2){
        //     cout << "non informative partition : "  << threshold << " " << numberOfSuspectPostion << " " << snps.size()<< endl;
        //     partitions[p1].print();
        // }
        
        if (partitions[p1].number() > threshold && partitions[p1].isInformative(meanDistance/2, true)){

            // cout << "informative partition 2 : " << endl;
            // partitions[p1].print();

            bool different = true;
            
            for (auto p2 = 0 ; p2 < listOfFinalPartitions.size() ; p2++){

                distancePartition dis = distance(listOfFinalPartitions[p2], partitions[p1], 2);
                if (dis.augmented){
                    Partition newPart = listOfFinalPartitions[p2];
                    newPart.mergePartition(partitions[p1], dis.phased);

                    //see if confidence is improved by merging the two partitions, meaning differences were shaky
                    if (dis.n01+dis.n10 < 0.1*(dis.n00+dis.n11) || newPart.compute_conf() > listOfFinalPartitions[p2].compute_conf()){
                        
                        // cout << endl << "Now merging " << dis.n00 << " " << dis.n01 << " " << dis.n10 << " " << dis.n11 << endl;
                        // listOfFinalPartitions[p2].print();
                        // partitions[p1].print();
                        // cout << "phasing : " << dis.phased << endl << endl; 

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

    if (listOfFinalPartitions.size() == 0){
        return vector<pair<pair<int,int>, vector<int>>> (0);
    }

    // cout << "selecting partitions" << endl;
    vector<Partition> listOfFinalPartitionsTrimmed = select_partitions(listOfFinalPartitions, numberOfReads, meanDistance/2);


    // cout << "threading clusters" << endl;
    //now aggregate all those binary partitions in one final partition. There could be up to 2^numberBinaryPartitions final groups
    vector<pair<pair<int,int>, vector<int>>> threadedClusters = threadHaplotypes(listOfFinalPartitionsTrimmed, numberOfReads, clusterLimits);

    // cout << "threaded clusters ! " << endl;
    // for (auto i = 0 ; i < threadedClusters.size() ; i++){cout << threadedClusters[i];}cout << endl;

    // //rescue reads that have not been assigned to a cluster
    // vector<int> finalClusters = rescue_reads(threadedClusters, snps, suspectPostitions);

    //print snps (just for debugging)
    // int step = 1;
    // int prop = 1; //1 for every base
    // int numberOfDisplayedReads = numberOfReads;
    // int start = 011;
    // int end = 10000000;
    // vector<string> reads (numberOfDisplayedReads);
    // string cons = "";
    // auto n = 0;
    // for (auto i : interestingParts){
        
    //     if (i < end && n%prop == 0){

    //         for (short n = 0 ; n < numberOfDisplayedReads*step ; n+= step){
    //             char c = '?';
    //             char justBefore = '?';
    //             char justAfter = '?';
    //             int ri = 0;
    //             for (auto r : snps[i].readIdxs){
    //                 if (r == n){
    //                     c = snps[i].content[ri];
    //                 }
    //                 ri ++;
    //             }
                
    //             ri = 0;
    //             for (auto r : snps[i-1].readIdxs){
    //                 if (r == n){
    //                     justBefore = snps[i-1].content[ri];
    //                 }
    //                 ri ++;
    //             }

    //             ri = 0;
    //             for (auto r : snps[i+1].readIdxs){
    //                 if (r == n){
    //                     justAfter = snps[i+1].content[ri];
    //                 }
    //                 ri ++;
    //             }

    //             // if (justAfter != '-' && justAfter != '?'){
    //             //     justAfter += 32;
    //             // }
    //             // if (justBefore != '-' && justBefore != '?'){
    //             //     justBefore += 32;
    //             // }
    //             // reads[n/step] = reads[n/step] + "...";
    //             // reads[n/step] += justBefore;
    //             reads[n/step] += c;
    //             // reads[n/step] += justAfter;
    //         }
    //         // for (short insert = 0 ; insert < min(9999,numberOfInsertionsHere[i]) ; insert++ ){
    //         //     int snpidx = insertionPos[10000*i+insert];
    //         //     for (short n = 0 ; n < numberOfReads*step ; n+= step){
    //         //         char c = '?';
    //         //         int ri = 0;
    //         //         for (auto r : snps[snpidx].readIdxs){
    //         //             if (r == n){
    //         //                 c = snps[snpidx].content[ri];
    //         //             }
    //         //             ri ++;
    //         //         }
    //         //         reads[n/step] += c;
    //         //     }
    //         // }
    //     }
    //     n++;
    // }
    // cout << "Here are the aligned reads : " << endl;
    // int index = 0;
    // for (auto neighbor : reads){
    //     if (neighbor[4] != '?'){
    //         cout << neighbor << " " << index  << endl;
    //     }
    //     index++;
    // }

    // return finalClusters;
    return threadedClusters;
}

//input : one partition and one list of chars
//output : is the list of chars close to the partition ? If so, augment the partition. Return the chi-square of the difference in bases
distancePartition distance(Partition &par1, Column &par2){

    /*
    when computing the distance, there is not 5 letters but 2 : the two alleles, which are the two most frequent letters
    */
    distancePartition res;
    res.augmented = true;
    vector <int> idxs1 = par1.getReads();
    vector<short> part1 = par1.getPartition();
    vector<float> confs1 = par1.getConfidence();

    vector <int> idxs2 = par2.readIdxs;
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
    char mostFrequent2 = 'A';
    int maxFrequence2 = content2[0];
    char secondFrequent2 = 'C';
    int secondFrequence2 = content2[1];
    if (content2[0] < content2[1]){
        mostFrequent2 = 'C';
        maxFrequence2 = content2[1];
        secondFrequent2 = 'A';
        secondFrequence2 = content2[0];
    }
    for (auto i = 2 ; i < 5 ; i++){
        if (content2[i] > maxFrequence2){
            secondFrequence2 = maxFrequence2;
            secondFrequent2 = mostFrequent2;
            maxFrequence2 = content2[i];
            mostFrequent2 = "ACGT-"[i];
        }
        else if (content2[i] > secondFrequence2){
            secondFrequent2 = "ACGT-"[i];
            secondFrequence2 = content2[i];
        }
    }

    // cout << "Two most frequent : " << mostFrequent2 << "," << secondFrequent2 << " : ";
    // for (auto i : par2){cout << i;} cout << endl;

    float scores [2] = {0,0}; //the scores when directing mostFrequent on either mostfrequent2 or secondFrequent2
    float bestScore = 0;
    //remember all types of matches for the chi square test
    int matches00[2] = {0,0};
    int matches01[2] = {0,0};
    int matches10[2] = {0,0};
    int matches11[2] = {0,0};
    Column newPartitions [2];

    newPartitions[0].readIdxs = par2.readIdxs;
    newPartitions[1].readIdxs = par2.readIdxs;
    for (auto c : par2.content){
        if (c == mostFrequent2){
            newPartitions[0].content.push_back('A');
            newPartitions[1].content.push_back('a');
        }
        else if (c == secondFrequent2){
            newPartitions[0].content.push_back('a');
            newPartitions[1].content.push_back('A');
        }
        else{
            newPartitions[0].content.push_back('0');
            newPartitions[1].content.push_back('0');
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

            if (par2.content[n2] == mostFrequent2){

                if (part1[n1] == 1){
                    scores[0] += conf;
                    scores[1] -= conf;
                    bestScore += conf;
                    matches11[0] += 1;
                    matches10[1] += 1;
                }
                else if (part1[n1] == -1){
                    scores[0] -= conf;
                    scores[1] += conf;
                    bestScore += conf;
                    matches01[0] += 1;
                    matches00[1] += 1;
                }
            }
            else if (par2.content[n2] == secondFrequent2){

                if (part1[n1] == 1){
                    scores[0] -= conf;
                    scores[1] += conf;
                    bestScore += conf;
                    matches10[0] += 1;
                    matches11[1] += 1;
                }
                else if (part1[n1] == -1){
                    scores[0] += conf;
                    scores[1] -= conf;
                    bestScore += conf;
                    matches00[0] += 1;
                    matches01[1] += 1;
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

    //now look at the best scores

    auto maxScore = scores[0];
    auto maxScoreIdx = 0; //can be either 0 or 1
    for (auto i = 0 ; i < 2 ; i++){
        //cout << scores[i] << " , ";
        if (scores[i]>=maxScore){
            maxScore = scores[i];
            maxScoreIdx = i;
        }
    }

    res.n00 = matches00[maxScoreIdx];
    res.n01 = matches01[maxScoreIdx];
    res.n10 = matches10[maxScoreIdx];
    res.n11 = matches11[maxScoreIdx];
    res.score = scores[maxScoreIdx]/bestScore;
    res.phased = -2*maxScoreIdx + 1;
    res.partition_to_augment = newPartitions[maxScoreIdx];
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

//input : a list of partitions
//output : only the partitions that look very sure of themselves
vector<Partition> select_partitions(vector<Partition> &listOfFinalPartitions, int numberOfReads, float errorRate){
    vector<int> frequenceOfPart (pow(2, listOfFinalPartitions.size()));

    vector<vector<int>> allIdxs;
    for (auto i : listOfFinalPartitions){
        allIdxs.push_back(i.getReads());
    }
    vector<vector<short>> allPartitions;
    for (auto i : listOfFinalPartitions){
        allPartitions.push_back(i.getPartition());
    }
    vector<vector<float>> allConfidences;
    for (auto i : listOfFinalPartitions){
        allConfidences.push_back(i.getConfidence());
    }

    //as a first step, we'll try to throw away the partitions that do not look very sure of themselves
    /*to to that, we'll look at each read, see if it's well clustered in one partition. 
    If yes, all partitions where it's badly clustered are suspicious
    Indeed, we expect the different errors to compensate each other
    */
    

    vector<bool> readsClassified (numberOfReads, true);
    for (int par = 0 ; par < listOfFinalPartitions.size() ; par++){
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
    for (int par = 0 ; par < listOfFinalPartitions.size() ; par++){
        //compute a score evaluating the certainty of the partition
        listOfFinalPartitions[par].compute_conf();
    }

    struct {
        bool operator()(Partition a, Partition b) const { return a.get_conf() > b.get_conf(); }
    } customLess;
    std::sort(listOfFinalPartitions.begin(), listOfFinalPartitions.end(), customLess);

    cout << "Here are all the partitions, sorted : " << endl;
    for (auto p : listOfFinalPartitions){
        p.print();
    }

    //draw a list of compatible partitions
    
    vector<bool> compatibles(listOfFinalPartitions.size(), false);
    for (int p = 0 ; p < listOfFinalPartitions.size() ; p++){
        bool compatible = true;
        for (int p2 = 0 ; p2 < p ; p2++){
            if (compatibles[p2]){ //for all more sturdy validated partitions, check if this one is compatible
                int compCode = compatible_partitions(listOfFinalPartitions[p], listOfFinalPartitions[p2]);
                if (compCode == 0){
                    compatible = false;
                    break;
                }
                else if (compCode == 2){ //means that they're the same partition
                    compatible = false;
                    listOfFinalPartitions[p2].mergePartition(listOfFinalPartitions[p]);
                    break;
                }
            }
        }
        if (compatible){
            compatibles[p] = true;
        }
    }
    
    vector<Partition> compatiblePartitions;
    cout << "compatible partitions : " << endl;
    for (auto p = 0 ; p < listOfFinalPartitions.size(); p++){
        if (compatibles[p]){
            compatiblePartitions.push_back(listOfFinalPartitions[p]);
            listOfFinalPartitions[p].print();
        }
    }
    
    //Now see if some partitions are just not sure enough of themselves

    vector<bool> trimmedListOfFinalPartitionBool (compatiblePartitions.size(), false); //true if a partition is kept, false otherwise

    allIdxs = {};
    for (auto i : compatiblePartitions){
        allIdxs.push_back(i.getReads());
    }
    allPartitions = {};
    for (auto i : compatiblePartitions){
        allPartitions.push_back(i.getPartition());
    }
    allConfidences = {};
    for (auto i : compatiblePartitions){
        allConfidences.push_back(i.getConfidence());
    }

    for (int par = 0 ; par < compatiblePartitions.size() ; par++){

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
        double errors [2] = {(1-means[0]/n0)/(2-means[0]/n0-means[1]/n1) * errorRate*4 + errorRate 
                            , (1-means[1]/n1)/(2-means[0]/n0-means[1]/n1) * errorRate*4 + errorRate }; //tolerate on average errorRate*3
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

        cout << "filtering partition, tolerating " << errors[0] << "," << errors[1] << endl;
        compatiblePartitions[par].print();
        cout << "There are " << numberOfUnsureReads << " unsure reads, (" << numberOfUnsureReads0 << "+" << numberOfUnsureReads1 << ") out of " 
            << numberOfPartitionnedRead << " partitionned reads, and in terms of size, the het rate is " <<
        compatiblePartitions[par].number() << " for a length of " << (compatiblePartitions[par].get_right()-compatiblePartitions[par].get_left()) << endl;

        if (float(numberOfUnsureReads)/numberOfPartitionnedRead < 0.15 //then the partition is sure enough of itself 
            && compatiblePartitions[par].number() > 20/allIdxs[par].size()*100){ //and big enough

            trimmedListOfFinalPartitionBool[par] = true;
        }

    }

    vector<Partition> trimmedListOfFinalPartition;
    for (auto p = 0 ; p < listOfFinalPartitions.size() ; p++){
        if (trimmedListOfFinalPartitionBool[p]){
            trimmedListOfFinalPartition.push_back(compatiblePartitions[p]);
            // cout << "remaining partition : " << endl;
            // listOfFinalPartitions[p].print();
        }
    }

    return trimmedListOfFinalPartition;
}

//input : a list of all binary partitions found in the reads
//output : a vector defining intervals with their two bounds and their attached local partition
vector<pair<pair<int,int>, vector<int>> >threadHaplotypes(vector<Partition> &compatiblePartitions, int numberOfReads, unordered_map <int, std::pair<int,int>> &clusterLimits){


    //now go through the reference and give the repartition of reads at each position

    //sort the partitions by position left
    std::sort(compatiblePartitions.begin(), compatiblePartitions.end(), [](Partition& lhs, Partition& rhs) {
                    return lhs.get_left() < rhs.get_left();});

    vector<pair<pair<int,int>, set<int>>> intervals; //list of intervals, with all partitions on each

    vector<int> borders;
    for (auto p : compatiblePartitions){
        borders.push_back(p.get_left());
        borders.push_back(p.get_right());
    }
    std::sort(borders.begin(), borders.end());

    //list the limits of the intervals
    int lastRight = 0;
    for (auto b : borders){
        set <int> s;
        intervals.push_back(make_pair(make_pair (lastRight, b), s ));
        lastRight = b;
    }


    //now list what partition is present on which interval
    int n = 0;
    for (auto p : compatiblePartitions){
        int n2 = 0;
        for (auto interval : intervals){
            if (p.get_left() <= interval.first.first && p.get_right() >= interval.first.second){
                intervals[n2].second.emplace(n);
            }
            n2++;
        }
        n++;
    }

    
    //for each interval, determine the partition of reads
    vector<pair<pair<int,int>, vector<int>>> res;
    for (auto interval : intervals){
        vector<Partition> localPartitions;
        for (auto lp : interval.second){
            localPartitions.push_back(compatiblePartitions[lp]);
        }
        vector<int> localPartition = threadHaplotypes_in_interval(localPartitions, numberOfReads);

        res.push_back(make_pair(make_pair(interval.first.first, interval.first.second), localPartition));
    }

    cout << "here are the intervals : " << endl;
    n = 0;
    for (auto interval : intervals){
        cout << interval.first.first << " <-> " << interval.first.second << " : ";
        for (auto p : interval.second){cout << p << ",";} cout << endl;
        for (auto r : res[n].second){cout << r+1;} cout << endl;
        n++;
    }

    return res;
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
    auto idxs2 = p2.getReads();
    auto content2 = p2.getPartition();

    int n = 0;
    float numberOf1s = 0;
    float numberOf0s = 0;

    int one_in_1 = 0;
    int one_in_0 = 0;
    int zero_in_1 = 0;
    int zero_in_0 = 0;

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

    if (zero_in_0/numberOf0s > 0.9){
        p1.flipPartition();
        compatible += 1;
    }
    else if(zero_in_1/numberOf0s > 0.9){
        p1.flipPartition();
        p2.flipPartition();
        compatible += 1;
    }

    if (one_in_0/numberOf1s > 0.9){
        compatible += 1;
    }
    else if (one_in_1/numberOf1s > 0.9){
        p2.flipPartition();
        compatible += 1;
    }

    cout << "Compatibility : " << numberOf0s << " " << zero_in_0 << " " << zero_in_1 << " ; " << numberOf1s << " " << one_in_0 << " " << one_in_1 << " ; " << compatible << endl;
    return compatible;
}

//input : a list of all binary partitions found in the reads
//output : a single partition, with several number corresponding to several clusters
vector<int> threadHaplotypes_in_interval(vector<Partition> &listOfFinalPartitions, int numberOfReads){


    vector<vector<short>> allPartitions;
    for (auto i : listOfFinalPartitions){
        allPartitions.push_back(i.getPartition());
    }
    vector<vector<float>> allConfidences;
    for (auto i : listOfFinalPartitions){
        allConfidences.push_back(i.getConfidence());
    }

    int numberOfAssignedReads=0;

    vector<int> clusters (numberOfReads, -1); //this will contain only high-confidence reads
    vector<int> clustersAll (numberOfReads, -1); //this will contain all reads

    for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
        int c = 0;
        int lastRead = 0;
        for (auto read : listOfFinalPartitions[binary].getReads()){

            for (int absentread = lastRead ; absentread < read ; absentread++) { //means that the read is absent in this interval 
                clusters[absentread] = -pow(2, listOfFinalPartitions.size())-2;
            }
            lastRead = read+1;

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
            c++;
        }

        //the last stretch is also devoid of reads
        for (int absentread = lastRead ; absentread < numberOfReads ; absentread++) { //means that the read is absent in this interval 
            clusters[absentread] = -pow(2, listOfFinalPartitions.size())-2;
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
    
    //to make sure we don't over-estimate the number of clusters, we'll make the assumption that all haplotypes have haplotype-specific mutation => the number of final cluster cannot be higher than the number of partitions
    vector<float> proportionOf1;
    for (auto p : listOfFinalPartitions){
        proportionOf1.push_back(p.proportionOf1());  
    }
    vector<float> listOfLikelihoods;
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
        // cout << "group " << group.first << " is expected " << numberOfAssignedReads*proba << " times and arrives " << group.second << " times " << endl;
        listOfLikelihoods.push_back(float(frequenceOfPart[group] -numberOfAssignedReads*proba) / frequenceOfPart[group]);
    }

    float minLikelihood = -10;
    // if (listOfLikelihoods.size() > listOfFinalPartitions.size() && listOfFinalPartitions.size()>1){
    //     //vector <float> minElements (listOfLikelihoods.size()-listOfFinalPartitions.size());
    //     vector <float> minElements (1); //DEBUG
    //     std::partial_sort_copy(listOfLikelihoods.begin(),  listOfLikelihoods.end(), minElements.begin(), minElements.end());
    //     minLikelihood = minElements[minElements.size()-1];
    //     // cout << "likelihoods : " << endl;
    //     // for (auto a : listOfLikelihoods) {cout << a << ",";}cout << endl;
    //     // cout << "minElements : " << endl;
    //     // for (auto m : minElements) {cout << m << ",";} cout << endl;
    // }

    
    // cout << "min likelihood : " << minLikelihood << endl;

    //establish a list of possible clusters
    std::set <int> listOfGroups;
    int indexOfCount = 0;
    for (auto group : count){
        if (listOfLikelihoods[indexOfCount] > minLikelihood && frequenceOfPart[group] > 5){ //>5 because we want at least 5 reads per haplotype
            listOfGroups.emplace(group);
        }
        indexOfCount++;
    }

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

    //now most reads should be assigned to a cluster. Rescue those that have not been assigned or assigned to a rare cluster

    for (auto read = 0 ; read < clusters.size() ; read++){

        if (clusters[read] < 0){
            clusters[read] = -1;
        }

        if (listOfGroups.find(clusters[read]) == listOfGroups.end() && clusters[read] >  -1){ //this means the read needs to be rescued 

            //cout << "rescuing "; for(auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++) {cout << allPartitions[binary][read];} cout << endl;
            //iterate through the list of groups and choose the one that can be explained by the less mistakes
            int bestGroup = -1;
            int bestGroupScore = 0;

            //compute a score for assigning the read to each existing group, then assign the read to the best group
            for (auto group : listOfGroups){

                int score = 0;
                int groupDecomposition = group;
                int readDecomposition = clustersAll[read];

                for (int binary = listOfFinalPartitions.size() -1 ; binary > -1  ; binary--){
                    int expectedAssignation = groupDecomposition%2;
                    int actualAssignation = readDecomposition%32;

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
            //     cout << allMores[binary][read]<< "/" << allLess[binary][read] << " ";
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

    cout << "Rescuing reads\r" << endl;

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









