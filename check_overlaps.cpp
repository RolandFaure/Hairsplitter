#include "check_overlaps.h"
#include "edlib.h"
// #include "spoa/spoa.hpp"
#include "cluster_graph.h"

#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>

#include "robin_hood.h"
#include "input_output.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::min;
using std::max;
using std::begin;
using std::end;
using std::pair;
using std::make_pair;

//definition of a small struct that will be useful later
struct distPart{
    float distance = 1;
    short phased = 1;
};
bool comp (distPart i, distPart j){
    return i.distance < j.distance;
}

//input : the set of all overlaps and the backbone reads
//output : a partition for all backbone reads. All reads also have updated backbone_seqs, i.e. the list of backbone reads they are leaning on
void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, std::vector<vector<short>> &partitions, bool assemble_on_assembly) {

    // vector<char> seq1 = {'A', 'A', 'A', 'A', 'A', 'A', '?', 'A', 'A', 'A', 'C', 'A', 'C', '?', 'A', 'A', 'A', 'C', 'A', 'C', '?'};
    // vector<char> seq2 = {'A', 'A', 'A', 'A', 'A', 'C', 'C', 'A', 'A', 'A', 'A', 'C', 'A', 'C', 'A', 'A', 'A', 'A', 'C', 'A', 'C'};

    // Partition par1(seq1);
    // vector <pair<Partition, int>> partitions = {make_pair(par1, 1)};

    // distance(partitions[0].first, seq2, 0, partitions);
    // partitions[0].first.print();

    //main loop : for each backbone read, build MSA (Multiple Sequence Alignment) and separate the reads
    int index = 0;
    for (unsigned long int read : backbones_reads){
        //524, 67
        if (allreads[read].neighbors_.size() > 0 && index <= 10 /*&& (allreads[read].size() == 13059 || allreads[read].size() == 6786)*/){
            vector<vector<char>> snps (allreads[read].size(), vector<char>(allreads[read].neighbors_.size()+1-int(assemble_on_assembly), '?')); //vector containing list of position, with SNPs at each position
            //first build an MSA
            // cout << "Generating MSA" << endl;
            Partition truePar(0); //for debugging
            float meanDistance = generate_msa(read, allOverlaps, allreads, snps, partitions.size(), truePar, assemble_on_assembly);

            // auto numberOf0s = 0;
            // for (auto i : truePar.getPartition()) {if (i==-1){numberOf0s++;}}
            // if (numberOf0s > 10){
            cout << "Looking at backbone read number " << index << " out of " << backbones_reads.size() << endl;

            //then separate the MSA
            // cout << "Separating reads" << endl;
            auto par = separate_reads(read, allOverlaps, allreads, snps, meanDistance);
            cout << "True partition : " << endl;
            truePar.print();
            // auto par = truePar.getPartition();//DEBUG
            // for (auto i = 0 ; i < par.size() ; i++) {par[i]++; }
            cout << "Proposed partition : " << endl;
            for (auto i = 0 ; i < par.size() ; i++){cout << par[i];}cout << endl;
            cout << endl;

            partitions.push_back(par);
    
        }
        index++;
    }

}

//input: a read with all its neighbor
//outputs : the precise alignment of all reads against input read in the form of matrix snps, return the mean editDistance/lengthOfAlginment
float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, int backboneReadIndex, Partition &truePar, bool assemble_on_assembly){

    //cout << "neighbors of read " << read << " : " << allreads[read].neighbors_.size() << endl;
    //go through the neighbors of the backbone read and align it

    //keep count of the distance between two reads to know the mean distance
    float totalDistance = 0;
    int totalLengthOfAlignment = 0;

     //while all the alignments are computed, build the positions
    robin_hood::unordered_map<int, int> insertionPos;
    vector<int> numberOfInsertionsHere (allreads[read].size(), 0);
    
    //first enter the sequence of the read in snps, if the backbone read is a read, not if it's a contig
    if (!assemble_on_assembly){
        string readSequence = allreads[read].sequence_.str();
        for (auto i = 0 ; i < readSequence.size() ; i++){
            snps[i][snps[i].size()-1] = readSequence[i];
        }
    }

    //to remember on what backbone the backbone read is leaning -> itself
    allreads[read].new_backbone(make_pair(backboneReadIndex, allreads[read].neighbors_.size()), allreads[read].neighbors_.size()+1);


    //small loop to compute truePartition DEBUG
    vector<char> truePartition; //for debugging
    for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){

        long int neighbor = allreads[read].neighbors_[n];
        Overlap overlap = allOverlaps[neighbor];
        if (overlap.sequence1 == read){
            // cout << "name : " << allreads[overlap.sequence2].name << endl;
            if (allreads[overlap.sequence2].name[1] == '1'){
                truePartition.push_back('A');
            }
            else{
                truePartition.push_back('C');
            }
            // cout << "name : " << allreads[overlap.sequence2].name << " " << allreads[overlap.sequence2].name[1] << " " << truePartition[truePartition.size()-1] << endl;
        }
        else{
            if (allreads[overlap.sequence1].name[1] == '1'){
                truePartition.push_back('A');
            }
            else{
                truePartition.push_back('C');
            }
            // cout << "name : " << allreads[overlap.sequence1].name << " " << allreads[overlap.sequence1].name[1] << " " << truePartition[truePartition.size()-1] << endl;
        }    
    }
    //do not forget the last read itself
    if (allreads[read].name[1] == '1'){
        truePartition.push_back('A');
    }
    else{
        truePartition.push_back('C');
    }
    // cout << "name : " << allreads[read].name << " " << allreads[read].name[1] << " " << truePartition[truePartition.size()-1] << endl;
    truePar = Partition(truePartition);

    // for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){
    //     long int neighbor = allreads[read].neighbors_[n];
    //     Overlap overlap = allOverlaps[neighbor];
    //     if (overlap.sequence1 == read){
    //         allreads[overlap.sequence2].new_backbone(make_pair(backboneReadIndex,n), allreads[read].neighbors_.size()+1);
    //     }
    //     else {
    //         allreads[overlap.sequence1].new_backbone(make_pair(backboneReadIndex,n), allreads[read].neighbors_.size()+1);
    //     }
    // }

    // /* compute only true partition
    // cout << "Aligning on backbone of length " << snps.size() << endl;

    for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){

        cout << "Aligned " << n << " reads out of " << allreads[read].neighbors_.size() << " on the backbone\r";

        long int neighbor = allreads[read].neighbors_[n];
        Overlap overlap = allOverlaps[neighbor];

        //cout << "The overlap I'm looking at looks like this: " << overlap.sequence1 << " " << overlap.sequence2 << " " << overlap.strand << endl;

        string toBeAlgined1;
        string toBeAlgined2;
        int true_begin; //true_begin and true_end are there because input alignments are not expected to go to the end of the reads
        int true_end;
        int start1 = 0;
        if (overlap.sequence1 == read){
            if (overlap.strand){ //if the two reads are on the same strand
                //now adjust the positions to go to the end of the reads with true_begin and true_end
                true_begin = min(overlap.position_1_1, overlap.position_2_1);
                true_end = min(allreads[overlap.sequence1].sequence_.size()-overlap.position_1_2, allreads[overlap.sequence2].sequence_.size()-overlap.position_2_2);
                start1 = overlap.position_1_1-true_begin;
                toBeAlgined1 = allreads[overlap.sequence1].sequence_.subseq(start1, overlap.position_1_2-overlap.position_1_1+true_begin+true_end).str();
                toBeAlgined2 = allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1-true_begin, overlap.position_2_2-overlap.position_2_1+true_begin+true_end).str();
                
            }
            else {
                true_begin = min(static_cast<size_t>(overlap.position_1_1), allreads[overlap.sequence2].sequence_.size()-overlap.position_2_2);
                true_end = min(allreads[overlap.sequence1].sequence_.size()-overlap.position_1_2, size_t(overlap.position_2_1));
                start1 = overlap.position_1_1-true_begin;
                toBeAlgined1 = allreads[overlap.sequence1].sequence_.subseq(overlap.position_1_1-true_begin, overlap.position_1_2-overlap.position_1_1+true_begin+true_end).str();
                toBeAlgined2 = allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1-true_end, overlap.position_2_2-overlap.position_2_1+true_end+true_begin).reverse_complement().str();
            }
            true_begin = overlap.position_1_1-true_begin;

            //to remember on what backbone this read is leaning
            allreads[overlap.sequence2].new_backbone(make_pair(backboneReadIndex,n), allreads[read].neighbors_.size()+1);
        }
        else {
            if (overlap.strand){ //if the two reads are on the same strand
                true_begin = min(overlap.position_1_1, overlap.position_2_1);
                true_end = min(allreads[overlap.sequence1].sequence_.size()-overlap.position_1_2, allreads[overlap.sequence2].sequence_.size()-overlap.position_2_2);
                start1 = overlap.position_2_1-true_begin;
                toBeAlgined1 = allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1-true_begin, overlap.position_2_2-overlap.position_2_1+true_begin+true_end).str();
                toBeAlgined2 = allreads[overlap.sequence1].sequence_.subseq(overlap.position_1_1-true_begin, overlap.position_1_2-overlap.position_1_1+true_begin+true_end).str();
            }
            else {
                true_begin = min(size_t(overlap.position_1_1), allreads[overlap.sequence2].sequence_.size()-overlap.position_2_2);
                true_end = min(allreads[overlap.sequence1].sequence_.size()-overlap.position_1_2, size_t(overlap.position_2_1));
                start1 = overlap.position_2_1-true_end;
                toBeAlgined1 = allreads[overlap.sequence2].sequence_.subseq(overlap.position_2_1-true_end, overlap.position_2_2-overlap.position_2_1+true_begin+true_end).str();
                toBeAlgined2 = allreads[overlap.sequence1].sequence_.subseq(overlap.position_1_1-true_begin, overlap.position_1_2-overlap.position_1_1+true_begin+true_end).reverse_complement().str();
            }
            true_begin = overlap.position_2_1-true_begin;

            //to remember on what backbone this read is leaning
            allreads[overlap.sequence1].new_backbone(make_pair(backboneReadIndex,n), allreads[read].neighbors_.size()+1);
        }

        //cout << "Aligning " << toBeAlgined1.c_str() << " and " << toBeAlgined2.c_str() << endl;

        EdlibAlignResult result = edlibAlign(toBeAlgined1.c_str(), toBeAlgined1.size(), toBeAlgined2.c_str(), toBeAlgined2.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
        //cout << "Aligned ! " << result.editDistance << endl;

        totalDistance += result.editDistance;
        totalLengthOfAlignment += result.alignmentLength;

        //cout << overlap.sequence1 << " " << overlap.sequence2 << " " << toBeAlgined1 << " " << toBeAlgined2 << endl;
        //cout << toBeAlgined1 << endl << toBeAlgined2 << endl;
        //cout << "alignment of read " << read << " : " << result.editDistance << " " << result.alignmentLength << ", " << edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD) << endl;  
        //list_of_CIGARs.push_back(edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD));

        // string aln (reinterpret_cast<char*>(result.alignment));

        //a loop going through the CIGAR and modifyning snps
        char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
        int indexQuery = start1; //query corresponds to read
        int indexTarget = 0; //target corresponds to the other read of the overlap
        int numberOfInsertionsThere = 0;
        int nmatches = 0;
        int nmismatches = 0;
        for (int i = 0; i < result.alignmentLength; i++) {
            if (moveCodeToChar[result.alignment[i]] == '=' || moveCodeToChar[result.alignment[i]] == 'X'){
                snps[indexQuery][n] = toBeAlgined2[indexTarget];
                indexQuery++;
                indexTarget++;
                numberOfInsertionsThere = 0;
                if (indexQuery > 900 && indexQuery < 1100){
                    if (moveCodeToChar[result.alignment[i]] == '=' ){
                        nmatches++;
                    }
                    if (moveCodeToChar[result.alignment[i]] == 'X'){
                        nmismatches++;
                    }
                }
            }
            else if (moveCodeToChar[result.alignment[i]] == 'I'){
                snps[indexQuery][n] = '-';
                indexQuery++;
                numberOfInsertionsThere = 0;
            }
            else if (moveCodeToChar[result.alignment[i]] == 'D'){ //hardest one
                if (numberOfInsertionsHere[indexQuery] < 99) {

                    if (numberOfInsertionsThere >= numberOfInsertionsHere[indexQuery]) { //i.e. this is a new column
                        insertionPos[100*indexQuery+numberOfInsertionsHere[indexQuery]] = snps.size();
                        numberOfInsertionsHere[indexQuery] += 1;
                        vector<char> newInsertedPos(snps[0].size(), '-');
                        for (auto re = 0 ; re < newInsertedPos.size() ; re++){
                            if (snps[indexQuery][re] == '?'){
                                newInsertedPos[re] = '?';
                            }
                        }
                        snps.push_back(newInsertedPos);
                        snps[snps.size()-1][n] = toBeAlgined2[indexTarget];
                    }
                    else{
                        snps[insertionPos[100*indexQuery+numberOfInsertionsThere]][n] = toBeAlgined2[indexTarget];
                    }
                    numberOfInsertionsThere ++;
                }
                indexTarget++;
            }
            //cout << "out" << endl;
        }

        // int nbread = 0;
        // if (n == nbread){
        //     cout << "Added the alignment to snps : " << float(result.editDistance)/result.alignmentLength << endl;
        //     cout << nmatches << " " << nmismatches << endl;

        //     for(auto i = 900 ; i< min(snps.size(), size_t(1100)); i++){
        //         cout << snps[i][nbread];
        //     }
        //     cout << " " << nbread << endl;

        //     for(auto i = 900 ; i< min(snps.size(), size_t(1100)); i++){
        //         cout << snps[i][numberOfNeighbors];
        //     }
        //     cout << " " << nbread << endl << endl;
        // }

        edlibFreeAlignResult(result);

    }

    // cout << "Finished aligning                                                                      " << endl;


    //print snps (just for debugging)
    // int step = 2;
    // int numberOfReads = 10;
    // int start = 10000;
    // int end = 10400;
    // vector<string> reads (min(int(snps[0].size()), numberOfReads));
    // for (unsigned short i = start ; i < end; i++){
        
    //     for (short n = 0 ; n < min(int(snps[0].size()), numberOfReads*step) ; n+= step){
    //         reads[n/step] += snps[i][n];
    //     }
    //     for (short insert = 0 ; insert < min(99,numberOfInsertionsHere[i]) ; insert++ ){
    //         auto snpidx = insertionPos[100*i+insert];
    //         for (short n = 0 ; n < min(int(snps[0].size()), numberOfReads*step) ; n+= step){
    //             reads[n/step] += snps[snpidx][n];
    //         }
    //     }
    // }
    // cout << "Here are the aligned reads : " << endl;
    // for (auto neighbor : reads){
    //     cout << neighbor << endl;
    // }
    // for (unsigned short i = start ; i < end; i++){cout << snps[i][1];} cout << endl;
    // for (unsigned short i = start ; i < end; i++){cout << snps[i][snps[0].size()-1];} cout << endl;

    // cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;
    return totalDistance/totalLengthOfAlignment;
    ///
}

//input : a set of reads aligned to read in matrix snps
//output : reads separated by their region of origin
vector<short> separate_reads(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, float meanDistance){

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

    for (int position = 0 ; position < snps.size() ; position++){

        //first look at the position to see if it is suspect
        int content [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for -
        int numberOfReads = 0;
        for (short n = 0 ; n < snps[position].size() ; n++){
                char base = snps[position][n];
                if (base != '?' && bases2content.contains(base)){
                    content[bases2content[base]] += 1;
                    numberOfReads += 1;
                }
        }

        float threshold = 1 + numberOfReads*meanDistance/2 + 3*sqrt(numberOfReads*meanDistance/2*(1-meanDistance/2));
        //DEBUG
        //threshold = 3;
        if (*std::max_element(content, content+5) < numberOfReads-threshold){ //this position is suspect
            //cout << threshold << " " << position << " ;bases : " << content[0] << " " << content[1] << " " << content[2] << " " << content[3] << " " << content[4] << endl;
            suspectPostitions.push_back(position);
            //go through the partitions to see if this suspicious position looks like smt we've seen before
            vector<distPart> distances (distanceBetweenPartitions.size());
            bool found = false;
            for (auto p = 0 ; p < partitions.size() ; p++){

                distancePartition dis = distance(partitions[p], snps[position]);
                if (float(dis.n01+dis.n10)/(dis.n00+dis.n11+dis.n01+dis.n10) <= meanDistance && dis.augmented ){
                    // cout << "merging..." << endl;;
                    // partitions[p].print();
                    // for (auto i : dis.partition_to_augment) {cout << i;} cout << endl;
                    partitions[p].augmentPartition(dis.partition_to_augment);
                    found = true;
                    break;
                }
                auto comparable = min(dis.n00,dis.n11) + dis.n01 + dis.n10;
                if (comparable > 5){
                    //distances[p].distance = max(float(0), (10-dis.chisquare)/10) ;
                    distances[p].distance = float(dis.n01+dis.n10)/(dis.n00 + dis.n11 + dis.n01 + dis.n10);
                    distances[p].phased = dis.phased;
                }
                else{
                    distances[p].distance = 1;
                }
                
            }

            if (!found){
                partitions.push_back(Partition(snps[position]));        
                distanceBetweenPartitions.push_back(distances);
            }
            
            numberOfSuspectPostion += 1;

            //two suspect positions next to each other can be artificially correlated through alignement artefacts
            position += 5;
        }
    }


    if (partitions.size() == 0){ //there are no position of interest
        return vector<short> (snps[0].size(), 1);
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

    float threshold = max(4.0, min(0.1*numberOfSuspectPostion, 0.001*snps.size()));

    // cout << "Here are all the partitions : " << endl;
    // for (auto p : partitions){
    //     p.print();
    // }
    // cout << "threshold : " << threshold << endl;

    vector<Partition> listOfFinalPartitions;
    for (auto p1 = 0 ; p1 < partitions.size() ; p1++){
        
        if (partitions[p1].number() > threshold && partitions[p1].isInformative(meanDistance/2, true)){

            // cout << "informative partition 2 : " << endl;
            // partitions[p1].print();
            bool different = true;
            for (auto p2 = 0 ; p2 < listOfFinalPartitions.size() ; p2++){

                bool same = distance(listOfFinalPartitions[p2], partitions[p1], 9, 2);
                if (same){ //no need to merge the two partitions, they have already been merged by distance()
                    // adjMatrix[p1][p2] = dis.chisquare-9;
                    // adjMatrix[p2][p1] = dis.chisquare-9;
                    // cout << "those two partitions are the same : " << endl;
                    // partitions[p1].print();
                    // listOfFinalPartitions[p2].print();
                    different = false;
                }    
                else {
                    // cout << "those two partitions are not the same : " << endl;
                    // partitions[p1].print();
                    // listOfFinalPartitions[p2].print();
                }
            }
            if (different){
                listOfFinalPartitions.push_back(partitions[p1]);
            }
        }
    }

    //now we have the list of final partitions : there may be several, especially if there are more than two copies

    if (listOfFinalPartitions.size() == 0){
        Partition p(snps[0].size());
        listOfFinalPartitions.push_back(p);
    }

    cout << "I end up with " << listOfFinalPartitions.size() << " partitions, deduced from " << partitions.size() << endl; 
    cout  << "I have " << numberOfSuspectPostion << " suspect positions, and " << listOfFinalPartitions[0].getPartition().size() << " reads " << endl;
    
    //now aggregate all those binary partitions in one final partition. There could be up to 2^numberBinaryPartitions final groups
    vector<short> finalClusters = threadHaplotypes(listOfFinalPartitions);


    return finalClusters;

}

//input : one partition and one list of chars
//output : is the list of chars close to the partition ? If so, augment the partition. Return the chi-square of the difference in bases
distancePartition distance(Partition &par1, vector<char> &par2){

    /*
    when computing the distance, there is not 5 letters but 2 : the two alleles, which are the two most frequent letters
    */
    distancePartition res;
    res.nonComparable = 0;
    res.augmented = true;
    vector<short> part1 = par1.getPartition();
    vector<float> confs1 = par1.getConfidence();

    robin_hood::unordered_flat_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['-'] = 4;
    bases2content['?'] = 5;
    
    int content2 [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for *, 5 for '?'
    int maxFrequence = 0;
    float numberOfBases = 0;
    for (int c = 0 ; c < part1.size() ; c++){
        if (part1[c] != 0 && par2[c] != '?'){
            if (part1[c] == 1){
                maxFrequence += 1;
            }
            numberOfBases++;
        }
        content2[bases2content[par2[c]]] += 1;
    }

    if (numberOfBases < 10){ //not comparable
        res.n00 = 0;
        res.n01 = 0;
        res.n10 = 0;
        res.n11 = 0;
        res.augmented = false;
        res.nonComparable = par2.size();
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
    vector<short> newPartitions [2];

    for (auto c = 0 ; c < par2.size() ; c++){

        float conf = 2*confs1[c]-1;

        if (par2[c] == mostFrequent2){
            newPartitions[0].push_back(1);
            newPartitions[1].push_back(-1);

            if (part1[c] == 1){
                scores[0] += conf;
                scores[1] -= conf;
                bestScore += conf;
                matches11[0] += 1;
                matches10[1] += 1;
            }
            else if (part1[c] == -1){
                scores[0] -= conf;
                scores[1] += conf;
                bestScore += conf;
                matches01[0] += 1;
                matches00[1] += 1;
            }
            else {
                res.nonComparable += 1;
            }
        }
        else if (par2[c] == secondFrequent2){
            newPartitions[0].push_back(-1);
            newPartitions[1].push_back(1);

            if (part1[c] == 1){
                scores[0] -= conf;
                scores[1] += conf;
                bestScore += conf;
                matches10[0] += 1;
                matches11[1] += 1;
            }
            else if (part1[c] == -1){
                scores[0] += conf;
                scores[1] -= conf;
                bestScore += conf;
                matches00[0] += 1;
                matches01[1] += 1;
            }
            else {
                res.nonComparable += 1;
            }
        }
        else{ //if this is a '?' or an error
            newPartitions[0].push_back(0);
            newPartitions[1].push_back(0);
            res.nonComparable += 1;
        }
        // if (c == 22){
        //     cout << "positions at 22 : " << par2[c] << " " << mostFrequent2 << " " << secondFrequent2 << " " << part1[c] << endl;
        // }

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

float computeChiSquare(distancePartition dis){

    int n = dis.n00 + dis.n01 + dis.n10 + dis.n11;
    if (n == 0){
        return 0;
    }
    float pmax1 = float(dis.n10+dis.n11)/n;
    float pmax2 = float(dis.n01+dis.n11)/n;
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
bool distance(Partition &par1, Partition &par2, float thresholdChi, int threshold_p){
    /*
    Two metrics are used to compare two partition : the chi to see if the two partition correlate when they are small
                                                    the p when the partitions are bigger and you can do probabilities on them
    */

    float chi = 0;

    vector<short> part1 = par1.getPartition();
    vector<short> part2 = par2.getPartition();

    vector<int> more1 = par1.getMore();
    vector<int> less1 = par1.getLess();

    vector<int> more2 = par2.getMore();
    vector<int> less2 = par2.getLess();
    
    int maxFrequence = 0;
    int maxFrequence2 = 0;
    float numberOfBases = 0;
    for (int c = 0 ; c < part1.size() ; c++){
        if (part1[c] != 0 && part2[c] != 0){
            if (part1[c] == 1){
                maxFrequence += 1;
            }
            if (part2[c] == 1){
                maxFrequence2 += 1;
            }
            numberOfBases++;
        }
    }

    if (numberOfBases < 10){ //not comparable
        return false;
    }

    int scores [2] = {0,0}; //the scores when directing mostFrequent on either mostfrequent2 or secondFrequent2
    short ndivergentPositions[2] = {0,0}; //number of positions where these two partitions could not have been so different by chance
    //remember all types of matches for the chi square test
    int matches00[2] = {0,0};
    int matches01[2] = {0,0};
    int matches10[2] = {0,0};
    int matches11[2] = {0,0};

    for (auto c = 0 ; c < part2.size() ; c++){

        float threshold1 = 0.5*(more1[c]+less1[c]) + 3*sqrt((more1[c]+less1[c])*0.5*(1-0.5)); //to check if we deviate significantly from the "random read", that is half of the time in each partition
        float threshold2 = 0.5*(more2[c]+less2[c]) + 3*sqrt((more2[c]+less2[c])*0.5*(1-0.5)); //to check if we deviate significantly from the "random read", that is half of the time in each partition

        if (part2[c] == 1){

            if (part1[c] == 1){
                scores[0] += 1;
                scores[1] -= 1;
                matches11[0] += 1;
                matches10[1] += 1;

                //if both positions are certain, this may be bad
                if (more1[c] > threshold1 && more2[c] > threshold2){
                    ndivergentPositions[1] += 1;
                }
            }
            else if (part1[c] == -1){
                scores[0] -= 1;
                scores[1] += 1;
                matches01[0] += 1;
                matches00[1] += 1;

                //if both positions are certain, this may be bad
                if (more1[c] > threshold1 && more2[c] > threshold2){
                    ndivergentPositions[0] += 1;
                }
            }
        }
        else if (part2[c] == -1){

            if (part1[c] == 1){
                scores[0] -= 1;
                scores[1] += 1;
                matches10[0] += 1;
                matches11[1] += 1;

                //if both positions are certain, this may be bad
                if (more1[c] > threshold1 && more2[c] > threshold2){
                    ndivergentPositions[0] += 1;
                }
            }
            else if (part1[c] == -1){
                scores[0] += 1;
                scores[1] -= 1;
                matches00[0] += 1;
                matches01[1] += 1;

                //if both positions are certain, this may be bad
                if (more1[c] > threshold1 && more2[c] > threshold2){
                    ndivergentPositions[1] += 1;
                }
            }
        }

    }

    //check if there are too many unenxplainable positions
    if (ndivergentPositions[0] >= threshold_p && ndivergentPositions[1] >= threshold_p){
        // cout << "Should not merge those two partitions ! " << endl;
        return false;
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

    int occurences[2] = {maxFrequence2, numberOfBases-int(maxFrequence2)};
 
    //now compute the chi square
    int n = numberOfBases;
    float pmax1 = maxFrequence/numberOfBases;
    float pmax2 = float(occurences[maxScoreIdx]) / numberOfBases;
    if (pmax1*pmax2*(1-pmax1)*(1-pmax2) == 0){ //if there is i.e. only one base in one partition, it can't be compared
        chi = 0;
    }
    //chi square test with 1 degree of freedom
    chi = pow((matches00[maxScoreIdx]-(1-pmax1)*(1-pmax2)*n),2)/((1-pmax1)*(1-pmax2)*n)
                        + pow((matches01[maxScoreIdx]-(1-pmax1)*pmax2*n),2)/((1-pmax1)*pmax2*n)
                        + pow((matches10[maxScoreIdx]-pmax1*(1-pmax2)*n),2)/(pmax1*(1-pmax2)*n)
                        + pow((matches11[maxScoreIdx]-pmax1*pmax2*n),2)/(pmax1*pmax2*n);

    //cout << "going chisquare : " << chi << endl;
    bool same = (chi > thresholdChi);

    if (same){
        // cout << "Let's merge, baby !:" << endl;
        // par1.print();
        // par2.print();
        par1.mergePartition(par2, -maxScoreIdx*2+1);
    }
    else {
        //cout << "Should not merge those two partitions ? " << chi << endl;
    }

    return same ;
}

//input : a list of all binary partitions found in the reads
//output : a single partition, with several number corresponding to several clusters
vector<short> threadHaplotypes(vector<Partition> &listOfFinalPartitions){

    vector<int> frequenceOfPart (pow(2, listOfFinalPartitions.size()));

    vector<short> res (listOfFinalPartitions[0].size(), -1);

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
    vector<vector<int>> allLess;
    for (auto i : listOfFinalPartitions){
        allLess.push_back(i.getLess());
    }

    //as a first step, we'll try to throw away the partitions that do not look very sure of themselves
    /*to to that, we'll look at each read, see if it's well clustered in one partition. 
    If yes, all partitions where it's badly clustered are suspicious
    Indeed, we expect the different errors to compensate each other
    */
    vector<Partition> trimmedListOfFinalPartition;

    vector<bool> readsClassified (listOfFinalPartitions[0].size(), false);
    for (int par = 0 ; par < listOfFinalPartitions.size() ; par++){
        for (int read = 0 ; read < listOfFinalPartitions[0].size() ; read++){
            if (allPartitions[par][read] != 0 && allConfidences[par][read] >= 0.7){
                readsClassified[read] = true;
            }
        }
    }

    //now we know what reads are not well classified anywhere. Now see if some partitions are just not sure enough of themselves
    for (int par = 0 ; par < listOfFinalPartitions.size() ; par++){

        //because of the ref, the two haplotypes are not strictly equivalent : try to compensate
        int numberOfPartitionnedRead = 0;
        float means [2] = {0,0};
        float n1 = 0;
        float n0 = 0;
        for (int read = 0 ; read < listOfFinalPartitions[0].size() ; read++){

            if (allPartitions[par][read] != 0){
                numberOfPartitionnedRead++;
                if (allPartitions[par][read] == 1){
                    means[1] += allConfidences[par][read];
                    n1++;
                }
                else if (allPartitions[par][read] == -1){
                    means[0] += allConfidences[par][read];
                    n0++;
                }
            }
        }

        double errors [2] = {(1-means[0]/n0)/(2-means[0]/n0-means[1]/n1) * 0.6 , (1-means[1]/n1)/(2-means[0]/n0-means[1]/n1) * 0.6 };
        //cout << "the center is : " << means[0]/n0 << " " << means[1]/n1 << endl;

        int numberOfUnsureReads = 0;
        for (int read = 0 ; read < listOfFinalPartitions[0].size() ; read++){

            if (allPartitions[par][read] != 0 && allConfidences[par][read] < 1-errors[(allPartitions[par][read]+1)/2] && readsClassified[read]){ //oops
                numberOfUnsureReads++;
            }
        }

        //cout << "Here is the number of unsure reads " << numberOfUnsureReads << " "<< numberOfPartitionnedRead << endl;

        if (float(numberOfUnsureReads)/numberOfPartitionnedRead < 0.1){
            trimmedListOfFinalPartition.push_back(listOfFinalPartitions[par]);
        }

    }
    listOfFinalPartitions = trimmedListOfFinalPartition;

    allPartitions = {};
    for (auto i : listOfFinalPartitions){
        allPartitions.push_back(i.getPartition());
    }
    allConfidences = {};
    for (auto i : listOfFinalPartitions){
        allConfidences.push_back(i.getConfidence());
    }

    // for (auto p : listOfFinalPartitions){
    //     cout << "informative partition : ";
    //     p.print();
    //     double conf = 1;
    //     int numberReads = 0;
    //     auto mores = p.getMore();
    //     auto confidences = p.getConfidence();
    //     for (auto c = 0 ; c < confidences.size() ; c++) {
    //         if (mores[c] > 1){
    //             conf*=confidences[c];
    //             numberReads++;
    //         }   

    //     }
    //     cout << numberReads << " " << exp(log(conf)/numberReads)  << endl;
    // }

    //as a second step, we'll thread the haplotypes between the partitions

    robin_hood::unordered_flat_map <int, int> count; //a map counting how many times a cluster appears
    int numberOfAssignedReads=0;

    for (int pos=0 ; pos < listOfFinalPartitions[0].size() ; pos++){
        int id = 0;
        for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
            id *= 2;
            auto camp = allPartitions[binary][pos];
            if (camp == 0 || allConfidences[binary][pos] < 0.7){ //0.7 to be pretty confident about the reads we separate we have (more than 70% on all partitions)
                id = -1;
                break;
            }
            else if (camp == 1 ){
                id += 1;
            }
        }
        if (id>=0){
            frequenceOfPart[id] += 1;
            res[pos] = id;
            count[id] += 1;
            numberOfAssignedReads ++;
        }
    }

    
    //to make sure we don't over-estimate the number of clusters, we'll make the assumption that all haplotypes have haplotype-specific mutation => the number of final cluster cannot be higher that the number of partitions
    vector<float> proportionOf1;
    for (auto p : listOfFinalPartitions){
        proportionOf1.push_back(p.proportionOf1());  
    }
    vector<float> listOfLikelihoods;
    vector<int> listOfAbundances;
    for (auto group : count){
        int id = group.first;
        double proba = 1;
        for (int binary = listOfFinalPartitions.size()-1 ; binary > -1 ; binary--){
            if (id%2 == 1){
                proba *= proportionOf1[binary];
            }
            else{
                proba *= 1-proportionOf1[binary];
            }
            id = int(id / 2);
        }
        // cout << "group " << group.first << " is expected " << numberOfAssignedReads*proba << " times and arrives " << group.second << " times " << endl;
        listOfLikelihoods.push_back(float(group.second -numberOfAssignedReads*proba) /group.second);
        listOfAbundances.push_back(group.second);  
    }
    int minLikelihood = -10;
    if (listOfFinalPartitions.size() != 1){
        vector <int> minElements (listOfAbundances.size()-listOfFinalPartitions.size());
        std::partial_sort_copy(listOfLikelihoods.begin(),  listOfLikelihoods.end(), minElements.begin(), minElements.end());
        minLikelihood = minElements[minElements.size()-1];
    }

    // cout << "minNumber of reads : " << minNumberOfReads << endl;
    // for (auto a : listOfAbundances) {cout << a << ",";}cout << endl;

    //establish a list of possible clusters
    std::set <int> listOfGroups;
    int indexOfCount = 0;
    for (auto group : count){
        if (listOfLikelihoods[indexOfCount] >= minLikelihood){
            listOfGroups.emplace(group.first);
        }
        indexOfCount++;
    }

    // cout << "possible clusters : "; for (auto g : listOfGroups){cout << g << " ";} cout << endl;
    // cout << "Res : ";
    // for (auto i : res){
    //     cout << i << ",";
    // }
    // cout << endl;

    //now most reads should be assigned to a cluster. Rescue those that have not been assigned or assigned to a singleton cluster

    for (auto read = 0 ; read < res.size() ; read++){

        if (listOfGroups.find(res[read]) == listOfGroups.end()){ //this means the read needs to be rescued 

            bool atLeastOnePartitionContainThisRead = false;
            for (auto p : listOfFinalPartitions){
                if (p.getPartition()[read] != 0){
                    atLeastOnePartitionContainThisRead = true;
                }
            }

            if (atLeastOnePartitionContainThisRead){
                //cout << "rescuing "; for(auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++) {cout << allPartitions[binary][read];} cout << endl;
                //iterate through the list of groups and choose the one that can be explained by the less mistakes
                int bestGroup = -1;
                int bestGroupScore = 1000000;

                //compute a score for each existing group, then assign the read to that group
                for (auto group : listOfGroups){

                    int score = 0;
                    int groupDecomposition = group;

                    for (int binary = listOfFinalPartitions.size() -1 ; binary > -1  ; binary--){
                        int expectedAssignation = groupDecomposition%2;

                        if ((allPartitions[binary][read]+1)/2 != expectedAssignation && allPartitions[binary][read] != 0){
                            score += allMores[binary][read]-allLess[binary][read];
                        }

                        groupDecomposition /= 2;
                    }
                    if (score < bestGroupScore){
                        bestGroup = group;
                        bestGroupScore = score;
                    }
                }
                res[read] = bestGroup;

                // cout << "Rescuing ";
                // for (int binary = 0 ; binary < listOfFinalPartitions.size()  ; binary++){
                //     cout << allMores[binary][read]<< "/" << allLess[binary][read] << " ";
                // }
                // cout << " as : " << bestGroup << endl;
            }
        }
    }

    // cout << "Res : ";
    // for (auto i : res){
    //     cout << i << ",";
    // }
    // cout << endl;

    return res;

}