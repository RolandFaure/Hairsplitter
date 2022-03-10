#include "check_overlaps.h"
#include "edlib.h"
// #include "spoa/spoa.hpp"
#include "cluster_graph.h"

#include <cmath>
#include <algorithm>
#include <unordered_map>

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
void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, std::vector<vector<short>> &partitions) {

    // vector<char> seq1 = {'A', 'A', 'A', 'A', 'A', 'A', '?', 'A', 'A', 'A', 'C', 'A', 'C', '?', 'A', 'A', 'A', 'C', 'A', 'C', '?'};
    // vector<char> seq2 = {'A', 'A', 'A', 'A', 'A', 'C', 'C', 'A', 'A', 'A', 'A', 'C', 'A', 'C', 'A', 'A', 'A', 'A', 'C', 'A', 'C'};

    // Partition par1(seq1);
    // vector <pair<Partition, int>> partitions = {make_pair(par1, 1)};

    // distance(partitions[0].first, seq2, 0, partitions);
    // partitions[0].first.print();

    //main loop : for each backbone read, build MSA (Multiple Sequence Alignment) and separate the reads
    int index = 0;
    for (unsigned long int read : backbones_reads){
        if (allreads[read].neighbors_.size() > 0 /*&& index ==2*/){
            cout << "Looking at backbone read number " << index << " out of " << backbones_reads.size() << endl;
            vector<vector<char>> snps (allreads[read].size(), vector<char>(allreads[read].neighbors_.size()+1, '?')); //vector containing list of position, with SNPs at each position
            //first build an MSA
            cout << "Generating MSA" << endl;
            Partition truePar(0); //for debugging
            float meanDistance = generate_msa(read, allOverlaps, allreads, snps, partitions.size(), truePar);
            // for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){
            //     if (partitions[allreads[i].backbone_seq[b].first].getPartition().size() <= allreads[i].backbone_seq[b].second){
            //         throw std::logic_error("BIBBU");
            //     }
            // }

            //then separate the MSA
            cout << "Separating reads" << endl;
            auto par = separate_reads(read, allOverlaps, allreads, snps, meanDistance);
            // cout << "True partition : " << endl;
            // truePar.print();
            partitions.push_back(par);
        }
        index++;
    }

}

//input: a read with all its neighbor
//outputs : the precise alignment of all reads against input read in the form of matrix snps, return the mean editDistance/lengthOfAlginment
float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, int backboneReadIndex, Partition &truePar){

    vector<char> truePartition; //for debugging
    //cout << "neighbors of read " << read << " : " << allreads[read].neighbors_.size() << endl;
    //go through the neighbors of the backbone read and align it
    unsigned short numberOfNeighbors = allreads[read].neighbors_.size();

    //keep count of the distance between two reads to know the mean distance
    float totalDistance = 0;
    int totalLengthOfAlignment = 0;

     //while all the alignments are computed, build the positions
    robin_hood::unordered_map<int, int> insertionPos;
    vector<int> numberOfInsertionsHere (allreads[read].size(), 0);
    
    //first enter the sequence of the read in snps
    string readSequence = allreads[read].sequence_.str();
    for (auto i = 0 ; i < readSequence.size() ; i++){
        snps[i][numberOfNeighbors] = readSequence[i];
    }

    //to remember on what backbone the backbone read is leaning -> itself
    allreads[read].new_backbone(make_pair(backboneReadIndex, allreads[read].neighbors_.size()), allreads[read].neighbors_.size()+1);

    for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){
        long int neighbor = allreads[read].neighbors_[n];
        Overlap overlap = allOverlaps[neighbor];
        if (overlap.sequence1 == read){
            //cout << "name : " << allreads[overlap.sequence2].name << endl;
            if (allreads[overlap.sequence2].name[1] == '0'){
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
            // cout << "name : " << allreads[overlap.sequence1].name << " " << allreads[overlap.sequence1].name[1] << " " << truePartition[truePartition.size()-1] << endl
        }    
    }
    truePar = Partition(truePartition);

    for (auto n = 0 ; n<allreads[read].neighbors_.size() ; n++){

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
                        snps.push_back(vector<char> (numberOfNeighbors+1, '-'));
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

    //cout << "Finished aligning" << endl;

    //get rid of the * that are on inserted columns but between two -
    for (size_t n = 0 ; n < numberOfNeighbors+1 ; n++){

        for (auto i = 0 ; i < numberOfInsertionsHere.size(); i++){
            //cout << "two : " << i << " " << snps.size() << endl; 
            if (snps[i][n] == '?'){
                //cout << "one : " << i << " "<< numberOfInsertionsHere.size() << " " << snps.size() << endl;
                for (short insert = 0 ; insert < min(99,numberOfInsertionsHere[i]) ; insert++){
                    snps[insertionPos[100*i+insert]][n] = '?';
                }
            }
        }
    }

    //print snps (just for debugging)
    // int step = 2;
    // int numberOfReads = 10;
    // int start = 300;
    // int end = 400;
    // vector<string> reads (min(numberOfNeighbors+1, numberOfReads));
    // for (unsigned short i = start ; i < end; i++){
        
    //     for (short n = 0 ; n < min(numberOfNeighbors+1, numberOfReads*step) ; n+= step){
    //         reads[n/step] += snps[i][n];
    //     }
    //     for (short insert = 0 ; insert < min(99,numberOfInsertionsHere[i]) ; insert++ ){
    //         auto snpidx = insertionPos[100*i+insert];
    //         for (short n = 0 ; n < min(numberOfNeighbors+1, numberOfReads*step) ; n+= step){
    //             reads[n/step] += snps[snpidx][n];
    //         }
    //     }
    // }
    // cout << "Here are the aligned reads : " << endl;
    // for (auto neighbor : reads){
    //     cout << neighbor << endl;
    // }
    // for (unsigned short i = start ; i < end; i++){cout << snps[i][1];} cout << endl;
    // for (unsigned short i = start ; i < end; i++){cout << snps[i][numberOfNeighbors];} cout << endl;

    // cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;
    return totalDistance/totalLengthOfAlignment;
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
        //threshold = 20;
        if (*std::max_element(content, content+5) < numberOfReads-threshold){ //this position is suspect
            //cout << threshold << " " << position << " ;bases : " << content[0] << " " << content[1] << " " << content[2] << " " << content[3] << " " << content[4] << endl;
            //go through the partitions to see if this suspicious position looks like smt we've seen before
            vector<distPart> distances (distanceBetweenPartitions.size());
            bool found = false;
            for (auto p = 0 ; p < partitions.size() ; p++){

                distancePartition dis = distance(partitions[p], snps[position], meanDistance/2);
                auto comparable = dis.n00 + dis.n01 + dis.n10; //do not compare 11, because it is too easy to match
                if (comparable > 5){
                    //distances[p].distance = float(dis.n01+dis.n10)/comparable;
                    distances[p].distance = float(dis.n01+dis.n10)/(comparable+dis.n11);
                    distances[p].phased = dis.phased;
                }
                else{
                    distances[p].distance = 1;
                }
                if (float(dis.nmismatch)/(dis.nmismatch+dis.nmatch) <= meanDistance ) {//same test as in distance
                    found = true;
                    break;
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

    cout << "Outputting the graph" << endl;

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
    // vector<vector <float>> adj (distanceBetweenPartitions.size(), vector<float> (distanceBetweenPartitions.size(), 0));
    // //each node keeps only the links to very close elements
    // int numberOfNeighborsKept = 5;
    // for (auto i = 0 ; i < distanceBetweenPartitions.size() ; i++){
    //     if (distanceBetweenPartitions[i].size() > numberOfNeighborsKept){
    //         vector <distPart> minElements (numberOfNeighborsKept);
    //         std::partial_sort_copy(distanceBetweenPartitions[i].begin(), distanceBetweenPartitions[i].end(), minElements.begin(), minElements.end(), comp);
    //         float maxDiff = minElements[numberOfNeighborsKept-1].distance;

    //         //cout << "maxdiff of partition " << i << " : " << maxDiff << endl;
    //         if (i == 1){
    //             // cout << "1 : " << endl;
    //             // for(auto j:adj[i]){cout<<j << " ";}
    //             // cout << endl;
    //         }
            
    //         int numberOf1s = 0;
    //         for (int j = 0 ; j < distanceBetweenPartitions[i].size() ; j++){
    //             if (distanceBetweenPartitions[i][j].distance > maxDiff || distanceBetweenPartitions[i][j].distance == 1 || numberOf1s >= numberOfNeighborsKept){
                    
    //             }
    //             else {
    //                 if (distanceBetweenPartitions[i][j].distance < meanDistance*2){
    //                     // cout << "distance : " << distanceBetweenPartitions[i][j].distance << endl;
    //                     adj[i][j] = 1;
    //                     adj[i][j] = 1;
    //                     numberOf1s++;
    //                 }
    //             }
    //         }
    //     }
    // }

    // for (auto i = 0 ; i < adj.size() ; i++){
    //     auto count = 0;
    //     for (auto j = 0 ; j < adj.size() ; j++)
    //     {
    //         if (adj[j][i] == 1){
    //             count += 1;
    //         }
    //     }
    //     if (count > 100){
    //         cout << "here is the wonderful partition : " << endl;
    //         partitions[i].print();
    //     }
    // }

    /*now cluster the different partitions, to achieve ( 1 partition = 1 actual division of the reads)
     -> start by chisquaring two partitions to see if they correlate
     -> then test if it is possible that "they are so different by chance"
    */

    // vector<int> clusters (partitions.size());

    // cluster_graph_chinese_whispers(adj, clusters);
    // outputGraph(adj, clusters, "graph.gdf");

    //filter out too small clusters
    // vector <int> sizeOfCluster;
    // for (auto p = 0 ; p < partitions.size() ; p++){
    //     while (clusters[p] >= sizeOfCluster.size() ){ 
    //         sizeOfCluster.push_back(0);
    //     }
    //     sizeOfCluster[clusters[p]] += 1;
    // }
    // vector<Partition> listOfFinalPartitions(sizeOfCluster.size(), Partition(partitions[0].size()));

    // //now merge each cluster in one partition
    // for (auto p = 0 ; p < partitions.size() ; p++){

    //     if (sizeOfCluster[clusters[p]] > 3){
    //         // if (sizeOfCluster[clusters[p]] == 41){
    //         //     partitions[p].print();
    //         // }
    //         listOfFinalPartitions[clusters[p]].mergePartition(partitions[p]);
    //     }
    // }

    // //filter out empty partitions
    // vector<Partition> listOfFinalPartitions2;
    // for (auto p : listOfFinalPartitions){
    //     if (p.number()> 0 && p.isInformative(meanDistance/2, true)) {
    //         listOfFinalPartitions2.push_back(p);
    //     }
    // }
    // listOfFinalPartitions = listOfFinalPartitions2;

    float threshold = min(0.1*numberOfSuspectPostion, 0.001*snps.size());

    // cout << "Here are all the partitions : " << endl;
    // for (auto p : partitions){
    //     p.print();
    // }
    // cout << "threshold : " << threshold << endl;

    vector<Partition> listOfFinalPartitions;
    for (auto p1 = 0 ; p1 < partitions.size() ; p1++){
        
        if (partitions[p1].number() > threshold && partitions[p1].isInformative(meanDistance/2, true)){

            // cout << "informative partition : ";
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

    //now we have the list of final partitions : there may be several if there are more than two copies

    if (listOfFinalPartitions.size() == 0){
        Partition p(snps[0].size());
        listOfFinalPartitions.push_back(p);
    }

    cout << "I end up with " << listOfFinalPartitions.size() << " partitions, deduced from " << partitions.size() << ". On average, one snps contribute to " << float(numberOfNeighbors)/numberOfSuspectPostion << " partitions" << endl; 
    cout << "Length of backbone read " << snps.size() << ". I have " << numberOfSuspectPostion << " suspect positions, and " << listOfFinalPartitions[0].getPartition().size() << " reads " << endl;
    for (auto p : listOfFinalPartitions){
        // if (p.number() == 33 || p.number() == 41){
            p.print();
        // }
        // cout << p.getMore()[0] << "/" << p.getLess()[0] << " " << p.getMore()[1] << "/" << p.getLess()[1] << endl;
    }


    //now aggregate all those binary partitions in one final partition. There could be up to 2^numberBinaryPartitions final groups
    vector<short> finalClusters = threadHaplotypes(listOfFinalPartitions);


    return finalClusters;

}

//input : one partition and one list of chars
//output : is the list of chars close to the partition ? If so, augment the partition. Return the chi-square of the difference in bases
distancePartition distance(Partition &par1, vector<char> &par2, float errorRate){

    /*
    when computing the distance, there is not 5 letters but 2 : the two alleles, which are the two most frequent letters
    */
    distancePartition res;
    res.nonComparable = 0;
    vector<short> part1 = par1.getPartition();

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
            mostFrequent2 = "ACGT*"[i];
        }
        else if (content2[i] > secondFrequence2){
            secondFrequent2 = "ACGT*"[i];
            secondFrequence2 = content2[i];
        }
    }
    // cout << "Two most frequent : " << mostFrequent2 << "," << secondFrequent2 << " : ";
    // for (auto i : par2){cout << i;} cout << endl;

    float scores [2] = {0,0}; //the scores when directing mostFrequent on either mostfrequent2 or secondFrequent2
    //remember all types of matches for the chi square test
    int matches00[2] = {0,0};
    int matches01[2] = {0,0};
    int matches10[2] = {0,0};
    int matches11[2] = {0,0};
    vector<short> newPartitions [2];

    for (auto c = 0 ; c < par2.size() ; c++){

        if (par2[c] == mostFrequent2){
            newPartitions[0].push_back(1);
            newPartitions[1].push_back(-1);

            if (part1[c] == 1){
                scores[0] += 1;
                scores[1] -= 1;
                matches11[0] += 1;
                matches10[1] += 1;
            }
            else if (part1[c] == -1){
                scores[0] -= 1;
                scores[1] += 1;
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
                scores[0] -= 1;
                scores[1] += 1;
                matches10[0] += 1;
                matches11[1] += 1;
            }
            else if (part1[c] == -1){
                scores[0] += 1;
                scores[1] -= 1;
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
    string seq1 = "";
    string seq1bis = "";
    string seq2 = "";
    string seq2bis = "";
    for (int c = 0 ; c < part1.size() ; c++){
        if (par2[c] != '?' && part1[c] != 0){
            if (par2[c] == mostFrequent2 || par2[c] == secondFrequent2){
                seq1 += '1'+part1[c];
                seq2 += par2[c];
            }
        }
        seq1bis += '1'+part1[c];
        seq2bis += par2[c];
    }
    auto maxScore = scores[0];
    auto maxScoreIdx = 0; //can be either 0 or 1
    for (auto i = 0 ; i < 2 ; i++){
        //cout << scores[i] << " , ";
        if (scores[i]>=maxScore){
            maxScore = scores[i];
            maxScoreIdx = i;
        }
    }

    res.nmatch = (par1.size()-res.nonComparable+maxScore)/2;
    res.nmismatch = (par1.size()-res.nonComparable-maxScore)/2;
    res.n00 = matches00[maxScoreIdx];
    res.n01 = matches01[maxScoreIdx];
    res.n10 = matches10[maxScoreIdx];
    res.n11 = matches11[maxScoreIdx];
    res.phased = -2*maxScoreIdx + 1;
    //cout << "Computing..." << maxScore << " " << par1.size()-res.nonComparable << " " << res.nmismatch << endl;

    //compute number of expected matches by chance
    int occurences[2] = {maxFrequence2, secondFrequence2};
    // res.easyMatches = maxFrequence*occurences[maxScoreIdx]*(res.nmismatch+res.nmatch)/numberOfBases/numberOfBases 
    //             + (numberOfBases-maxFrequence)*(numberOfBases-occurences[maxScoreIdx])*(res.nmismatch+res.nmatch)/numberOfBases/numberOfBases;
    
    //now compute the chi square
    int n = res.nmatch + res.nmismatch;
    float pmax1 = maxFrequence/numberOfBases;
    float pmax2 = float(occurences[maxScoreIdx]) / (maxFrequence2+secondFrequence2);
    if (pmax1*pmax2*(1-pmax1)*(1-pmax2) == 0){ //if there is only one base in one partition, it can't be compared
        res.chisquare = 0;
        return res;
    }
    //chi square test with 1 degree of freedom
    res.chisquare = pow((matches00[maxScoreIdx]-(1-pmax1)*(1-pmax2)*n),2)/((1-pmax1)*(1-pmax2)*n)
                        + pow((matches01[maxScoreIdx]-(1-pmax1)*pmax2*n),2)/((1-pmax1)*pmax2*n)
                        + pow((matches10[maxScoreIdx]-pmax1*(1-pmax2)*n),2)/(pmax1*(1-pmax2)*n)
                        + pow((matches11[maxScoreIdx]-pmax1*pmax2*n),2)/(pmax1*pmax2*n);
    
    if (float(res.nmismatch)/(res.nmatch+res.nmismatch) <= errorRate*2 ){
        // cout << "To compare : " << endl << seq1 << " ; " << seq1bis << endl << seq2 << " ; "<< seq2bis  << endl;
        // cout << "res : " << res.nmatch << "," << res.nmismatch << "," << res.nonComparable << ", easy matching : " << 
        //            chi_square << endl;

        par1.augmentPartition(newPartitions[maxScoreIdx]);
        // if (par1.number() > 70){
        //     cout << "Intergrating new partition  ";
        //     for (auto i : newPartitions[maxScoreIdx]){cout << i;}
        //     cout << endl;
        //     par1.print();
        // }
        return res;
    }
    else {
        return res;
    }
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

    int occurences[2] = {maxFrequence2, numberOfBases-maxFrequence2};
 
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

    robin_hood::unordered_flat_map <int, int> count; //a map counting how many times a cluster appears

    for (int pos=0 ; pos < listOfFinalPartitions[0].size() ; pos++){
        int id = 0;
        for (auto binary = 0 ; binary < listOfFinalPartitions.size() ; binary++){
            id *= 2;
            auto camp = allPartitions[binary][pos];
            if (camp == 0 || allConfidences[binary][pos] < 0.7){ //0.7 to be pretty confident about the reads we separate we have (more than 70% on all partitions)
                id = -1;
                break;
            }
            else if (camp ==1){
                id += 1;
            }
        }
        if (id>=0){
            frequenceOfPart[id] += 1;
            res[pos] = id;
            count[id] += 1;
        }
    }

    //establish a list of possible clusters
    vector<int> listOfGroups;
    for (auto group : count){
        if (group.second > 1){
            listOfGroups.push_back(group.first);
        }
    }
    //to make sure we don't over-estimate the number of clusters, we'll make the assumption that all haplotypes have haplotype-specific mutation => the number of final cluster cannot be higher that the number of partition
    vector <int> minElements (listOfFinalPartitions.size());
    std::partial_sort_copy(listOfGroups.begin(),  listOfGroups.end(), minElements.begin(), minElements.end());
    int minNumberOfReads = minElements[listOfFinalPartitions.size()-1];

    //now most reads should be assigned to a cluster. Rescue those that have not been assigned or assigned to a singleton cluster

    for (auto read = 0 ; read < res.size() ; read++){

        if (res[read] == -1 || count[res[read]] < max(minNumberOfReads, 2)){ //this means the read needs to be rescued 

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


    cout << "Res : ";
    for (auto i : res){
        cout << i << ",";
    }
    cout << endl;

    return res;

}