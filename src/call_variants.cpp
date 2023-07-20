#include "call_variants.h"

#include <iostream>
#include <fstream> //for reading files
#include <sstream> //for parsing strings
#include <string>
#include <omp.h> //for efficient parallelization

#include "input_output.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ifstream;
using robin_hood::unordered_map;
using std::pair;
using std::make_pair;


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
    double totalLengthOfAlignment = 0;

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
        if (i != previous_char){
            previous_previous_previous_char = previous_previous_char;
            previous_previous_char = previous_char;
            previous_char = i;
        }
        newRef += (unsigned char) ('!' + 5*ACGT.find(previous_previous_previous_char) + ACGT.find(previous_previous_char) + 25*ACGT.find(previous_char));
    }
    newref = newRef;

    /*
    //print snps (just for debugging)
    cout << "Printing SNPs, in split_read: cicizzx" << endl;
    int step = 1; //porportions of reads
    int prop = 1; //proportion of positions
    int firstRead = 0;
    int lastRead = polishingReads.size();
    int numberOfReads = lastRead-firstRead;
    int start = 0;
    int end = 500;
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
            cout << neighbor << " " << index << " " << allreads[allOverlaps[allreads[bbcontig].neighbors_[index]].sequence1].name << endl;
        }
        index+= step;
    }
    int n =start;
    for(unsigned char i : consensus.substr(start, end-start)){
        cout << newRef[n];
        n+=prop;
    } cout << endl;
    // cout << "meanDistance : " << totalDistance/totalLengthOfAlignment << endl;
    // exit(1);
    */
    
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
void call_variants(
    std::vector<Column> &snps, 
    std::vector <Read> &allreads,
    std::vector <Overlap> &allOverlaps, 
    long int contig,
    std::string &ref,
    std::vector<size_t> &suspectPostitions, 
    float &meanError,
    std::string &tmpFolder,
    std::string &outputFile,
    bool DEBUG){

    std::ofstream out(outputFile, std::ios_base::app);

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
        unordered_map<char, int> content;
        // bool here496 = false;
        for (short n = 0 ; n < snps[position].content.size() ; n++){
            char base = snps[position].content[n];
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
        vector<pair<char, int>> content_sorted;
        for (auto it = content.begin() ; it != content.end() ; it++){
            content_sorted.push_back(make_pair(it->first, it->second));
        }
        std::sort(content_sorted.begin(), content_sorted.end(), [](const pair<char, int>& a, const pair<char, int>& b) {return a.second > b.second;});

        if (content_sorted[1].second > minimumNumberOfReadsToBeConsideredSuspect //there is a "frequent" base other than the ref base
            && (content_sorted[1].second > content_sorted[2].second * 5 || minimumNumberOfReadsToBeConsideredSuspect == 2) //this other base occurs much more often than the third most frequent base (could it really be chance ?)
            && content_sorted[0].first%5 != content_sorted[1].first%5 //the central base differs
            && position - posoflastsnp > 5){ //the snp is not too close to the previous one
            
            posoflastsnp = position;
            suspectPositions.push_back(position);

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

    #pragma omp critical
    {
        auto coverage = depthOfCoverage / ref.size();
        //output the suspect positions
        out << endl << "CONTIG\t" << allreads[contig].name << "\t" << allreads[contig].sequence_.size() << "\t" << coverage << "\n";
        for (long int n : allreads[contig].neighbors_){
            out << "READ\t" << allreads[ allOverlaps[n].sequence1 ].name 
                << "\t" << allOverlaps[n].position_1_1 << "\t" << allOverlaps[n].position_1_2 
                << "\t" << allOverlaps[n].position_2_1 << "\t" << allOverlaps[n].position_2_2
                << "\t" << allOverlaps[n].strand << "\n";
        }
        //output the list of reads
        for (auto c : suspiciousColumns){
            out << "SNPS\t" << c.pos << "\t" << c.ref_base << "\t" << c.second_base << "\t:";
            int idx = 0;
            for (auto r = 0 ; r < c.readIdxs.size() ; r++){
                while (idx < c.readIdxs[r]){
                    out << " ";
                    idx+= 1;
                }
                out << c.content[r];
                idx+=1;
            }
            out << "\n";
        }
    }
}

int main(int argc, char *argv[])
{
    std::vector <Read> allreads; 
    robin_hood::unordered_map<std::string, unsigned long int> indices;
    vector<unsigned long int> backbone_reads;
    std::vector <Overlap> allOverlaps;
    vector <Link> allLinks;

    if (argc < 9){
        std::cout << "Usage: ./call_variants <gfa_file> <reads_file> <sam_file> <num_threads> <tmpDir> <error_rate_out> <DEBUG> <file_out>\n";
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
    //erase the output file
    std::ofstream out(file_out);
    out.close();

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

    cout << " - Calling variants on each contig using basic pileup\n";

    int index = 0;
    float totalErrorRate = 0;

    omp_set_num_threads(num_threads);
    #pragma omp parallel
    {
    #pragma omp for
        for (int contigIndex = 0 ; contigIndex < backbone_reads.size() ; contigIndex++){

            long int contig = backbone_reads[contigIndex];
            if (allreads[contig].name != "edge_8@00"){ //for debugging purposes

                // if (DEBUG){
                //     #pragma omp critical
                //     {
                //         cout << "Looking at contig number " << index << " out of " << backbone_reads.size() << " (" << allreads[contig].name << ")" << ". By thread " 
                //             << omp_get_thread_num() << ", " << allreads[contig].neighbors_.size() << " reads align here.\n";
                //     }
                // }
                
                parse_reads_on_contig(readsFile, contig, allOverlaps, allreads);
    
                vector<Column> snps;  //vector containing list of position, with SNPs at each position
                //build a pileup of the reads on the contig

                robin_hood::unordered_map<int, int> insertionPositions;
                string ref3mers;
                std::unordered_map <int, std::vector<std::pair<int,int>>> readLimits;
                
                float meanDistance = generate_msa(contig, allOverlaps, allreads, snps, insertionPositions, 
                    contigIndex, readLimits, ref3mers, tmpFolder, DEBUG);
                
                #pragma omp critical
                {
                    totalErrorRate += meanDistance;
                }
                

                //call variants
                vector<size_t> suspectPostitions;
                call_variants(snps, allreads, allOverlaps, contig, ref3mers, suspectPostitions, meanDistance, tmpFolder, file_out, DEBUG);

                //output the columns with variants in a file

                //free memory by deleting snps
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
            index++;
        }
    }

    //output the errorRate
    std::ofstream errorRateFile;
    errorRateFile.open(error_rate_out);
    errorRateFile << totalErrorRate/backbone_reads.size() << endl;
    errorRateFile.close();

    return 0;

}

