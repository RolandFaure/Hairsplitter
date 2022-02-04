
#include <iostream>
#include <string>
#include <chrono>

#include "edlib.h"
#include "input_output.h"
#include "check_overlaps.h"

using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::vector;
using namespace std::chrono;

int main(int argc, char *argv[])
{

    string sequence1 = "ACTGGCTCGTTCGAAAGCTCGT";
    string sequence2 = "TTACTGGCTCATTCGAAACGCTCGT";
    string sequence3 = "GCTCGTTGAAAAGCTCGTTGGCT";

    // EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(),
    //                                     edlibNewAlignConfig(42, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    // cout << edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD) << endl;  

    std::vector <Read> allreads; 
    robin_hood::unordered_map<std::string, unsigned long int> indices;
    vector<unsigned long int> backbone_reads;

    auto t1 = high_resolution_clock::now();
    //parse_reads("/home/rfaure/Documents/these/overlap_filtering/mock2.fasta", allreads, indices);
    parse_reads("/home/rfaure/Documents/these/overlap_filtering/reads.fq", allreads, indices);

    std::vector <Overlap> allOverlaps;
    //parse_PAF("/home/rfaure/Documents/these/overlap_filtering/mock2_alignments.paf", allOverlaps, allreads, indices, backbone_reads);
    parse_PAF("/home/rfaure/Documents/these/overlap_filtering/alignments.paf", allOverlaps, allreads, indices, backbone_reads);

    checkOverlaps(allreads, allOverlaps, backbone_reads);

    auto t2 = high_resolution_clock::now();

    cout << "Finished in " << duration_cast<milliseconds>(t2-t1).count() << "ms"  << endl;


    return 0;
}
