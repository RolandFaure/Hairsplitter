
#include <iostream>
#include <string>
#include <chrono>

#include "edlib.h"
#include "input_output.h"
#include "check_overlaps.h"
#include "modify_gfa.h"
#include "clipp.h" //library to build command line interfaces
#include "phase_variants.h"

using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::vector;
using namespace std::chrono;
using namespace clipp;

//../../code/build/OverlapCheck -a alignments.paf -i alignments_on_polished.paf -r assembly_polished.fasta -o alignments_filtered.paf -f nanopore_medium.fq 

void check_dependancies(){

    
    auto perm = system("mkdir tmp 2> trash.txt");
    system("rm trash.txt 2> tmp/trash.txt");

    auto res = system("minimap2 --version > tmp/trash.txt 2> tmp/trash.txt");
    auto miniasm = system("miniasm -V > tmp/trash.txt 2> tmp/trash.txt");
    auto racon = system("racon --version > tmp/trash.txt 2> tmp/trash.txt");
    
    if (res != 0){
        cout << "MISSING DEPENDANCY: minimap2" << endl;
    }
    if (miniasm != 0){
        cout << "MISSING DEPENDANCY: miniasm" << endl;
    }
    if (racon != 0){
        cout << "MISSING DEPENDANCY: racon" << endl;
    }

    if (res != 0 || miniasm != 0 || racon != 0){
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{   
    check_dependancies();

    string fastqfile, allfile, refFile, alnOnRefFile, outputFile, samFile, vcfFile;
    bool polish = false;
    auto cli = (
            required("-f", "--fastq") & opt_value("fastqfile", fastqfile),
            //option("-a", "--all-vs-all").doc("PAF file of the input alignments") & opt_value("all-vs-all file", allfile),
            required("-g", "--gfa") & opt_value("GFA file", refFile),
            required("-i", "--aln-on-ref") & opt_value("aln on ref", alnOnRefFile),
            clipp::option("-p", "--polish").set(polish).doc("Use this option if the assembly is not polished"),
            // option("-s", "--sam") & opt_value("reads aligned on a reference", samFile),
            // option("-v", "--vcf") & opt_value("vcf file", vcfFile),
            required("-o", "--output") & opt_value("output", outputFile)
        );

    if(!parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << make_man_page(cli, argv[0]);
    }
    else {
        std::vector <Read> allreads; 
        robin_hood::unordered_map<std::string, unsigned long int> indices;
        vector<unsigned long int> backbone_reads;
        std::vector <Overlap> allOverlaps;
        vector <Link> allLinks;

        auto t1 = high_resolution_clock::now();
        
        cout << "Parsing reads..." << endl;
        parse_reads(fastqfile, allreads, indices);

        parse_assembly(refFile, allreads, indices, backbone_reads, allLinks);
        parse_PAF(alnOnRefFile, allOverlaps, allreads, indices, backbone_reads, false);

        std::unordered_map<unsigned long int, vector< pair<pair<int,int>, vector<int>> >> partitions;
        cout << "Checking overlaps" << endl;
        std::unordered_map <int, std::pair<int,int>> clusterLimits;
        std::unordered_map <int, vector<pair<int,int>>> readLimits;
        checkOverlaps(allreads, allOverlaps, backbone_reads, partitions, true, clusterLimits, readLimits, polish);
        cout << "Finished checking, now outputting" << endl;

        modify_GFA(refFile, allreads, backbone_reads, allOverlaps, partitions, outputFile, allLinks, clusterLimits, readLimits);
        output_GFA(allreads, backbone_reads, outputFile, allLinks);

        auto t2 = high_resolution_clock::now();

        cout << "Finished in " << duration_cast<seconds>(t2-t1).count() << "s"  << endl;
    }

    // string sequence1 = "ACTGGCTCGTTCGAAAGCTCGT";
    // string sequence2 = "TTACTGGCTCATTCGAAACGCTCGT";
    // string sequence3 = "GCTCGTTGAAAAGCTCGTTGGCT";

    // EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(),
    //                                     edlibNewAlignConfig(42, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    // cout << edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD) << endl;  

    
    return 0;
}
