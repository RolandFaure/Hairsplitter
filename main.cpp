
#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

#include "edlib.h"
#include "input_output.h"
#include "check_overlaps.h"
#include "modify_gfa.h"
#include "clipp.h" //library to build command line interfaces
//#include "phase_variants.h"

using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::vector;
using namespace std::chrono;
using namespace clipp;

string MINIMAP;
string RACON;
string GRAPHUNZIP;
bool DEBUG;

//../../code/build/OverlapCheck -a alignments.paf -i alignments_on_polished.paf -r assembly_polished.fasta -o alignments_filtered.paf -f nanopore_medium.fq 

void check_dependancies(){

    
    auto perm = system("mkdir tmp 2> trash.txt");
    system("rm trash.txt 2> tmp/trash.txt");

    string com = " --version > tmp/trash.txt 2> tmp/trash.txt";
    auto command = MINIMAP + com;
    auto res = system(command.c_str());

    // com = " -V > tmp/trash.txt 2> tmp/trash.txt";
    // command = MINIASM + com;
    // auto miniasm = system(command.c_str());

    command = "python3 --version > tmp/trash.txt 2> tmp/trash.txt";
    auto python3 = system(command.c_str());

    com = " --version > tmp/trash.txt 2> tmp/trash.txt";
    command = RACON + com;
    auto racon = system(command.c_str());

    com = " -h > tmp/trash.txt 2> tmp/trash.txt";
    command = GRAPHUNZIP + com;
    auto graphunzip = system(command.c_str());

    command = "awk --help > tmp/trash.txt 2> tmp/trash.txt";
    auto awk = system(command.c_str());
    
    if (res != 0){
        cout << "MISSING DEPENDANCY: minimap2" << endl;
    }
    // if (miniasm != 0){
    //     cout << "MISSING DEPENDANCY: miniasm" << endl;
    // }
    if (racon != 0){
        cout << "MISSING DEPENDANCY: racon" << endl;
    }

    if (python3 != 0){
        cout << "MISSING DEPENDANCY: python3" << endl;
    }
    else if (graphunzip != 0){
        command = "GraphUnzip -h > tmp/trash.txt 2> tmp/trash.txt";
        auto graphunzip2 = system(command.c_str());
        if (graphunzip2 == 0){
            graphunzip = 0;
            GRAPHUNZIP = "GraphUnzip";
        }
        else{
            cout << "MISSING DEPENDANCY: could not run graphunzip. Make sure the path to graphunzip is correct and that you have python3, numpy and scipy installed" << endl;
        }
    }

    if (awk != 0){
        cout << "WARNING: awk could not be found on this computer. Make sure to align the reads on the assembly by yourself and then run Hairsplitter with option -a" << endl;
    }


    if (res != 0 || racon != 0 || graphunzip != 0 || python3 != 0){
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{   
    string fastqfile, refFile, outputFile, samFile, vcfFile;
    string alnOnRefFile = "no_file";
    string outputGAF = "tmp/outout.gaf";
    string path_minimap = "minimap2";
    string path_miniasm = "miniasm";
    string path_racon = "racon";

    //try linking to the local copy of GraphUnzip that comes with HairSplitter
    string path_graphunzip_path = argv[0];
    string path_graphunzip = path_graphunzip_path.substr(0, path_graphunzip_path.size()-12) + "GraphUnzip/graphunzip.py";

    int num_threads = 1;
    bool polish = false;
    bool dont_simplify = false;
    DEBUG = false;
    auto cli = (
            required("-f", "--fastq").doc("Sequencing reads") & value("raw reads", fastqfile),
            required("-i", "--assembly").doc("Original assembly in GFA or FASTA format") & value("assembly", refFile),
            // option("-s", "--sam") & opt_value("reads aligned on a reference", samFile),
            // option("-v", "--vcf") & opt_value("vcf file", vcfFile),
            required("-o", "--outputGFA").doc("Output assembly file, same format as input") & value("output assembly", outputFile),
            clipp::option("-a", "--aln-on-asm").doc("Reads aligned on assembly (PAF format)") & value("aligned reads", alnOnRefFile),
            // clipp::option("-q", "--outputGAF").doc("Output GAF file") & value("output GAF", outputGAF),
            clipp::option("-p", "--polish").set(polish).doc("Use this option if the assembly is not polished"),
            clipp::option("-t", "--threads").doc("Number of threads") & value("threads", num_threads),
            clipp::option("-s", "--dont_simplify").set(dont_simplify).doc("Don't rename the contigs and don't merge them"),
            clipp::option("--path-to-minimap2").doc("Path to the executable minimap2 (if not in PATH)") & value("path to minimap2", path_minimap),
            // clipp::option("--path-to-miniasm").doc("Path to the executable miniasm (if not in PATH)") & value("path to miniasm", path_miniasm),
            clipp::option("--path-to-racon").doc("Path to the executable racon (if not in PATH)") & value("path to racon", path_racon),
            clipp::option("--path-to-graphunzip").doc("Path to graphunzip.py (if not in PATH)") & value("path to graphunzip", path_graphunzip),
            clipp::option("-d", "--debug").set(DEBUG)
        );

    if(!parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << make_man_page(cli, argv[0]);
    }
    else {

        MINIMAP = path_minimap;
        RACON = path_racon;
        GRAPHUNZIP = path_graphunzip;

        system("mkdir tmp/ 2> trash.txt");
        check_dependancies();

        std::vector <Read> allreads; 
        robin_hood::unordered_map<std::string, unsigned long int> indices;
        vector<unsigned long int> backbone_reads;
        std::vector <Overlap> allOverlaps;
        vector <Link> allLinks;

        auto t1 = high_resolution_clock::now();

        string format = "";
        if (refFile.substr(refFile.size()-4, 4) == ".gfa"){
            format = "gfa";
        }
        else if (refFile.substr(refFile.size()-3, 3) == ".fa" || refFile.substr(refFile.size()-6, 6) == ".fasta"){
            format = "fasta";
        }
        else{
            cout << "ERROR: Unrecognized extension for input assembly. Authorized file extensions are .gfa (for GFA) or .fa/.fasta (for fasta)" << endl;
            exit(EXIT_FAILURE);
        }

        cout << "\n\t******************\n\t*                *\n\t*  Hairsplitter  *\n\t*    Welcome!    *\n\t*                *\n\t******************\n\n";
        cout << "-- Please note that details on what Hairsplitter does will be jotted down in file output.txt --\n";

        //generate the paf file if not already generated
        if (alnOnRefFile == "no_file"){

            cout <<  "\n===== STAGE 1: Aligning reads on the reference\n\n";

            alnOnRefFile = "reads_aligned_on_assembly.paf";

            if (format == "gfa"){
                string fastaFile = "tmp/"+refFile.substr(0, refFile.size()-4)  +".fa";
                string command = "awk '/^S/{print \">\"$2\"\\n\"$3}' " + refFile + " > " + fastaFile;
                auto awk = system(command.c_str());
                if (awk != 0){
                    cout << "DEPENDANCY ERROR: Hairsplitter needs awk to run without using option -a. Please install awk or use option -a." << endl;
                    exit(EXIT_FAILURE);
                }
                command = MINIMAP + " " + fastaFile + " " + fastqfile + " -x map-ont --secondary=no > " + alnOnRefFile + " 2> tmp/logminimap.txt";
                cout << " - Running minimap with command line:\n     " << command << "\n   The output of minimap2 is dumped on tmp/logminimap.txt\n";
                system(command.c_str());
            }
            else{
                string command = MINIMAP + " " + refFile + " " + fastqfile + " -t " + num_threads + " -x map-ont --secondary=no > " + alnOnRefFile + " 2> tmp/logminimap.txt"; 
                cout << " - Running minimap with command line:\n     " << command << "\n   The log of minimap2 is being dumped on tmp/logminimap.txt\n";
                system(command.c_str());
            }
        }
        else{
            cout <<  "\n===== STAGE 1: Aligning reads on the reference\n\n";
            cout << "  Skipped because alignments already inputted with option -a.\n";
        }

        cout << "\n===== STAGE 2: Loading data\n\n";

        cout << " - Loading all reads from " << fastqfile << " in memory\n";
        parse_reads(fastqfile, allreads, indices);

        cout << " - Loading all contigs from " << refFile << " in memory\n";
        parse_assembly(refFile, allreads, indices, backbone_reads, allLinks, format);

        cout << " - Loading alignments of the reads on the contigs from " << alnOnRefFile << "\n";
        parse_PAF(alnOnRefFile, allOverlaps, allreads, indices, backbone_reads, false);

        cout << "\n===== STAGE 3: Checking every contig and separating reads when necessary\n\n";
        cout << " For each contig I am going to:\n  - Align all reads precisely on the contig\n  - See at what positions there seem to be many reads disagreeing\n";
        cout << "  - See if the reads disagreeing seem always to be the same ones\n  - If they are always the same ones, they probably come from another haplotype, so separate!\n\n";
        cout << " *To see in more details what contigs have been separated, check the output.txt*\n";

        std::unordered_map<unsigned long int, vector< pair<pair<int,int>, vector<int>> >> partitions;
        std::unordered_map <int, vector<pair<int,int>>> readLimits;
        checkOverlaps(allreads, allOverlaps, backbone_reads, partitions, true, readLimits, polish, num_threads);

        //output GAF, the path of all reads on the new contigs
        if (format == "gfa"){
            // cout << "Outputting GAF" << endl;
            output_GAF(allreads, backbone_reads, allLinks, allOverlaps, partitions, outputGAF);
        }
        
        cout << "\n===== STAGE 4: Creating and polishing all the new contigs\n\n This can take time, as we need to polish every new contig using Racon\n";
        if (format == "gfa"){
            modify_GFA(refFile, allreads, backbone_reads, allOverlaps, partitions, allLinks, readLimits, num_threads);
            string zipped_GFA = "zipped_gfa.gfa";
            output_GFA(allreads, backbone_reads, zipped_GFA, allLinks);

            //now "unzip" the assembly, improving the contiguity where it can be improved
            cout << "\n===== STAGE 5: Linking all the new contigs that have been produced (maybe bridging repeated regions)\n\n";
            string simply = "";
            if (dont_simplify){
                simply = " --dont_merge -r";
            }
            string com = " unzip -l " + outputGAF + " -g " + zipped_GFA + simply + " -o " + outputFile + " >tmp/logGraphUnzip.txt 2>tmp/trash.txt";
            string command = GRAPHUNZIP + com;
            cout << " - Running GraphUnzip with command line:\n     " << command << "\n   The output of GraphUnzip is dumped on tmp/logGraphUnzip.txt\n";
            system(command.c_str());

            cout << "\n *To see in more details what supercontigs were created with GraphUnziop, check the output.txt*\n";
            string output = "output.txt";
            std::ofstream o(output, std::ios_base::app);//appending to the file
            o << "\n\n *****Linking the created contigs***** \n\nLeft, the name of the produced supercontig. Right, the list of new contigs with a suffix -0, -1...indicating the copy of the contig, linked with _ \n\n";
            o.close();
            command = "cat output.txt supercontigs.txt > output2.txt 2> tmp/trash.txt";
            system(command.c_str());
            command =  "mv output2.txt output.txt & rm supercontigs.txt 2> tmp/trash.txt";
            system(command.c_str());
        }
        else if (format == "fasta"){
            modify_FASTA(refFile, allreads, backbone_reads, allOverlaps, partitions, readLimits, num_threads);
            output_FASTA(allreads, backbone_reads, outputFile);
        }
        
        auto t2 = high_resolution_clock::now();

        cout << "\n\nFinished in " << duration_cast<seconds>(t2-t1).count() << "s. If you experienced any trouble running Hairsplitter, contact us via an issue on github, github.com/RolandFaure/hairsplitter :-)\n\n"  << endl;
    }
    
    return 0;
}
