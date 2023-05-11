#include <algorithm>
#include <iostream>
#include <string>
#include <chrono>
#include <fstream>

#include "edlib.h"
#include "input_output.h"
#include "split_reads.h"
#include "modify_gfa.h"
#include "reassemble_unaligned_reads.h"
#include "clipp.h" //library to build command line interfaces
#include "tools.h"
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
string HAIRSPLITTER;
string WTDBG2;
string SAMTOOLS;
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
    
    if (res != 0){
        cout << "MISSING DEPENDANCY: minimap2" << endl;
    }

    com = " --version > tmp/trash.txt 2> tmp/trash.txt";
    command = WTDBG2 + com;
    auto wtdbg2 = system(command.c_str());
    if (wtdbg2 != 0){
        cout << "MISSING DEPENDANCY: wtdbg2. Proceeding without re-assembling unaligned reads. (was trying command line " << command << ")" << endl;
        WTDBG2 = "no_wtdbg2";
    }

    com = " --version > tmp/trash.txt 2> tmp/trash.txt";
    command = SAMTOOLS + com;
    auto samtools = system(command.c_str());
    if (samtools != 0){
        cout << "MISSING DEPENDANCY: samtools. Proceeding without re-assembling unaligned reads. (was trying command line " << command << ")" << endl;
        WTDBG2 = "no_wtdbg2";
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
            cout << "MISSING DEPENDANCY: could not run graphunzip. Make sure the path to graphunzip is correct and that you have python3, numpy and scipy installed. " 
            << "I was looking for GraphUnzip as: " << GRAPHUNZIP << endl;
        }
    }


    if (res != 0 || racon != 0 || graphunzip != 0 || python3 != 0){
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{   

    string fastqfile, refFile, outputFolder, samFile, vcfFile;
    string alnOnRefFile = "no_file";
    string path_minimap = "minimap2";
    string path_miniasm = "miniasm";
    string path_racon = "racon";
    string path_wtdbg2 = "no_wtdbg2";
    string path_samtools = "samtools";

    //try linking to the local copy of GraphUnzip that comes with HairSplitter
    string path_path = argv[0];
    string path_graphunzip = path_path.substr(0, path_path.size()-18) + "GraphUnzip/graphunzip.py";
    string path_cut_gfa = path_path.substr(0, path_path.size()-18) + "GraphUnzip/cut_gfa.py";

    int num_threads = 1;
    bool polish = false;
    bool dont_simplify = false;
    bool output_read_groups = false;
    bool version = false;
    DEBUG = false;
    bool force = false;
    auto cli = (
            required("-f", "--fastq").doc("Sequencing reads") & value("raw reads", fastqfile),
            required("-i", "--assembly").doc("Original assembly in GFA or FASTA format") & value("assembly", refFile),
            // option("-s", "--sam") & opt_value("reads aligned on a reference", samFile),
            // option("-v", "--vcf") & opt_value("vcf file", vcfFile),
            required("-o", "--output").doc("Output directory") & value("output directory", outputFolder),
            clipp::option("-a", "--aln-on-asm").doc("Reads aligned on assembly (SAM format)") & value("aligned reads", alnOnRefFile),
            clipp::option("-q", "--output-read-groups").set(output_read_groups).doc("Output read groups (txt format)"),
            // clipp::option("-q", "--outputGAF").doc("Output GAF file") & value("output GAF", outputGAF),
            clipp::option("-p", "--polish").set(polish).doc("Use this option if the assembly is not polished"),
            clipp::option("-t", "--threads").doc("Number of threads") & value("threads", num_threads),
            clipp::option("-s", "--dont_simplify").set(dont_simplify).doc("Don't rename the contigs and don't merge them"),
            clipp::option("--path-to-minimap2").doc("Path to the executable minimap2 (if not in PATH)") & value("path to minimap2", path_minimap),
            // clipp::option("--path-to-miniasm").doc("Path to the executable miniasm (if not in PATH)") & value("path to miniasm", path_miniasm),
            clipp::option("--path-to-racon").doc("Path to the executable racon (if not in PATH)") & value("path to racon", path_racon),
            clipp::option("--path-to-wtdbg2").doc("Path to wtdbg2 (if not in PATH)") & value("path to wtdbg2 executable (empty path to disable assembly of unaligned reads)", path_wtdbg2),
            clipp::option("--path-to-samtools").doc("Path to samtools (if not in PATH)") & value("path to samtools", path_samtools),
            clipp::option("-F", "--force").set(force).doc("Force overwrite of output folder if it exists"),
            clipp::option("-v", "--version").set(version),
            clipp::option("-d", "--debug").set(DEBUG)
        );

    if(!parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << make_man_page(cli, argv[0]);
    }
    else {

        if (version){
            cout << "HairSplitter version 1.1.0" << endl;
            exit(EXIT_SUCCESS);
        }

        MINIMAP = path_minimap;
        RACON = path_racon;
        GRAPHUNZIP = path_graphunzip;
        WTDBG2 = path_wtdbg2;
        HAIRSPLITTER = path_path.substr(0, path_path.size()-18);
        SAMTOOLS = path_samtools;

        //strip the last / if it exists in the output folder
        if (outputFolder[outputFolder.size()-1] == '/'){
            outputFolder = outputFolder.substr(0, outputFolder.size()-1);
        }

        string outputGAF = outputFolder+"/tmp/outout.gaf";
        string outputFile = outputFolder+"/hairsplitter_assembly.gfa";

        string mkdir_out = "mkdir "+ outputFolder + " && mkdir "+ outputFolder + "/tmp/";
        auto mkdir = system(mkdir_out.c_str());
        if (mkdir != 0 && !force){
            cout << "Could not run command line " << mkdir_out << endl;
            cout << "ERROR: could not create output folder \"" << outputFolder << "\" make sure it does not exist already, or use option -F." << endl;
            exit(EXIT_FAILURE);
        }

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
            //convert the fasta file to gfa
            cout << "Converting the fasta file to gfa" << endl;

            string newRefFile = outputFolder+"/tmp/assembly.gfa";
            convert_FASTA_to_GFA(refFile, newRefFile);
            refFile = newRefFile;
            format = "fasta";
        }
        else{
            cout << "ERROR: Unrecognized extension for input assembly. Authorized file extensions are .gfa (for GFA) or .fa/.fasta (for fasta)" << endl;
            exit(EXIT_FAILURE);
        }

        cout << "\n\t******************\n\t*                *\n\t*  Hairsplitter  *\n\t*    Welcome!    *\n\t*                *\n\t******************\n\n";
        cout << "-- Please note that details on what Hairsplitter does will be jotted down in file "+outputFolder+"/hairsplitter_summary.txt --\n";

        //generate the paf file if not already generated
        if (alnOnRefFile == "no_file"){

            cout <<  "\n===== STAGE 1: Aligning reads on the reference\n\n";

            alnOnRefFile = outputFolder+"/tmp/reads_aligned_on_assembly.sam";

            //cut the gfa in small contigs to speed up the computation
            string command_cut = "python " + path_cut_gfa + " -a " + refFile + " -l 100000 -o "+outputFolder+"/tmp/assembly_cut.gfa > "+outputFolder+"/tmp/logcut.txt";
            auto res_cut = system(command_cut.c_str());
            if (res_cut != 0){
                cout << "ERROR while running " << command_cut << endl;
                cout << "ERROR: Hairsplitter needs python3 to run. Make sure to have python3 installed. If this is not \
                    the problem, you may circumvent the error by using option -a." << endl;
                exit(EXIT_FAILURE);
            }

            refFile = outputFolder+"/tmp/assembly_cut.gfa";
            string fastaFile = outputFolder+"/tmp/assembly.fa";

            convert_GFA_to_FASTA(refFile, fastaFile);

            string command = MINIMAP + " " + fastaFile + " " + fastqfile + " -ax map-ont --secondary=no -t "+ std::to_string(num_threads) +" > " + alnOnRefFile + " 2> "+outputFolder+"/tmp/logminimap.txt";
            cout << " - Running minimap with command line:\n     " << command << "\n   The log of minimap2 can be found at "+outputFolder+"/tmp/logminimap.txt" << endl;
            auto res_minimap = system(command.c_str());
            if (res_minimap != 0){
                cout << "ERROR while running " << command << endl;
                cout << "ERROR: minimap2 could not run properly, check "+outputFolder+"/tmp/logminimap.txt" << endl;
                exit(EXIT_FAILURE);
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
        parse_assembly(refFile, allreads, indices, backbone_reads, allLinks);

        cout << " - Loading alignments of the reads on the contigs from " << alnOnRefFile << "\n";
        if (alnOnRefFile.substr(alnOnRefFile.size()-4,4) == ".paf"){
            cout << "ERROR: please provide a .sam file as input for the alignments of the reads on the contigs." << endl;
            exit(EXIT_FAILURE);
            // parse_PAF(alnOnRefFile, allOverlaps, allreads, indices, backbone_reads, false);
        }
        else if (alnOnRefFile.substr(alnOnRefFile.size()-4,4) == ".sam"){
            parse_SAM(alnOnRefFile, allOverlaps, allreads, indices);
        }
        else{
            cout << "ERROR: the file containing the alignments on the assembly should be either .paf or .sam" << endl;
            exit(EXIT_FAILURE);
        }

        if (WTDBG2 != "no_wtdbg2"){
            cout << " - Re-assembling unaligned reads (to skip this step use option --path-to-wtdbg2 0)\n";
            reassemble_unaligned_reads(allreads, allOverlaps, fastqfile, backbone_reads, outputFolder, num_threads, indices, allLinks);
        }

        cout << "\n===== STAGE 3: Checking every contig and separating reads when necessary\n\n";
        cout << " For each contig I am going to:\n  - Align all reads precisely on the contig\n  - See at what positions there seem to be many reads disagreeing\n";
        cout << "  - See if the reads disagreeing seem always to be the same ones\n  - If they are always the same ones, they probably come from another haplotype, so separate!\n\n";
        cout << " *To see in more details what contigs have been separated, check out the hairsplitter_summary.txt in the output folder*\n";

        std::unordered_map<unsigned long int, vector< pair<pair<int,int>, pair<vector<int>, std::unordered_map<int, string>>  > >> partitions;
        std::unordered_map <int, vector<pair<int,int>>> readLimits;
        string tmpFolder = outputFolder+"/tmp/";
        float errorRate;
        split_contigs(fastqfile, allreads, allOverlaps, backbone_reads, partitions, true, readLimits, polish, num_threads, tmpFolder, errorRate);

        //output GAF, the path of all reads on the new contigs
        // cout << "Outputting GAF" << endl;
        output_GAF(allreads, backbone_reads, allLinks, allOverlaps, partitions, outputGAF);
        
        cout << "\n===== STAGE 4: Creating and polishing all the new contigs\n\n This can take time, as we need to polish every new contig using Racon\n";

        string readGroupsFile = tmpFolder + "read_groups.txt";
        cout << " - Outputting how reads are partitionned into groups in file " << readGroupsFile << "\n";

        output_readGroups(readGroupsFile, allreads, backbone_reads, partitions, allOverlaps);
        
        
        modify_GFA(fastqfile, allreads, backbone_reads, allOverlaps, partitions, allLinks, readLimits, num_threads, tmpFolder, errorRate);
        string zipped_GFA = outputFolder+"/tmp/zipped_gfa.gfa";
        output_GFA(allreads, backbone_reads, zipped_GFA, allLinks);

        //now "unzip" the assembly, improving the contiguity where it can be improved
        cout << "\n===== STAGE 5: Linking all the new contigs that have been produced (maybe bridging repeated regions)\n\n";
        string simply = "";
        if (dont_simplify){
            simply = " --dont_merge -r";
        }
        string com = " unzip -l " + outputGAF + " -g " + zipped_GFA + simply + " -o " + outputFile + " 2>"+outputFolder+"/tmp/logGraphUnzip.txt >"+outputFolder+"/tmp/trash.txt";
        string command = GRAPHUNZIP + com;
        cout << " - Running GraphUnzip with command line:\n     " << command << "\n   The output of GraphUnzip is dumped on "+outputFolder+"/tmp/logGraphUnzip.txt\n";
        int resultGU = system(command.c_str());
        if (resultGU != 0){
            cout << "ERROR: GraphUnzip failed. Please check the output of GraphUnzip in "+outputFolder+"/tmp/logGraphUnzip.txt" << endl;
            exit(EXIT_FAILURE);
        }

        cout << "\n *To see in more details what supercontigs were created with GraphUnzip, check the hairsplitter_summary.txt*\n";
        string output = "output.txt";
        std::ofstream o(output, std::ios_base::app);//appending to the file
        o << "\n\n *****Linking the created contigs***** \n\nLeft, the name of the produced supercontig. Right, the list of new contigs with a suffix -0, -1...indicating the copy of the contig, linked with _ \n\n";
        o.close();
        command = "cat output.txt supercontigs.txt > output2.txt 2> "+outputFolder+"/tmp/trash.txt";
        system(command.c_str());
        command =  "mv output2.txt "+outputFolder+"/hairsplitter_summary.txt && rm supercontigs.txt output.txt 2> "+outputFolder+"/tmp/trash.txt";
        system(command.c_str());
        
        //convert the output to fasta
        string fasta_name = outputFile.substr(0, outputFile.size()-4) + ".fasta";
        convert_GFA_to_FASTA(outputFile, fasta_name);
        
        auto t2 = high_resolution_clock::now();

        cout << "\n\nFinished in " << duration_cast<seconds>(t2-t1).count() << "s. If you experienced any trouble running Hairsplitter, contact us via an issue on github, github.com/RolandFaure/hairsplitter :-)\n\n"  << endl;
    }
    
    return 0;
}
