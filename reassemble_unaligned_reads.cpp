#include "reassemble_unaligned_reads.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

extern string MINIMAP; //path to the minimap executable
extern string WTDBG2;
extern string SAMTOOLS;

/**
 * @brief Reassemble unaligned reads using wtdbg2
 * 
 * @param allreads 
 * @param allOverlaps 
 * @param fileReads 
 * @param backbones_reads 
 * @param outputFolder 
 * @param num_threads 
 */
void reassemble_unaligned_reads(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::string &fileReads,
                std::vector<unsigned long int> &backbones_reads, string outputFolder, int num_threads, 
                robin_hood::unordered_map<std::string, unsigned long int> &indices, std::vector<Link> &allLinks){

    //output all unaligned reads to a file
    ofstream unalignedReadsFile(outputFolder + "/tmp/unaligned_reads.fasta");
    if (!unalignedReadsFile){
        cout << "problem opening " << outputFolder + "/tmp/unaligned_reads.fasta" << ". Do you have the right permissions ?" << endl;
        throw std::invalid_argument( "File could not be written" );
    }

    ifstream fastqfile(fileReads);
     
    for (auto read : allreads){
        if (read.neighbors_.size()==0){

            //load the sequence
            fastqfile.seekg(read.get_position_in_file());
            string line;
            getline(fastqfile, line);
            unalignedReadsFile << read.name + "\n" + line + "\n";
        }
    }
    unalignedReadsFile.close();

    string reads_file = outputFolder + "/tmp/unaligned_reads.fasta";
    assemble_with_wtdbg2(reads_file, outputFolder, num_threads);

    //load the contigs
    string contigs_file = outputFolder + "/tmp/wtdbg2.fa";
    string new_contigs_file = outputFolder + "/tmp/wtdbg2.gfa";
    convert_FASTA_to_GFA(contigs_file, new_contigs_file);
    string sam_file = outputFolder + "/tmp/unaligned_reads_aligned.sam";
    parse_assembly(new_contigs_file, allreads, indices, backbones_reads, allLinks);
    parse_SAM(sam_file, allOverlaps, allreads, indices);

}

/**
 * @brief Assemble and a polish a file of reads using wtdbg2
 * 
 * @param fileReads 
 * @param outputFolder 
 * @param num_threads 
 */
void assemble_with_wtdbg2(std::string &fileReads, std::string outputFolder, int num_threads){

    string wtdbg2_folder = WTDBG2.substr(0, WTDBG2.find_last_of("/"));

    string comAsm = wtdbg2_folder+ "/wtdbg2 -A -e 5 -l 1000 -L 5000 -S 1 -R -o " + outputFolder + "/tmp/wtdbg2 -i " + fileReads + " 2>tmp/trash.txt";
    auto res = system(comAsm.c_str());
    if (res != 0){
        cout << "ERROR wtdbg2 failed, while running " << comAsm << endl;
        exit(1);
    }

    string cons_wtdbg2 = wtdbg2_folder+"/wtpoa-cns -i " + outputFolder + "/tmp/wtdbg2.ctg.lay.gz -fo " + outputFolder + "/tmp/dbg.raw.fa 2>tmp/trash.txt";
    int res_wtdbg2 = system(cons_wtdbg2.c_str());

    // polish consensus, not necessary if you want to polish the assemblies using other tools
    string comMap = MINIMAP+" -ax map-pb -r2k " + outputFolder + "/tmp/dbg.raw.fa " + outputFolder + "/tmp/reads.fasta 2>tmp/trash.txt | "+SAMTOOLS+" sort >tmp/dbg.bam 2>tmp/trash.txt";
    res = system(comMap.c_str());
    if (res != 0){
        cout << "ERROR minimap2 failed, while running " << comMap << endl;
        exit(1);
    }

    string comSamtools = SAMTOOLS+" view -F0x900 " + outputFolder + "/tmp/dbg.bam 2>tmp/trash.txt | "+wtdbg2_folder+"/wtpoa-cns -d " + outputFolder + "/tmp/dbg.raw.fa -i - -fo " + outputFolder + "/tmp/dbg.cns.fa 2>tmp/trash.txt";
    res = system(comSamtools.c_str());
    if (res != 0){
        cout << "ERROR samtools failed, while running " << comSamtools << endl;
        exit(1);
    }

    string new_contigs_file = outputFolder + "/tmp/wtdbg2.fa";
    string comUnfold = "awk '{if(\">\" == substr($1,1,1)){ printf \"\\n\"; print;} else printf $1;}' " + outputFolder + "/tmp/dbg.cns.fa >" + new_contigs_file + " 2>" + outputFolder + "/tmp/trash.txt";
    res = system(comUnfold.c_str());
    if (res != 0){
        cout << "ERROR awk failed, while running " << comUnfold << endl;
        exit(1);
    }

    //rename the reads of the new assembly to be sure not to conflict with the original reads
    string prefix = "reassembled_by_HairSplitter_";
    rename_reads(new_contigs_file, prefix);

    //align the unaligned reads on the new assembly
    string comMap2 = MINIMAP+" -ax map-pb -r2k " + outputFolder + "/tmp/wtdbg2.fa " + fileReads +" > "+ outputFolder +"/tmp/unaligned_reads_aligned.sam 2>" + outputFolder + "/tmp/trash.txt ";
    res = system(comMap2.c_str());
    if (res != 0){
        cout << "ERROR minimap2 failed, while running " << comMap2 << endl;
        exit(1);
    }
}







