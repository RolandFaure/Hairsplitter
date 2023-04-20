#include "reassemble_unaligned_reads.h"
#include <iostream>
#include <fstream>
#include <chrono>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;
using namespace std::chrono;

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
    string outTmp = outputFolder + "/tmp/";
    string id = "0";
    string ref = "";
    assemble_with_wtdbg2(reads_file, outTmp, ref, id);

    //align the unaligned reads on the new assembly
    string comMap2 = MINIMAP+" -ax map-pb -r2k " + outTmp + "wtdbg2_"+id+".fa " + fileReads +" > "+ outTmp +"unaligned_reads_aligned.sam 2>" + outTmp + "trash.txt ";
    auto res = system(comMap2.c_str());
    if (res != 0){
        cout << "ERROR minimap2 failed vt, while running " << comMap2 << endl;
        exit(1);
    }

    //load the contigs
    string contigs_file = outputFolder + "/tmp/wtdbg2_"+id+".fa";
    string new_contigs_file = outputFolder + "/tmp/wtdbg2_"+id+".gfa";
    convert_FASTA_to_GFA(contigs_file, new_contigs_file);
    string sam_file = outputFolder + "/tmp/unaligned_reads_aligned.sam";
    parse_assembly(new_contigs_file, allreads, indices, backbones_reads, allLinks);
    parse_SAM(sam_file, allOverlaps, allreads, indices);

}

/**
 * @brief Assemble and a polish a file of reads using wtdbg2
 * 
 * @param fileReads All the reads to assemble
 * @param outputFolder 
 * @param ref File to the reference
 * @param num_threads 
 */
void assemble_with_wtdbg2(std::string &fileReads, std::string outputFolder, std::string &ref, std::string id){

    auto lastslash = WTDBG2.find_last_of("/");
    if (lastslash == string::npos){
        lastslash = 0;
    }
    string wtdbg2_folder = WTDBG2.substr(0, lastslash);

    auto t0 = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();

    if (ref == ""){
        string comAsm = wtdbg2_folder+ "/wtdbg2 -A -e 1 -l 200 -L 0 -S 100 --no-read-clip --no-chainning-clip --ctg-min-length 200 --ctg-min-nodes 0 -R -o "
                                + outputFolder + "wtdbg2_"+id+" -i " + fileReads + " 2>"+outputFolder+"trash.txt";
        auto res = system(comAsm.c_str());
        if (res != 0){
            cout << "ERROR wtdbg2 failed, while running " << comAsm << endl;
            exit(1);
        }


        t1 = high_resolution_clock::now();

        string cons_wtdbg2 = wtdbg2_folder+"/wtpoa-cns -t 1 -i " + outputFolder + "wtdbg2_"+id+".ctg.lay.gz -fo " + outputFolder + "dbg_"+id+".raw.fa 2>tmp/trash.txt";
        int res_wtdbg2 = system(cons_wtdbg2.c_str());
        t2 = high_resolution_clock::now();
        ref = outputFolder + "dbg_"+id+".raw.fa";
    }

    // polish consensus, not necessary if you want to polish the assemblies using other tools
    // string comMap = MINIMAP+" -ax map-pb -r2k " + outputFolder + "dbg_"+id+".raw.fa " + fileReads + " 2>" + outputFolder + "trash.txt | "
    //                     +SAMTOOLS+" sort >" + outputFolder + "dbg_"+id+".bam 2>" + outputFolder + "trash.txt";
    string comMap = MINIMAP+" -ax map-pb " + ref + " " + fileReads + " 2>" + outputFolder + "trash.txt | "
                        + SAMTOOLS+" sort >" + outputFolder + "dbg_"+id+".bam 2>" + outputFolder + "trash.txt";
    auto res = system(comMap.c_str());
    if (res != 0){
        cout << "ERROR minimap2 failed tt, while running " << comMap << endl;
        exit(1);
    }
    auto t3 = high_resolution_clock::now();

    // string comSamtools = SAMTOOLS+" view -F0x900 " + outputFolder + "dbg_"+id+".bam 2>" + outputFolder + "trash.txt | "
    //                 +wtdbg2_folder+"/wtpoa-cns -d " + outputFolder + "dbg_"+id+".raw.fa -i - -fo " + outputFolder + "dbg_"+id+".cns.fa 2>" + outputFolder + "trash.txt";
    string comSamtools = SAMTOOLS+" view -F0x900 " + outputFolder + "dbg_"+id+".bam 2>" + outputFolder + "trash.txt | "
                    +wtdbg2_folder+"/wtpoa-cns -t 1 -d " + ref + " -i - -fo " + outputFolder + "dbg_"+id+".cns.fa 2>" + outputFolder + "trash.txt";
    
    res = system(comSamtools.c_str());
    if (res != 0){
        cout << "ERROR samtools failed, while running " << comSamtools << endl;
        exit(1);
    }
    // cout << "Samtools command dici " << comSamtools << endl;
    auto t4 = high_resolution_clock::now();

    string new_contigs_file = outputFolder + "wtdbg2_"+id+".fa";
    string comUnfold = "awk '{if(\">\" == substr($1,1,1)){ printf \"\\n\"; print;} else printf $1;}' " + outputFolder + "dbg_"+id+".cns.fa >" + new_contigs_file + " 2>" + outputFolder + "trash.txt";
    res = system(comUnfold.c_str());
    if (res != 0){
        cout << "ERROR awk failed, while running " << comUnfold << endl;
        exit(1);
    }
    auto t5 = high_resolution_clock::now();

    //rename the contigs of the new assembly to be sure not to conflict with the original contigs
    string prefix = "reassembled_by_HairSplitter_";
    rename_reads(new_contigs_file, prefix);

    auto t6 = high_resolution_clock::now();

    // cout << "THe times are ggtvy : " <<  duration_cast<milliseconds>(t1-t0).count() << " " <<  duration_cast<milliseconds>(t2-t1).count() <<
    //     " " <<  duration_cast<milliseconds>(t3-t2).count() << " " <<  duration_cast<milliseconds>(t4-t3).count() << " " <<  duration_cast<milliseconds>(t5-t4).count() <<
    //     " " <<  duration_cast<milliseconds>(t6-t5).count() << endl;
}







