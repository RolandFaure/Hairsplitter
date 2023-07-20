#include "tools.h"
#include "reassemble_unaligned_reads.h"
#include "edlib/include/edlib.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::stoi;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::min;
using std::max;

//input : a CIGAR
//output : an alignment string where 1 letter = 1 base. i.e. 5M1D1M -> MMMMMIM
std::string convert_cigar(std::string &cigar){

    string res;
    string num = "";
    for (auto c : cigar){
        if ((int)c - '0' >= 0 && (int)c - '0' <= 9){
            num += c;
        }
        else{
            int n = std::stoi(num);
            for (auto i = 0 ; i < n ; i++){
                res += c;
            }
            num = "";
        }
    }

    return res;
}

string reverse_complement(string &seq){
    string res = "";
    for (auto c : seq){
        if (c == 'A'){
            res = 'T' + res;
        }
        else if (c == 'T'){
            res = 'A' + res;
        }
        else if (c == 'C'){
            res = 'G' + res;
        }
        else if (c == 'G'){
            res = 'C' + res;
        }
        else{
            res = c + res;
        }
    }
    return res;
}

void print_alignment(std::string &ref, std::string &read, std::string &cigar, int start, int end){

    string ref_aligned;
    string read_aligned;

    int ref_pos = start;
    int read_pos = 0;

    for (auto i = 0 ; i < start ; i++){
        ref_aligned += '-';
    }

    for (auto i = 0 ; i < cigar.size() ; i++){
        if (cigar[i] == 'M' || cigar[i] == 'X' || cigar[i] == '='){
            ref_aligned += ref[ref_pos];
            read_aligned += read[read_pos];
            ref_pos++;
            read_pos++;
        }
        else if (cigar[i] == 'D'){
            ref_aligned += '-';
            read_aligned += read[read_pos];
            read_pos++;
        }
        else if (cigar[i] == 'I'){
            ref_aligned += ref[ref_pos];
            read_aligned += '-';
            ref_pos++;
        }
    }

    // for (auto i = ref_pos ; i < end-1 ; i++){
    //     ref_aligned += '-';
    // }

    cout << ref_aligned << endl;
    cout << read_aligned << endl;

}

int compute_edit_distance(std::string &cigar, std::string &ref, std::string &read, int start, int end){
    
        int res = 0;
    
        int ref_pos = start;
        int read_pos = 0;
    
        for (auto i = 0 ; i < cigar.size() ; i++){
            if (cigar[i] == 'M' || cigar[i] == 'X' || cigar[i] == '='){
                if (ref[ref_pos] != read[read_pos]){
                    res++;
                }
                ref_pos++;
                read_pos++;
            }
            else if (cigar[i] == 'D'){
                res++;
                read_pos++;
            }
            else if (cigar[i] == 'I'){
                res++;
                ref_pos++;
            }
        }
    
        return res;
}

/**
 * @brief Convert a GFA file to a FASTA file
 * 
 * @param gfa_file 
 * @param fasta_file 
 */
void convert_GFA_to_FASTA(std::string &gfa_file, std::string &fasta_file){
    ifstream gfa(gfa_file);
    ofstream fasta(fasta_file);

    string line;
    while (getline(gfa, line)){
        if (line[0] == 'S'){

            //the name is the second field
            string name = line.substr(2, line.find('\t', 2) - 2);
            name = name.substr(0, name.find(' '));

            //the sequence is the third field
            string seq = line.substr(line.find('\t', 2) + 1, line.find('\t', line.find('\t', 2) + 1) - line.find('\t', 2) - 1);
            if (seq.size() == 0){
                continue;
            }
            fasta << ">" << name << endl;
            fasta << seq << endl;
        }
    }
    gfa.close();
    fasta.close();
}

/**
 * @brief Convert a FASTA file to a GFA file
 * 
 * @param fasta_file 
 * @param gfa_file 
 */
void convert_FASTA_to_GFA(std::string &fasta_file, std::string &gfa_file){
    ifstream fasta(fasta_file);
    if (!fasta.is_open()){
        cout << "ERROR : could not open " << fasta_file << endl;
        exit(1);
    }
    ofstream gfa(gfa_file);
    if (!gfa.is_open()){
        cout << "ERROR : could not open " << gfa_file << endl;
        exit(1);
    }

    string line;
    string seq;
    string name;
    int i = 0;
    while (getline(fasta, line)){
        if (line[0] == '>'){
            if (i != 0){
                gfa << "S\t" << name << "\t" << seq << endl;
            }
            string name = line.substr(1);
        }
        else{
            seq += line;
        }
        i++;
    }
    if (i != 0){
        gfa << "S\t" << name << "\t" << seq << endl;
    }
    fasta.close();
    gfa.close();
}

/**
 * @brief Function that uses racon to polish a sequence
 * 
 * @param backbone sequence to be polished
 * @param polishingReads list of reads to polish it
 * @param overhang length of the unpolished ends to be used
 * @param id an id (typically a thread id) to be sure intermediate files do not get mixed up with other threads 
 * @param techno the sequencing technology used to generate the reads (ont, pacbio, hifi)
 * @param MINIMAP path to the minimap2 executable
 * @param RACON path to the racon executable
 * @return polished sequence 
 */
string consensus_reads(string const &backbone, vector <string> &polishingReads, string &id, string &outFolder, string& techno, string &MINIMAP, string &RACON){
    
    if (polishingReads.size() == 0){
        return backbone;
    }

    //check if the last character of outFolder is a /
    if (outFolder[outFolder.size()-1] != '/'){
        outFolder += "/";
    }

    system("mkdir tmp/ 2> trash.txt");
    std::ofstream outseq(outFolder+"unpolished_"+id+".fasta");
    outseq << ">seq\n" << backbone;
    outseq.close();

    std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        polishseqs << ">read"+std::to_string(read)+"\n" << polishingReads[read] << "\n";
    }
    polishseqs.close();

    // assemble de novo using wtdbg2
    /*
    string comAsm = "/home/rfaure/Documents/software/wtdbg2/wtdbg2 -e 5 -l 1000 -L 3000 -S 1 -R -o tmp/wtdbg2"+id+" -i tmp/reads_"+id+".fasta 2>tmp/trash.txt";
    system(comAsm.c_str());

    string cons_wtdbg2 = "/home/rfaure/Documents/software/wtdbg2/wtpoa-cns -i tmp/wtdbg2"+id+".ctg.lay.gz -fo tmp/dbg.raw"+id+".fa 2>tmp/trash.txt";
    int res_wtdbg2 = system(cons_wtdbg2.c_str());

    if (res_wtdbg2){
        // polish consensus, not necessary if you want to polish the assemblies using other tools
        string comMap = "minimap2 -ax map-pb -r2k tmp/dbg.raw"+id+".fa tmp/reads_"+id+".fasta 2>tmp/trash.txt | samtools sort >tmp/dbg"+id+".bam 2>tmp/trash.txt";
        system(comMap.c_str());

        string comSamtools = "samtools view -F0x900 tmp/dbg"+id+".bam 2>tmp/trash.txt | /home/rfaure/Documents/software/wtdbg2/wtpoa-cns -d tmp/dbg.raw"+id+".fa -i - -fo tmp/dbg.cns"+id+".fa 2>tmp/trash.txt";
        system(comSamtools.c_str());

        string comUnfold = "awk '{if(\">\" == substr($1,1,1)){ printf \"\\n\"; print;} else printf $1;}' tmp/dbg.cns"+id+".fa > tmp/consensus"+id+".fa  2>tmp/trash.txt";
        system(comUnfold.c_str());

        //now read the consensus and return it
        std::ifstream in("tmp/consensus"+id+".fa");
        string consensus;
        string line;
        while (std::getline(in, line)){
            if (line[0] != '>'){
                consensus += line;
            }
        }
        in.close();
        return consensus;
    }

    else{
    */

    string technoFlag;
    if (techno == "ont"){
        technoFlag = " -x map-ont ";
    }
    else if (techno == "pacbio"){
        technoFlag = " -x map-pb ";
    }
    else if (techno == "hifi"){
        technoFlag = " -x map-hifi ";
    }

    string com = " -t 1 "+ technoFlag + " " + outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".paf 2>"+ outFolder +"trash.txt";
    string commandMap = MINIMAP + com; 
    system(commandMap.c_str());

    com = " -w 500 -e 1 -t 1 "+ outFolder +"reads_"+id+".fasta "+ outFolder +"mapped_"+id+".paf "+ outFolder +"unpolished_"+id+".fasta > "+ outFolder +"polished_"+id+".fasta 2>"+ outFolder +"trash.txt";
    string commandPolish = RACON + com;
    system(commandPolish.c_str());

    //polish using CONSENT
    // string CONSENT = "/home/rfaure/Documents/software/CONSENT/CONSENT-polish";
    // string com = " -S 30 --contigs tmp/unpolished_"+id+".fasta --reads tmp/reads_"+id+".fasta --out tmp/polished_"+id+".paf >tmp/log_CONSENT.txt 2>tmp/trash.txt";
    // string commandPolish = CONSENT + com;
    // system(commandPolish.c_str());

    std::ifstream polishedRead(outFolder +"polished_"+id+".fasta");
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
    auto before_size = min(size_t(300), backbone.size());
    auto after_size = min(size_t(200), consensus.size());

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

    result = edlibAlign(after_end.c_str(), after_end.size(), before_end.c_str(), before_end.size(),
                                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    int end_pos = result.endLocations[0]+1;
    string additional_end_seq = before_end.substr(end_pos , before_end.size()-end_pos);
    edlibFreeAlignResult(result);
    
    // cout << "consensus : " << endl << (additional_start_seq + consensus + additional_end_seq).substr(100,150) << endl;
    // cout << backbone.substr(100,150) << endl;
    
    return additional_start_seq + consensus + additional_end_seq;

}

/**
 * @brief renames all the reads of the fasta file adding the prefix
 * 
 * @param fasta_file 
 * @param prefix 
 */
void rename_reads(std::string &fasta_file, std::string &prefix){
    std::ifstream in(fasta_file);
    std::ofstream out(fasta_file+".tmp");
    std::string line;
    int i = 0;
    while (std::getline(in, line)){
        if (line[0] == '>'){
            out << ">"+prefix+std::to_string(i) << std::endl;
            i++;
        }
        else{
            out << line << std::endl;
        }
    }
    in.close();
    out.close();
    std::remove(fasta_file.c_str());
    std::rename((fasta_file+".tmp").c_str(), fasta_file.c_str());
}

/**
 * @brief assebmles the reads using wtdbg2
 * 
 * @param backbone 
 * @param polishingReads 
 * @param id 
 * @param outFolder 
 * @param techno 
 * @return std::string 
 */
std::string consensus_reads_wtdbg2(
    std::string const &backbone, 
    std::vector <std::string> &polishingReads, 
    std::string &id,
    std::string &outFolder,
    std::string &techno,
    string &MINIMAP, 
    string &RACON
){

    //aligns all the reads on backbone to assemble only the interesting parts
    if (polishingReads.size() == 0){
        return backbone;
    }

    system("mkdir tmp/ 2> trash.txt");
    std::ofstream outseq(outFolder+"unpolished_"+id+".fasta");
    outseq << ">seq\n" << backbone;
    outseq.close();

    std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        polishseqs << ">"+std::to_string(read)+"\n" << polishingReads[read] << "\n";
    }
    polishseqs.close();

    string technoFlag = "";
    if (techno == "pacbio"){
        technoFlag = " -x map-pb ";
    }
    else if (techno == "nanopore"){
        technoFlag = " -x map-ont ";
    }
    else if (techno == "hifi"){
        technoFlag = " -x map-hifi ";
    }
    else{
        cout << "unknown sequencing technology" << endl;
        exit(1);
    }

    string com = " -t 1 "+ technoFlag = " "+ outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".paf 2>"+ outFolder +"trash.txt";
    string commandMap = MINIMAP + com; 
    system(commandMap.c_str());

    vector<string> clippedReads;

    //parse the paf file and clip all the reads
    ifstream in(outFolder +"mapped_"+id+".paf");
    string line;
    while(getline(in, line)){

        int read_start = 0;
        int read_end = 0;
        int contig_start = 0;
        int contig_end = 0;
        int contig_length = 0;
        int read = -1;
        string orientation = "+";

        short fieldnumber = 0;
        string field;
        std::istringstream line2(line);
        while (getline(line2, field, '\t'))
        {
            if (fieldnumber == 0){
                read = std::atoi(field.c_str());
            }
            if (fieldnumber == 2){
                read_start = std::atoi(field.c_str());
            }
            if (fieldnumber == 3){
                read_end = std::atoi(field.c_str());
            }
            if (fieldnumber == 4){
                orientation = field;
            }
            if (fieldnumber == 6){
                contig_length = std::atoi(field.c_str());
            }
            if (fieldnumber == 7){
                contig_start = std::atoi(field.c_str());
            }
            if (fieldnumber == 8){
                contig_end = std::atoi(field.c_str());
            }
            fieldnumber+= 1;
        }

        int overhangLeft = contig_start;
        int overhangRight = contig_length - contig_end;
        if (orientation == "-"){
            overhangLeft = contig_length - contig_end;
            overhangRight = contig_start;
        }
        int left = max(0,read_start-overhangLeft);
        int length_of_subread = min(int(polishingReads[read].size())-left, read_end+overhangRight-left);
        clippedReads.push_back(polishingReads[read].substr(left , length_of_subread));
    }

    std::string newcontig = "";
    //create a file containing the reads
    string fileReads = outFolder+"reads_"+id+".fa";
    std::ofstream o(fileReads);
    for (auto r = 0 ; r < clippedReads.size() ; r++){
        o << ">" << r <<"\n" << clippedReads[r] << "\n";
    }
    o.close();

    string unpolished = outFolder +"unpolished_"+id+".fasta";
    assemble_with_wtdbg2(fileReads, outFolder, unpolished, id);
    //get the result
    string resultFile = outFolder+"/wtdbg2_"+id+".fa";
    std::ifstream assembled(resultFile);
    if (!assembled){
        cout << "problem reading files in modify_gfa dlfoc, while trying to read " << resultFile << endl;
        throw std::invalid_argument( "File could not be read" );
    }
    line = "";
    bool nextLine = false;
    while(getline(assembled, line)){
        if (line[0] != '>' && nextLine){
            newcontig = line; //there should be only one contig, which is what was just assembled
            break;
        }
        else if (line[0] == '>'){
            nextLine = true;
        }
    }

    return newcontig;
}



