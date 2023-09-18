#include "tools.h"
// #include "reassemble_unaligned_reads.h"
#include "edlib/include/edlib.h"

#include <iostream>
#include <fstream>
#include <algorithm>

using std::cout;
using std::endl;
using std::stoi;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::min;
using std::max;
using std::pair;

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

//output : a CIGAR
//input : an alignment string where 1 letter = 1 base. i.e. 5M1D1M -> MMMMMIM
std::string convert_cigar2(std::string &cigar){

    string res;
    int number = 0;
    char current = ' ';
    for (auto c : cigar){
        if (c == current || current == ' '){
            number++;
            current = c;
        }
        else{
            res += std::to_string(number) + current;
            number = 1;
            current = c;
        }
    }
    res += std::to_string(number) + current;

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
            //split line on the first space
            name = line.substr(1, line.find(' ')-1);
            if (i != 0){
                gfa << "S\t" << name << "\t" << seq << endl;
            }
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
 * @param full_backbone the full backbone, to be used if the reads do not align well on the backbone
 * @param start_pos_on_full_backbone the position of the backbone on the full backbone
 * @param sizeOfWindow the size of the window of the backbones
 * @param polishingReads list of reads to polish it
 * @param fullReads the full reads, to be used if the reads do not align well on the backbone
 * @param CIGARs list of CIGARs of the alignment of the read on the backbone
 * @param overhang length of the unpolished ends to be used
 * @param id an id (typically a thread id) to be sure intermediate files do not get mixed up with other threads 
 * @param techno the sequencing technology used to generate the reads (ont, pacbio, hifi)
 * @param MINIMAP path to the minimap2 executable
 * @param RACON path to the racon executable
 * @return polished sequence 
 */
string consensus_reads(
    string &backbone, 
    string &full_backbone, 
    int start_pos_on_full_backbone, 
    int sizeOfWindow, 
    vector <string> &polishingReads,
    vector <string> &fullReads,
    vector <pair<string,int>> &CIGARs,
    string &id, 
    string &outFolder, 
    string& techno, 
    string &MINIMAP, 
    string &RACON,
    string &path_to_python,
    std::string &path_src){
    
    if (polishingReads.size() == 0){
        return backbone;
    }

    //check if the last character of outFolder is a /
    if (outFolder[outFolder.size()-1] != '/'){
        outFolder += "/";
    }

    if (polishingReads.size() == 0){
        return backbone;
    }

    std::ofstream outseq(outFolder+"unpolished_"+id+".fasta");
    outseq << ">seq\n" << backbone << endl;
    outseq.close();

    std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        if (polishingReads[read].size() > 100){
            polishseqs << ">read"+std::to_string(read)+"\n" << polishingReads[read] << "\n";
        }
    }
    polishseqs.close();


    string technoFlag;
    if (techno == "ont" && false){ //minimap with setting ont can be absolutely terrible in the case of a long homopolymer
        technoFlag = " -x map-ont ";
    }
    else if (techno == "pacbio" || true){
        technoFlag = " -x map-pb ";
    }
    else if (techno == "hifi"){
        technoFlag = " -x map-hifi ";
    }

    // string com = " -t 1 "+ technoFlag + " " + outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".paf 2>"+ outFolder +"trash.txt";
    // string commandMap = MINIMAP + com; 
    // auto map = system(commandMap.c_str());
    // if (map != 0){
    //     cout << "ERROR minimap2 failed, while running " << commandMap << endl;
    //     exit(1);
    // }

    //create a sam file from the CIGARs and the reads
    std::ofstream sam(outFolder+"mapped_"+id+".sam");
    sam << "@HD\tVN:1.6\tSO:coordinate" << endl;
    sam << "@SQ\tSN:seq\tLN:" << backbone.size() << endl;
    for (int read = 0 ; read < polishingReads.size() ; read++){
        if (polishingReads[read].size() > 100){
            sam << "read" << read << "\t0\tseq\t"<< CIGARs[read].second <<"\t60\t" << CIGARs[read].first << "\t*\t0\t0\t" << polishingReads[read] << "\t*\tAS:i:0\tXS:i:0" << endl;
        }   
    }
    sam.close();

    //sort and index mapped_id.sam
    string command = "samtools sort "+ outFolder +"mapped_"+id+".sam > "+ outFolder +"mapped_"+id+".bam && samtools index "+ outFolder +"mapped_"+id+".bam";
    auto sort = system(command.c_str());
    if (sort != 0){
        cout << "ERROR samtools sort failed, while running " << command << endl;
        exit(1);
    }

    //run a basic consensus
    //use quasitools : "quasitools consensus -p 100 hs/tmp/mapped_0.bam hs/tmp/unpolished_0.fasta > hs/tmp/consensus_0.fasta && mv hs/tmp/consensus_0.fasta hs/tmp/polished_0.fasta"
    command =  "samtools consensus "+ outFolder +"mapped_"+id+".bam > "+ outFolder +"consensus_"+id+".fasta";
    // cout << "Running " << command << endl;
    auto res_cons = system(command.c_str());
    if (res_cons != 0){
        cout << "ERROR basic_consensus failed, while running " << command << endl;
        exit(1);
    }

    //then map all the reads on unpolsihed.fasta to obtain a new mapped.sam
    string com = " -a -t 1 "+ technoFlag + " " + outFolder +"consensus_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".sam 2>"+ outFolder +"trash.txt";
    command = MINIMAP + com;
    auto map2 = system(command.c_str());
    if (map2 != 0){
        cout << "ERROR minimap2 failed, while running " << command << endl;
        exit(1);
    }

    //check that the alignment of the reads were correct, otherwise change the backbone
    string nameOfFile = outFolder +"mapped_"+id+".sam";
    bool bbaligns = true;
    int alternativeBackbone = 0;
    while (!check_alignment(nameOfFile) && alternativeBackbone < polishingReads.size()){ //this means that no reads aligned really well in the first 5 reads
        // cout << "in toolsl qfhggfqsft taking another ref " << endl;
        bbaligns = false;
        //realign taking the first read as a backbone
        backbone = polishingReads[alternativeBackbone];
        alternativeBackbone++;

        std::ofstream outseq(outFolder+"consensus_"+id+".fasta");
        outseq << ">seq\n" << backbone;
        outseq.close();

        std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
        for (int read =0 ; read < fullReads.size() ; read++){
            polishseqs << ">read"+std::to_string(read)+"\n" << fullReads[read] << "\n";
        }
        polishseqs.close();

        string com = " -a -t 1 "+ technoFlag + " " + outFolder +"consensus_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".sam 2>"+ outFolder +"trash.txt";
        string commandMap = MINIMAP + com;
        auto map = system(commandMap.c_str());
        if (map != 0){
            cout << "ERROR minimap2 failed, while running " << commandMap << endl;
            exit(1);
        }
        outseq.close();

    }

    if (alternativeBackbone == polishingReads.size()){ //the group of reads is really weird and we cannot align them on anything coherently
        return backbone;
    }

    // cout << "minimap2 done in tools , ran command " << commandMap << endl;

    com = " -w 500 -e 1 -t 1 "+ outFolder +"reads_"+id+".fasta "+ outFolder +"mapped_"+id+".sam "+ outFolder +"consensus_"+id+".fasta > "+ outFolder +"polished_"+id+".fasta 2>"+ outFolder +"trash.txt";
    string commandPolish = RACON + com;
    auto polishres = system(commandPolish.c_str());
    if (polishres != 0){
        cout << "ERROR racon failed, while running " << commandPolish << endl;
        exit(1);
    }

    std::ifstream polishedRead(outFolder +"polished_"+id+".fasta");
    string consensus;
    string line;
    while(getline(polishedRead, line)){
        if (line[0] != '>'){
            consensus = line; 
        }
    }

    if (consensus == ""){ //that's typically if there are too few reads to consensus
        return backbone;
    }

    string new_seq;
    if (bbaligns){ //in the case the reads aligned well on the backbone
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

        new_seq = additional_start_seq + consensus + additional_end_seq;

    }
    else{ // in the case we had to reassemble, it is more tricky to find the limits of the seq, so it makes no sense to align the borders
        new_seq = consensus;
    }
    
    // cout << "consensus : " << endl << (additional_start_seq + consensus + additional_end_seq).substr(100,150) << endl;
    // cout << backbone.substr(100,150) << endl;

    //remove all the temporary files
    // std::remove((outFolder+"unpolished_"+id+".fasta").c_str());
    // std::remove((outFolder+"reads_"+id+".fasta").c_str());
    // std::remove((outFolder+"mapped_"+id+".sam").c_str());
    // std::remove((outFolder+"polished_"+id+".fasta").c_str());
    std::remove((outFolder+"trash.txt").c_str());
    
    return new_seq;

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
 * @brief creates a consensus sequence using medaka
 * 
 * @param backbone sequence to be polished
 * @param polishingReads reads to polish it
 * @param id thread id to be sure intermediate files do not get mixed up with other threads
 * @param outFolder tmp folder to store intermediate files
 * @return std::string 
 */
std::string consensus_reads_medaka(
    std::string const &backbone, 
    std::vector <std::string> &polishingReads, 
    std::string &id,
    std::string outFolder,
    std::string &MEDAKA,
    std::string &SAMTOOLS,
    std::string &path_to_python,
    std::string &path_src){

    outFolder += "/";
    //output all polishing reads in a file
    std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        polishseqs << ">"+std::to_string(read)+"\n" << polishingReads[read] << "\n";
    }
    //add the backbone to the reads
    polishseqs << ">seq\n" << backbone;
    polishseqs.close();

    //output the backbone in a file
    std::ofstream outseq(outFolder+"unpolished_"+id+".fasta");
    outseq << ">seq\n" << backbone;
    outseq.close();

    //create a bam file with the reads aligned on the backbone and index it
    string comMap = "minimap2 -ax map-pb -r2k "+ outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta 2>"+ outFolder +"trash.txt "
        +"| "+ SAMTOOLS + " sort >"+ outFolder +"mapped_"+id+".bam && "+ SAMTOOLS + " index "+ outFolder +"mapped_"+id+".bam 2>"+ outFolder +"trash.txt";
    auto map = system(comMap.c_str());
    if (map != 0){
        cout << "ERROR minimap2 yuq fd failed, while running " << comMap << endl;
        exit(1);
    }

    //now create a first rough polish using a basic pileup
    // Build the command to run haplodmf_count_freqs.py
    string command;
    command =  path_to_python + " " + path_src + "/haplodmf_count_freqs.py " + outFolder +"mapped_"+id+".bam " + outFolder +"consensus_"+id+".fasta 2>"+ outFolder +"trash.txt";    
    // cout << "Running " << command << endl;
    
    // Run the command
    int return_code = system(command.c_str());
    if (return_code != 0) {
        cout << "Error running command: " << command << endl;
        exit(1);
    }

    //run medaka on the consensus
    string com = MEDAKA + "_consensus -i "+ outFolder +"reads_"+id+".fasta -d "+ outFolder +"consensus_"+id+".fasta -o "+ outFolder +"medaka_"+id+" -t 1 -f -x 2>"+ outFolder +"trash.txt >" +outFolder +"trash.txt" ;
    // cout << "Running in tools.cpp yyxk" << com << endl;
    auto med_res = system(com.c_str());
    if (med_res != 0){        
        cout << "ERROR medaka failed, while running " << com << endl;
        exit(1);
    }

    //get the result
    std::ifstream consensus(outFolder +"medaka_"+id+"/consensus.fasta");
    if (!consensus){
        cout << "problem reading files in tools.cpp wiwio, while trying to read " << outFolder +"medaka_"+id+"/consensus.fasta" << endl;
        throw std::invalid_argument( "File could not be read" );
    }

    bool nextLine = false;
    string newcontig = "";
    string line;
    while(getline(consensus, line)){
        if (line[0] != '>' && nextLine){
            newcontig = line; //there should be only one contig, which is what was just assembled
            break;
        }
        else if (line[0] == '>'){
            nextLine = true;
        }
    }

    consensus.close();

    //remove all the temporary files
    std::remove((outFolder+"unpolished_"+id+".fasta").c_str());
    std::remove((outFolder+"reads_"+id+".fasta").c_str());
    //remove the medaka folder
    command = "rm -r "+ outFolder +"medaka_"+id;
    auto rm = system(command.c_str());
    if (rm != 0){
        cout << "ERROR rm failed, while running " << command << endl;
        exit(1);
    }

    // cout << "medaka done in tools.cpp\n" << newcontig << endl;

    return newcontig;
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
    string &RACON,
    string &SAMTOOLS,
    string &WTDBG2
){

    outFolder += "/";
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
    else if (techno == "ont"){
        technoFlag = " -x map-ont ";
    }
    else if (techno == "hifi"){
        technoFlag = " -x map-hifi ";
    }
    else{
        cout << "unknown sequencing technology: " << techno << endl;
        exit(1);
    }

    string com = " -t 1 "+ technoFlag + " "+ outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".paf 2>"+ outFolder +"trash.txt";
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
    assemble_with_wtdbg2(fileReads, outFolder, unpolished, id, WTDBG2, MINIMAP, SAMTOOLS);
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


/**
 * @brief Assemble and a polish a file of reads using wtdbg2
 * 
 * @param fileReads All the reads to assemble
 * @param outputFolder 
 * @param ref File to the reference
 * @param num_threads 
 */
void assemble_with_wtdbg2(std::string &fileReads, std::string outputFolder, std::string &ref, std::string id, std::string &WTDBG2, std::string &MINIMAP, std::string &SAMTOOLS){

    auto lastslash = WTDBG2.find_last_of("/");
    if (lastslash == string::npos){
        lastslash = 0;
    }
    string wtdbg2_folder = WTDBG2.substr(0, lastslash);

    if (ref == ""){
        string comAsm = wtdbg2_folder+ "/wtdbg2 -A -e 1 -l 200 -L 0 -S 100 --no-read-clip --no-chainning-clip --ctg-min-length 200 --ctg-min-nodes 0 -R -o "
                                + outputFolder + "wtdbg2_"+id+" -i " + fileReads + " 2>"+outputFolder+"trash.txt";
        auto res = system(comAsm.c_str());
        if (res != 0){
            cout << "ERROR wtdbg2 failed, while running " << comAsm << endl;
            exit(1);
        }

        string cons_wtdbg2 = wtdbg2_folder+"/wtpoa-cns -t 1 -i " + outputFolder + "wtdbg2_"+id+".ctg.lay.gz -fo " + outputFolder + "dbg_"+id+".raw.fa 2>tmp/trash.txt";
        int res_wtdbg2 = system(cons_wtdbg2.c_str());
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

    string new_contigs_file = outputFolder + "wtdbg2_"+id+".fa";
    string comUnfold = "awk '{if(\">\" == substr($1,1,1)){ printf \"\\n\"; print;} else printf $1;}' " + outputFolder + "dbg_"+id+".cns.fa >" + new_contigs_file + " 2>" + outputFolder + "trash.txt";
    res = system(comUnfold.c_str());
    if (res != 0){
        cout << "ERROR awk failed, while running " << comUnfold << endl;
        exit(1);
    }

    //rename the contigs of the new assembly to be sure not to conflict with the original contigs
    string prefix = "reassembled_by_HairSplitter_";
    rename_reads(new_contigs_file, prefix);

    // cout << "THe times are ggtvy : " <<  duration_cast<milliseconds>(t1-t0).count() << " " <<  duration_cast<milliseconds>(t2-t1).count() <<
    //     " " <<  duration_cast<milliseconds>(t3-t2).count() << " " <<  duration_cast<milliseconds>(t4-t3).count() << " " <<  duration_cast<milliseconds>(t5-t4).count() <<
    //     " " <<  duration_cast<milliseconds>(t6-t5).count() << endl;
}

/**
 * @brief Check if the reads aligned well on the backbone
 * 
 * @param paf_file 
 * @return true 
 * @return false 
 */
bool check_alignment(std::string &paf_file){
    ifstream in(paf_file);
    string line;
    int read = 0;
    int aligned = 0;
    //if it is a paf file 
    if (paf_file.substr(paf_file.size()-4,4) == ".paf"){ 
        while (getline(in, line)){
            //look at the 12th field of the line
            read += 1;
            short fieldnumber = 0;
            string field;
            std::istringstream line2(line);
            string quality;
            while (getline(line2, field, '\t'))
            {
                if (fieldnumber == 11){
                    quality = field;
                }
                fieldnumber+= 1;
            }
            if (quality == "60"){
                aligned += 1;
            }
            if (read == 6){
                break;
            }
        }
    }
    //if it is a sam file
    else if (paf_file.substr(paf_file.size()-4,4) == ".sam"){
        while (getline(in, line)){
            //if there are more than 70% M the read is aligned
            read += 1;
            short fieldnumber = 0;
            string field;
            std::istringstream line2(line);
            string cigar;
            while (getline(line2, field, '\t'))
            {
                if (fieldnumber == 5){
                    cigar = field;
                }
                fieldnumber+= 1;
            }
            int M = 0;
            string number = "";
            int cigar_size = 0;
            for (auto c : cigar){
                if (c == 'M'){
                    M += std::atoi(number.c_str());
                    cigar_size += std::atoi(number.c_str());
                    number = "";
                }
                else if (c == 'I' || c == 'D'){
                    cigar_size += std::atoi(number.c_str());
                    number = "";
                }
                else{
                    number += c;
                }
            }
            if (M > 0.7*cigar_size){
                aligned += 1;
            }
        }
    }

    return (aligned >= 2);
}

