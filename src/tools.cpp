#include "tools.h"
// #include "reassemble_unaligned_reads.h"
#include "edlib/include/edlib.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <unordered_set>

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
using std::unordered_map;

//input : a CIGAR
//output : an alignment string where 1 letter = 1 base. i.e. 5M1D1M -> MMMMMIM
std::string convert_cigar(std::string &cigar){

    if (cigar == "*"){
        return "";
    }

    string res;
    string num = "";
    for (auto c : cigar){
        if ((int)c - '0' >= 0 && (int)c - '0' <= 9){
            num += c;
        }
        else{
            int n = 0;
            try {
                n = std::stoi(num);
            }
            catch (...){
                cout << "ERROR : could not convert " << cigar << " to int" << endl;
                exit(1);
            }
            
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

void cut_GFA(const std::string& input_assembly, const std::string& output_assembly, int max_length) {
    std::unordered_map<std::string, int> length_of_contigs;
    std::vector<std::string> L_lines;

    std::ifstream in(input_assembly);
    std::ofstream out(output_assembly);

    if (!in || !out) {
        std::cerr << "Error opening file\n";
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line[0] =='S') {
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, '\t')) {
                tokens.push_back(token);
            }
            length_of_contigs[tokens[1]] = tokens[2].length();
            for (int chunk = 0; chunk * max_length < tokens[2].length(); ++chunk) {
                int end = std::min((chunk + 1) * max_length, (int) tokens[2].length());
                out << "S\t" << tokens[1] << "@" << chunk << "\t" << tokens[2].substr(chunk * max_length, end - chunk * max_length) << "\t";
                for (size_t i = 3; i < tokens.size(); ++i) {
                    out << tokens[i] << "\t";
                }
                out << "\n";
                if (chunk > 0) {
                    out << "L\t" << tokens[1] << "@" << chunk - 1 << "\t+\t" << tokens[1] << "@" << chunk << "\t+\t0M\n";
                }
            }
        } else if (line[0] == 'L') {
            L_lines.push_back(line);
        }
    }

    for (const auto& line : L_lines) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, '\t')) {
            tokens.push_back(token);
        }
        int index1 = (length_of_contigs[tokens[1]] - 1) / max_length;
        int index2 = (length_of_contigs[tokens[3]] - 1) / max_length;
        out << "L\t" << tokens[1] << "@" << index1 << "\t+\t" << tokens[3] << "@" << index2 << "\t+\t" << tokens[5] << "\n";
    }

    in.close();
    out.close();
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
 * @param racon_window_size the size of the window to be used in racon, if -1, then the size of the backbone, if 0 then defaults back to 500
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
    int racon_window_size,
    string &MINIMAP, 
    string &RACON,
    string &path_to_python,
    std::string &path_src){
    
    if (polishingReads.size() == 0){
        return backbone;
    }

    int window_size = racon_window_size;
    if (window_size == -1){
        window_size = backbone.size();
    }
    else if (window_size == 0){
        window_size = 500;
    }

    //check if the last character of outFolder is a /
    if (outFolder[outFolder.size()-1] != '/'){
        outFolder += "/";
    }

    std::ofstream outseq(outFolder+"unpolished_"+id+".fasta");
    outseq << ">seq\n" << backbone << endl;
    outseq.close();

    std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
    for (int read =0 ; read < polishingReads.size() ; read++){
        if (polishingReads[read].size() > 0){
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

    string com2 = " -a -t 1 "+ technoFlag + " " + outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".sam 2>"+ outFolder +"trash.txt";
    string commandMap = MINIMAP + com2; 
    auto map = system(commandMap.c_str());
    if (map != 0){
        cout << "ERROR minimap2 failed, while running " << commandMap << endl;
        exit(1);
    }

    //create a sam file from the CIGARs and the reads
    // std::ofstream sam(outFolder+"mapped_"+id+".sam");
    // sam << "@HD\tVN:1.6\tSO:coordinate" << endl;
    // sam << "@SQ\tSN:seq\tLN:" << backbone.size() << endl;
    // for (int read = 0 ; read < polishingReads.size() ; read++){
    //     // cout << "readeue " << read << " " << polishingReads[read].size() << endl;
    //     if (polishingReads[read].size() > 0){
    //         sam << "read" << read << "\t0\tseq\t"<< CIGARs[read].second <<"\t60\t" << CIGARs[read].first << "\t*\t0\t0\t" << polishingReads[read] << "\t*\tAS:i:0\tXS:i:0" << endl;
    //     }   
    // }
    // sam.close();

    //check that the alignment of the reads were correct, otherwise change the backbone
    string nameOfFile = outFolder +"mapped_"+id+".sam";
    int good_aln = check_alignment(nameOfFile);
    if (good_aln != 0){ //this means that no reads aligned really well, hence recompute a backbone

        string newbackbone;
        if (good_aln == 1){ //only small indels
            // cout << "small indelsicine" << endl;
            std::ofstream polishseqs(outFolder+"reads_"+id+".fasta");
            for (int read = 0 ; read < fullReads.size() ; read++){
                polishseqs << ">read"+std::to_string(read)+"\n" << fullReads[read] << "\n";
            }
            polishseqs.close();
            string com = " -a --secondary=no -t 1 "+ technoFlag + " " + outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".sam 2>"+ outFolder +"trash.txt";
            string commandMap = MINIMAP + com;
            auto map = system(commandMap.c_str());
            if (map != 0){
                cout << "ERROR minimap2 failed, while running " << commandMap << endl;
                exit(1);
            }

            newbackbone = alternative_backbone(nameOfFile, backbone);
        }
        else{ //bib problem in aln so reassemble
            newbackbone = basic_assembly(outFolder+"reads_"+id+".fasta", MINIMAP, outFolder, id);
        }

        if (newbackbone == ""){
            newbackbone  = backbone;
            if (newbackbone == "") { //very bizarre
                return backbone;
            }
        }

        // cout << "length of new backbone: " << newbackbone.size() << endl;

        std::ofstream outseq(outFolder+"unpolished_"+id+".fasta");
        outseq << ">seq\n" << newbackbone << endl;
        outseq.close();

        //realign all the reads on the new backbone
        string com2 = " -a --secondary=no -t 1 "+ technoFlag + " " + outFolder +"unpolished_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".sam 2>"+ outFolder +"trash.txt";
        string commandMap2 = MINIMAP + com2;
        auto map2 = system(commandMap2.c_str());
        if (map2 != 0){
            cout << "ERROR minimap2 failed, while running " << commandMap2 << endl;
            exit(1);
        }
        outseq.close();

    }

    //sort and index mapped_id.sam
    string command = "samtools sort "+ outFolder +"mapped_"+id+".sam > "+ outFolder +"mapped_"+id+".bam && samtools index "+ outFolder +"mapped_"+id+".bam";
    auto sort = system(command.c_str());
    if (sort != 0){
        cout << "ERROR samtools sort failed, while running " << command << endl;
        exit(1);
    }

    //run a basic consensus
    command = "samtools consensus -m simple -c 0 "+ outFolder +"mapped_"+id+".bam > "+ outFolder +"consensus_"+id+".fasta";
    auto res_cons = system(command.c_str());
    if (res_cons != 0){
        cout << "ERROR basic_consensus failed, while running " << command << endl;
        exit(1);
    }

    //then map all the reads on consensus.fasta to obtain a new mapped.sam
    string com = " -a --secondary=no -t 1 "+ technoFlag + " " + outFolder +"consensus_"+id+".fasta "+ outFolder +"reads_"+id+".fasta > "+ outFolder +"mapped_"+id+".sam 2>"+ outFolder +"trash.txt";
    command = MINIMAP + com;
    auto map2 = system(command.c_str());
    if (map2 != 0){
        cout << "ERROR minimap2 failed, while running " << command << endl;
        exit(1);
    }

    // if (alternativeBackbone == polishingReads.size()){ //the group of reads is really weird and we cannot align them on anything coherently
    //     return backbone;
    // }

    // cout << "minimap2 done in tools , ran command " << commandMap << endl;

    com = " -w " + std::to_string(window_size) + " -e 1 -t 1 "+ outFolder +"reads_"+id+".fasta "+ outFolder +"mapped_"+id+".sam "+ outFolder +"consensus_"+id+".fasta > "+ outFolder +"polished_"+id+".fasta 2>"+ outFolder +"trash.txt";
    string commandPolish = RACON + com;
    auto polishres = system(commandPolish.c_str());
    // cout << "command to polish: " << commandPolish << endl;
    if (polishres != 0){
        // cout << "WARNING racon failed, while running " << commandPolish << endl;
        return backbone;
        // exit(1);
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
    string alnfile = "mapped_"+id+".sam";
    if (check_alignment(alnfile) == 0 || true){ //in the case the reads aligned well on the backbone
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
        // cout << "nerwseq : " << new_seq.size() << endl;

    }
    else{ // in the case we had to reassemble, it is more tricky to find the limits of the seq, so it makes no sense to align the borders
        cout << "THE alignfdj does not check out DEBUG code 1033" << endl;
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
    command = SAMTOOLS + " consensus " + outFolder +"mapped_"+id+".bam > "+ outFolder +"consensus_"+id+".fasta";
    // command =  path_to_python + " " + path_src + "/haplodmf_count_freqs.py " + outFolder +"mapped_"+id+".bam " + outFolder +"acgt_"+id+".txt 2>"+ outFolder +"trash.txt";    
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
 * @return 0 (good alignment), 1 (few indels), 2 (breakpoints with S and H)
 */
int check_alignment(std::string &paf_file){
    
    int result_code = 0;
    
    ifstream in(paf_file);
    string line;
    int read = 0;
    int aligned = 0;
    string new_reference = "";
    int nb_reads = 0;
    std::unordered_map<int,int> putative_breakpoints;
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

            if (line[0] == '@'){
                continue;
            }
            read += 1;
            short fieldnumber = 0;
            string field;
            std::istringstream line2(line);
            string cigar;
            bool big_indel = false;
            int startPos;
            string name_of_read;
            while (getline(line2, field, '\t'))
            {
                if (fieldnumber == 0){
                    name_of_read = field;
                }
                if (fieldnumber == 5){
                    cigar = convert_cigar(field);
                }
                if (fieldnumber == 3){
                    startPos = std::atoi(field.c_str());
                }
                fieldnumber+= 1;
            }
            int pos_on_ref = startPos;
            int M = 0;
            int cigar_size = 0;
            int size_of_indel = 0;
            int pos_of_last_breakpoint = -1;
            bool new_indel = true;
            for (auto c : cigar){
                if (c == 'M'){
                    M += 1;
                    cigar_size += 1;
                    size_of_indel = 0;
                    pos_on_ref += 1;
                    new_indel = true;
                }
                else if (c == 'I' || c == 'D'){ //there are supposed to be no clips, as the reads were already tailored to the backbone
                    cigar_size += 1;
                    size_of_indel += 1;
                    if (size_of_indel >= 30 && new_indel){
                        new_indel = false;
                        big_indel = true;
                        int position_indel_approximate = pos_on_ref;
                        if (c == 'I'){
                            position_indel_approximate = int(position_indel_approximate / 100) * 100; //insertions can be positioned arbitrarily, so approximate
                        }
                        if (putative_breakpoints.find(position_indel_approximate) == putative_breakpoints.end()){
                            putative_breakpoints[position_indel_approximate] = 0;
                        }
                        putative_breakpoints[position_indel_approximate] += 1;
                        if (putative_breakpoints[position_indel_approximate] > 2){
                            result_code = 1;
                        }
                    }
                    if (c == 'I'){
                        pos_on_ref += 1;
                    }

                }
                else if (c == 'S' || c == 'H'){ //there are supposed to be no clips, as the reads were already tailored to the backbone
                    size_of_indel += 1;
                    if (size_of_indel > 30){
                        big_indel = true;
                        if (putative_breakpoints.find(pos_on_ref) == putative_breakpoints.end()){
                            putative_breakpoints[pos_on_ref] = 0;
                        }
                        if (pos_on_ref > pos_of_last_breakpoint){
                            putative_breakpoints[pos_on_ref] += 1;
                            pos_of_last_breakpoint = pos_on_ref +1;
                        }
                        if (putative_breakpoints[pos_on_ref] > 2){
                            return 2;
                        }
                    }
                }

            }

            nb_reads++;
            // cout << "fquodfj " << cigar << endl;
        }
    }

    // if (float(aligned)/nb_reads >= 0.4 && float(aligned)/nb_reads < 0.9){
    //     cout << "uarepu tools aligned " << aligned << " " << nb_reads << endl;
    //     exit(1);
    // }
    // cout << "aligned " << aligned << " " << nb_reads << endl;
    if (nb_reads < 2){
        return 2;
    }
    return result_code;
}

/**
 * @brief Compute a new backbone stitching together reads based on the alignment of the reads on the previous backbone
 * 
 * @param sam_file 
 * @param backbone 
 * @return std::string 
 */
std::string alternative_backbone(std::string &sam_file, std::string &backbone){

    ifstream in(sam_file);
    string line;

    vector<bool> replaced_bases(backbone.size(), false);
    vector<string> new_ref(backbone.size(), "");
    for (int c = 0 ; c < backbone.size() ; c++){
        new_ref[c] = backbone[c];
    }


    while (getline(in, line)){

        if (line[0] =='@'){
            continue;
        }

        string field;
        std::istringstream line2(line);
        string cigar;
        string seq = "";
        bool big_indel = false;
        int field_number = 0;
        int start_pos = 0;
        int flag = 0;
        while (getline(line2, field, '\t'))
        {
            if (field_number == 5){
                if (field != "*" && field.size()>2){
                    cigar = convert_cigar(field);
                }
            }
            else if (field_number == 1){
                flag = std::atoi(field.c_str());
            }
            else if (field_number == 9){
                if (flag%16 == 0){
                    seq = field;
                }
                else{
                    seq = reverse_complement(field);
                }
            }
            else if (field_number == 3){
                start_pos = std::atoi(field.c_str());
            }
            field_number++;
        }

        //go through the CIGAR and fill new_ref
        int pos_on_ref = start_pos-1;
        int pos_on_seq = 0;

        for (auto c : cigar){
            if (c == 'M'){
                if (!replaced_bases[pos_on_ref]){
                    new_ref[pos_on_ref] = backbone[pos_on_ref];
                }
                pos_on_ref++;
                pos_on_seq++;
            }
            else if (c == 'I'){
                if (!replaced_bases[pos_on_ref]){
                    new_ref[pos_on_ref-1] += seq[pos_on_seq];
                }
                pos_on_seq++;
            }
            else if (c == 'D'){
                if (!replaced_bases[pos_on_ref]){
                    new_ref[pos_on_ref] = "";
                }
                pos_on_ref++;
            }
            else if (c == 'S'){
                pos_on_seq++;
            }
            else if (c == 'H'){
                pos_on_seq++;
            }

            if (pos_on_ref > 0){
                replaced_bases[pos_on_ref-1] = true;
            }
        }
    }

    //now compute the new reference as a concatenation of new_ref
    string new_reference = "";
    int b = 0;
    for (auto s : new_ref){
        if (replaced_bases[b]){
            new_reference += s;
        }
        b++;
    }

    return new_reference;
}

/**
 * @brief Basic overlap-layout assembler using minimap2 to align the reads
 * 
 * @param read_file 
 * @param MINIMAP 
 * @param tmp_folder 
 * @param id 
 * @return std::string reads stitched together
 */
std::string basic_assembly(std::string read_file, string &MINIMAP, string &tmp_folder, string &id){

    //first do the all-vs-all alignment using minimap2
    string max_allowed_time = "600s";
    string com = "timeout " + max_allowed_time + " " + MINIMAP + " -t 1 -X -cx ava-ont " + read_file + " " + read_file + " > " + tmp_folder + "all_vs_all_"+id+".paf 2>" + tmp_folder + "trash.txt";
    auto res = system(com.c_str());
    if (res != 0){
        // cout << "ERROR minimap2 failed or time outed, while running (m)" << com << endl;
        // exit(1);
        return ""; //this may be a timeout (can happen on centromeric/telomeric regions)
    }
    // cout << "DEBUDJD TOOLS " << endl;
    string paf_file = tmp_folder + "all_vs_all_"+id+".paf";


    //then let's do a basic assembly
    //go trhough the paf file
    ifstream in3(paf_file);
    struct Overlap_minimap{
        string name1;
        int length1;
        int start1;
        int end1;
        string name2;
        int length2;
        int start2;
        int end2;
        bool strand;
    };
    std::unordered_map<string, pair<vector<Overlap_minimap>, vector<Overlap_minimap>>> reads; //associates a read to its neighbor left and right

    string line;
    while (getline(in3, line)){

        std::istringstream iss(line);
        string field;
        short fieldnumber = 0;
        Overlap_minimap overlap1;

        iss >> field;
        overlap1.name1 = field;
        iss >> field;
        overlap1.length1 = std::atoi(field.c_str());
        iss >> field;
        overlap1.start1 = std::atoi(field.c_str());
        iss >> field;
        overlap1.end1 = std::atoi(field.c_str());
        iss >> field;
        if (field == "+"){
            overlap1.strand = true;
        }
        else{
            overlap1.strand = false;
        }
        iss >> field;
        overlap1.name2 = field;
        iss >> field;
        overlap1.length2 = std::atoi(field.c_str());
        iss >> field;
        overlap1.start2 = std::atoi(field.c_str());
        iss >> field;
        overlap1.end2 = std::atoi(field.c_str());
        

        //adjust the overlap if it is just a few bases off the end of the read left
        if (overlap1.start1 <= 10){
            if (overlap1.strand == true && overlap1.start1 < overlap1.start2){
                overlap1.start2 -= overlap1.start1;
                overlap1.start1 = 0;
            }
            else if (overlap1.strand == false && overlap1.start1 < overlap1.length2-overlap1.end2){
                overlap1.end2 += overlap1.start1;
                overlap1.start1 = 0;
            }
        }

        //adjust the overlap if it is just a few bases off the end of the read right
        if (overlap1.end1 >= overlap1.length1-10){
            if (overlap1.strand == true && overlap1.length1 - overlap1.end1 < overlap1.length2-overlap1.end2){
                overlap1.end2 += overlap1.length1 - overlap1.end1;
                overlap1.end1 = overlap1.length1;
            }
            else if (overlap1.strand == false && overlap1.length1 - overlap1.end1 > overlap1.start2){
                overlap1.start2 -= overlap1.length1 - overlap1.end1;
                overlap1.end1 = overlap1.length1;
            }
        }

        //adjust the overlap if it is just a few bases off the end of the read left
        if (overlap1.start2 <= 10){
            if (overlap1.strand == true && overlap1.start2 < overlap1.start1){
                overlap1.start1 -= overlap1.start2;
                overlap1.start2 = 0;
            }
            else if (overlap1.strand == false && overlap1.start2 < overlap1.length1-overlap1.end1){
                overlap1.end1 += overlap1.start2;
                overlap1.start2 = 0;
            }
        }

        //adjust the overlap if it is just a few bases off the end of the read right
        if (overlap1.end2 >= overlap1.length2-10){
            if (overlap1.strand == true && overlap1.length2 - overlap1.end2 < overlap1.length1-overlap1.end1){
                overlap1.end1 += overlap1.length2 - overlap1.end2;
                overlap1.end2 = overlap1.length2;
            }
            else if (overlap1.strand == false && overlap1.length2 - overlap1.end2 < overlap1.start1){
                overlap1.start1 -= overlap1.length2 - overlap1.end2;
                overlap1.end2 = overlap1.length2;
            }
        }

        //now create overlap2
        Overlap_minimap overlap2;
        overlap2.start1 = overlap1.start2;
        overlap2.end1 = overlap1.end2;
        overlap2.start2 = overlap1.start1;
        overlap2.end2 = overlap1.end1;
        overlap2.strand = overlap1.strand;
        overlap2.length1 = overlap1.length2;
        overlap2.length2 = overlap1.length1;
        overlap2.name1 = overlap1.name2;
        overlap2.name2 = overlap1.name1;

        if (reads.find(overlap1.name1) == reads.end()){
            reads[overlap1.name1] = {vector<Overlap_minimap>(), vector<Overlap_minimap>()};
        }
        if (overlap1.start1 == 0 && ((overlap1.strand == true && overlap1.start2 != 0) || (overlap1.strand == false && overlap1.end2 != overlap1.length2))){
            reads[overlap1.name1].first.push_back(overlap1);
        }
        if (overlap1.end1 == overlap1.length1 && ((overlap1.strand == true && overlap1.end2 != overlap1.length2) || (overlap1.strand == false && overlap1.start2 != 0))){
            reads[overlap1.name1].second.push_back(overlap1);
        }

        if (reads.find(overlap2.name1) == reads.end()){
            reads[overlap2.name1] = {vector<Overlap_minimap>(), vector<Overlap_minimap>()};
        }
        if (overlap2.start1 == 0 && ((overlap2.strand == true && overlap2.start2 != 0) || (overlap2.strand == false && overlap2.end2 != overlap2.length2))){
            reads[overlap2.name1].first.push_back(overlap2);
        }
        if (overlap2.end1 == overlap2.length1 && ((overlap2.strand == true && overlap2.end2 != overlap2.length2) || (overlap2.strand == false && overlap2.start2 != 0))){
            reads[overlap2.name1].second.push_back(overlap2);
        }
    }
    in3.close();

    // cout << "ioufpoisdu parsed the paf file " << endl;

    //now that we parsed the paf, let's try to assemble by stitching reads together
    struct contig_parts{
        string read;
        int start;
        int end;
        bool strand;
    };
    vector <contig_parts> new_contig;

    //start with read 1 and extend it right as far as possible
    string current_read;
    for (auto r : reads){
        if (r.second.first.size() > 0 && r.second.second.size() > 0){ //take a contig that has neighbors left and right
            cout << "read right are ";
            for (auto rr : r.second.first){
                cout << rr.name2 << " ";
            }
            current_read = r.first;
            break;
        }
    }
    if (current_read == ""){
        // cout << "DEBUG did not manage to assemblef fdfdqcc " << endl;
        return ""; //did not manage to assemble
    }

    contig_parts first_read;
    first_read.read = current_read;
    first_read.start = 0;
    first_read.end = reads[current_read].first[0].length1;
    first_read.strand = true;

    new_contig.push_back(first_read);
    bool current_strand = true;

    std::unordered_set already_used_contigs = {current_read};

    cout << "starting with read " << current_read << endl;

    string new_current_read = current_read;
    //first extend to the right
    while (new_current_read != ""){
        current_read = new_current_read;
        new_current_read = ""; //means it could not extend
        //look for next read

        
        if (current_strand){

            for (int desperation_level = 0 ; desperation_level < 2 ; desperation_level++){ //first try to take only neighbors that are connected at the other end, else fall back to any neighbor

                for (Overlap_minimap overlap : reads[current_read].second){
                    if (already_used_contigs.find(overlap.name2) == already_used_contigs.end() &&
                        (desperation_level == 1 || reads[overlap.name2].first.size()>0 && reads[overlap.name2].second.size()>0 )){ //then let's extend with the next read
                        contig_parts next_read;
                        next_read.read = overlap.name2;
                        
                        if (overlap.strand){
                            next_read.start = overlap.end2;
                            next_read.end = overlap.length2;
                            next_read.strand = true;
                        }
                        else{
                            next_read.start = 0;
                            next_read.end = overlap.start2;
                            next_read.strand = false;
                        }

                        if (next_read.end - next_read.start > 0){
                            new_contig.push_back(next_read);
                            new_current_read = overlap.name2;
                            current_strand = next_read.strand;
                            already_used_contigs.emplace(current_read);

                        
                            cout << "movingd on tho " << new_current_read << " thanks to overlap : " << endl;
                            // cout << overlap.length1 << " " << overlap.name1 << " " << overlap.start1 << " " << overlap.end1 << " " << overlap.strand << endl;
                            // cout << overlap.length2 << " " << overlap.name2 << " " << overlap.start2 << " " << overlap.end2 << " " << overlap.strand << endl;

                            break;
                        }
                    }
                }
                if (new_current_read != ""){
                    break;
                }
            }
        }
        else{
            for (int desperation_level = 0 ; desperation_level < 2 ; desperation_level++){

                for (Overlap_minimap overlap : reads[current_read].first){
                    if (already_used_contigs.find(overlap.name2) == already_used_contigs.end() &&
                        (desperation_level == 1 || reads[overlap.name2].first.size()>0 && reads[overlap.name2].second.size()>0 )){ //let's extend
                        contig_parts next_read;
                        next_read.read = overlap.name2;
                        
                        if (overlap.strand){
                            next_read.start = 0;
                            next_read.end = overlap.start2;
                            next_read.strand = false;
                        }
                        else{
                            next_read.start = overlap.end2;
                            next_read.end = overlap.length2;
                            next_read.strand = true;
                        }

                        if (next_read.end - next_read.start > 0){
                            new_contig.push_back(next_read);
                            new_current_read = overlap.name2;
                            current_strand = next_read.strand;
                            already_used_contigs.emplace(current_read);
                            cout << "movingd on ztho " << new_current_read << endl;

                            break;
                        }
                    }
                }
                if (new_current_read != ""){
                    break;
                }
            }
        }
        
    }

    //now extend to the left
    // cout << "going legt now" << endl;
    current_read = new_contig[0].read;
    current_strand = new_contig[0].strand;
    vector<contig_parts> new_contig_left_reversed;
    new_current_read = current_read;
    already_used_contigs.clear();
    while (new_current_read != ""){
        current_read = new_current_read;
        new_current_read = ""; //means it could not extend
        //look for next read
        
        if (current_strand){

            for (int desperation_level=0 ; desperation_level < 2 ; desperation_level++){

                for (Overlap_minimap overlap : reads[current_read].first){
                    if (already_used_contigs.find(overlap.name2) == already_used_contigs.end()&&
                        (desperation_level == 1 || reads[overlap.name2].first.size()>0 && reads[overlap.name2].second.size()>0 )){ //then let's extend with the next read
                        contig_parts next_read;
                        next_read.read = overlap.name2;
                        
                        if (overlap.strand){
                            next_read.start = 0;
                            next_read.end = overlap.start2;
                            next_read.strand = true;
                        }
                        else{
                            next_read.start = overlap.end2;
                            next_read.end = overlap.length2;
                            next_read.strand = false;
                        }

                        if (next_read.end - next_read.start > 0){
                            new_contig_left_reversed.push_back(next_read);
                            new_current_read = overlap.name2;
                            current_strand = next_read.strand;
                            already_used_contigs.emplace(current_read);
                            cout << "movingd on Atho " << new_current_read << endl;
                            // cout << overlap.length1 << " " << overlap.name1 << " " << overlap.start1 << " " << overlap.end1 << " " << overlap.strand << endl;
                            // cout << overlap.length2 << " " << overlap.name2 << " " << overlap.start2 << " " << overlap.end2 << " " << overlap.strand << endl;

                            break;
                        }
                    }
                }
                if (new_current_read != ""){
                    break;
                }
            }
        }
        else{

            for (int desperation_level=0 ; desperation_level < 2 ; desperation_level++){

                for (Overlap_minimap overlap : reads[current_read].second){
                    if (already_used_contigs.find(overlap.name2) == already_used_contigs.end()&&
                        (desperation_level == 1 || reads[overlap.name2].first.size()>0 && reads[overlap.name2].second.size()>0 )){ //let's extend
                        contig_parts next_read;
                        next_read.read = overlap.name2;
                        
                        if (!overlap.strand){
                            next_read.start = 0;
                            next_read.end = overlap.start2;
                            next_read.strand = false;
                        }
                        else{
                            next_read.start = overlap.end2;
                            next_read.end = overlap.length2;
                            next_read.strand = true;
                        }

                        if (next_read.end - next_read.start > 0){
                            new_contig_left_reversed.push_back(next_read);
                            new_current_read = overlap.name2;
                            current_strand = next_read.strand;
                            already_used_contigs.emplace(current_read);
                            cout << "movingd on ptho " << new_current_read << endl;

                            break;
                        }
                    }
                }
                if (new_current_read != ""){
                    break;
                }
            }
        }

    }

    vector<contig_parts> new_contig_full;
    //reverse the left part
    for (int i = new_contig_left_reversed.size()-1 ; i >= 0 ; i--){
        new_contig_full.push_back(new_contig_left_reversed[i]);
    }
    //add the rigth part
    for (auto c : new_contig){
        new_contig_full.push_back(c);
    }


    //parse the read files to get the sequences
    std::unordered_map<string, string> reads_sequences;
    ifstream in2(read_file);
    string line2 = "";
    string current_read_name = "";
    string current_read_seq = "";
    while (getline(in2, line2)){
        if (line2[0] == '>'){
            if (current_read_name != ""){
                // cout << "reaqdsljf " << current_read_name << " " << current_read_seq.substr(0,10) << endl;
                reads_sequences[current_read_name] = current_read_seq;
            }
            current_read_name = line2.substr(1);
            current_read_seq = "";
        }
        else{
            current_read_seq += line2;
            // cout << "adding " << line2.substr(0,10) << " to " << current_read_name << endl;
        }
    }
    reads_sequences[current_read_name] = current_read_seq;
    in2.close();   

    // cout << "read1 (from) " << read_file << endl;
    // cout << reads_sequences["read1"] << endl;

    //now convert the contig parts into a string
    string new_contig_string = "";
    for (auto c : new_contig_full){
        cout << "contig part " << c.read << " " << c.start << " " << c.end << " " << c.strand << endl;
        if (c.strand){
            new_contig_string += reads_sequences[c.read].substr(c.start, c.end-c.start);
        }
        else{
            string s = reads_sequences[c.read].substr(c.start, c.end-c.start);
            new_contig_string += reverse_complement(s);
        }
    }

    return new_contig_string;
}

