#include "clean_graph.h"

#include <fstream>
#include <unordered_map>
#include <set>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stoi;
using std::unordered_map;
using std::pair;
using std::set;

extern string MINIMAP;

/**
 * @brief Supress small (parts of) contigs that are not linked properly in the GFA and that are found almost identically elsewhere
 * 
 * @param assemblyFile input assembly file (fasta or gfa)
 * @param outputFile cleaned assembly file to write (fasta or gfa)
 * @param num_threads number of threads to use
 * @param outFolder outputFolder (without the /tmp)
 */
void clean_graph(
    std::string& assemblyFile,
    std::string& outputFile,
    std::string& logFile,
    int num_threads,
    std::string &outFolder,
    int &nb_of_deleted_contigs,
    int &length_of_deleted_contigs){

    nb_of_deleted_contigs = 0;
    length_of_deleted_contigs = 0;
    string gfaFile;
    string fastaFile;
    string format;
    unordered_map<string, int> contigs;
    std::vector<int> contigLengths;
    std::vector<pair<set<string>, set<string>>> contigNeighbors;
    if (assemblyFile.substr(assemblyFile.size()-4, 4) == ".gfa"){
        format = "gfa";
        gfaFile = assemblyFile;
        //convert 
        fastaFile =  outFolder+"/tmp/assembly.fa";
        convert_GFA_to_FASTA(assemblyFile, fastaFile);

        //go through the GFA and inventoriate all the contigs, their length and the number of links they have
        std::ifstream gfa(gfaFile);
        std::string line;
        while (std::getline(gfa, line)){
            if (line[0] == 'S'){
                string field;
                std::istringstream line2(line);
                int fieldNumber = 0;
                string nameOfContig;
                int lengthOfContig;
                while (getline(line2, field, '\t')){
                    if (fieldNumber == 1){
                        nameOfContig = field;
                    }
                    else if (fieldNumber == 2){
                        lengthOfContig = field.size();
                    }
                    fieldNumber++;
                }
                contigs[nameOfContig] = contigLengths.size();
                contigLengths.push_back(lengthOfContig);
                contigNeighbors.push_back(pair<set<string>, set<string>>(set<string>(), set<string>()));
            }
            else if (line[0] == 'L'){
                string field;
                std::istringstream line2(line);
                int fieldNumber = 0;
                string nameOfContig1;
                string nameOfContig2;
                bool orientation1;
                bool orientation2;
                while (getline(line2, field, '\t')){
                    if (fieldNumber == 1){
                        nameOfContig1 = field;
                    }
                    else if (fieldNumber == 2){
                        if (field == "+"){
                            orientation1 = true;
                        }
                        else{
                            orientation1 = false;
                        }
                    }
                    else if (fieldNumber == 3){
                        nameOfContig2 = field;
                    }
                    else if (fieldNumber == 4){
                        if (field == "+"){
                            orientation2 = true;
                        }
                        else{
                            orientation2 = false;
                        }
                    }
                    fieldNumber++;
                }
                if (orientation1){
                    contigNeighbors[contigs[nameOfContig1]].first.insert(nameOfContig2);
                }
                else{
                    contigNeighbors[contigs[nameOfContig1]].second.insert(nameOfContig2);
                }
                if (orientation2){
                    contigNeighbors[contigs[nameOfContig2]].first.insert(nameOfContig1);
                }
                else{
                    contigNeighbors[contigs[nameOfContig2]].second.insert(nameOfContig1);
                }
            }
        }
    }
    else if (assemblyFile.substr(assemblyFile.size()-3, 3) == ".fa" || assemblyFile.substr(assemblyFile.size()-6, 6) == ".fasta"){
        //convert the fasta file to gfa to load it
        fastaFile = assemblyFile;
        gfaFile = outFolder+"/tmp/assembly.gfa";
        convert_FASTA_to_GFA(assemblyFile, gfaFile);
        format = "fasta";

        //go through the fasta and inventoriate all the contigs and their length 
        std::ifstream fasta(fastaFile);
        std::string line;

        string nameOfContig;

        while (std::getline(fasta, line)){
            if (line[0] == '>'){
                nameOfContig = line.substr(1, line.size()-1);
                contigs[nameOfContig] = contigLengths.size();
                contigLengths.push_back(0);
            }
            else{
                contigLengths[contigs[nameOfContig]] += line.size();
            }
        }
    }

    //align the assembly against itself using minimap2
    string samFile = outFolder+"/tmp/assembly_against_itself.sam";
    
    string command = MINIMAP+" -t "+std::to_string(num_threads)+" -ax asm5 --no-long-join -N 2 "+fastaFile+" "+fastaFile+" > "+samFile + " 2> "+outFolder+"/tmp/minimap2.log";
    auto minimaprun = system(command.c_str());
    if (minimaprun != 0){
        std::cerr << "Error running minimap2 in clean_graph.cpp: IOUX.\n" << command << std::endl;
        exit(1);
    }

    //parse the sam file and list the contigs that are mostly contained in another contig
    std::vector<std::string> contained_contigs;
    std::ofstream log(logFile, std::ios_base::app);
    log << "STAGE 1 : Deleting contigs that look like parts of haplotypes but that are disconnected from the rest of the assembly\n";
    std::ifstream sam(samFile);
    std::string line;

    while (std::getline(sam, line)){
        if (line[0] == '@'){
            continue;
        }
        string field;
        std::istringstream line2(line);
        string nameOfContig;
        string contigItMapsOn;
        string cigar;
        int fieldNumber = 0;
        int NM = 0;
        bool sameOrientation = true;
        while (getline(line2, field, '\t')){
            
            if (fieldNumber == 0){
                nameOfContig = field;
            }
            //check the flag to see if the contig is mapped in the same orientation
            else if (fieldNumber == 1){
                if (stoi(field) & 16){
                    sameOrientation = false;
                }
            }
            else if (fieldNumber == 2){
                contigItMapsOn = field;
            }
            else if (fieldNumber == 5){
                cigar = field;
            }
            //look for the NM tag
            else if (field.substr(0, 5) == "NM:i:"){
                NM = stoi(field.substr(5, field.size()-5));
            }
            fieldNumber++;

        }
        if (contigItMapsOn == nameOfContig){
            continue;
        }
        //count the number of S/H at the beginning and at the end of the CIGAR and the total length of the query, for example 30S50M45S -> 20,45, lengthOfContig = 125
        int cigarStart = 0;
        int cigarEnd = 0;
        int cigarNumber = 0;
        bool firstLetter = true;
        int lastLetter = 0;
        int lengthOfContig = 0;
        for (char c : cigar){
            if (c == 'S' || c == 'H'){
                if (cigarStart == 0){
                    cigarStart = std::stoi(cigar.substr(lastLetter, cigarNumber-lastLetter));
                    lengthOfContig += cigarStart;
                }
                else{
                    cigarEnd = std::stoi(cigar.substr(lastLetter, cigarNumber-lastLetter));
                    lengthOfContig += cigarEnd;
                }
                lastLetter = cigarNumber+1;
            }
            else if (c == 'M' || c == 'I' || c == 'D'){
                if (c == 'M' || c == 'I'){
                    lengthOfContig += std::stoi(cigar.substr(lastLetter, cigarNumber-lastLetter));
                }
                lastLetter = cigarNumber+1;
            }
            cigarNumber++;
        }
        // cout << cigar << " " << lengthOfContig << endl;
        //if the contig is mostly contained in another contig, add it to the list
        if (cigarStart < 0.1*lengthOfContig && cigarEnd < 0.1*lengthOfContig){
            //check if the contig is disconnected from the rest of the graph
            bool disconnected = false;
            if (contigNeighbors[contigs[nameOfContig]].first.size() == 0 && contigNeighbors[contigs[nameOfContig]].second.size() == 0){
                disconnected = true;
                contained_contigs.push_back(nameOfContig);
                nb_of_deleted_contigs++;
                length_of_deleted_contigs += lengthOfContig;
            }
            else if (contigNeighbors[contigs[nameOfContig]].first.size() == 0){
                //check if the contig has the same neihbor as the contig it maps on
                if (sameOrientation && contigNeighbors[contigs[nameOfContig]].second == contigNeighbors[contigs[contigItMapsOn]].second){
                    disconnected = true;
                    contained_contigs.push_back(nameOfContig);
                    nb_of_deleted_contigs++;
                    length_of_deleted_contigs += lengthOfContig;
                }
                else if (!sameOrientation && contigNeighbors[contigs[nameOfContig]].second == contigNeighbors[contigs[contigItMapsOn]].first){
                    disconnected = true;
                    contained_contigs.push_back(nameOfContig);
                    nb_of_deleted_contigs++;
                    length_of_deleted_contigs += lengthOfContig;
                }
            }
            else if (contigNeighbors[contigs[nameOfContig]].second.size() == 0){
                //check if the contig has the same neihbor as the contig it maps on
                if (sameOrientation && contigNeighbors[contigs[nameOfContig]].first == contigNeighbors[contigs[contigItMapsOn]].first){
                    disconnected = true;
                    contained_contigs.push_back(nameOfContig);
                    nb_of_deleted_contigs++;
                    length_of_deleted_contigs += lengthOfContig;
                }
                else if (!sameOrientation && contigNeighbors[contigs[nameOfContig]].first == contigNeighbors[contigs[contigItMapsOn]].second){
                    disconnected = true;
                    contained_contigs.push_back(nameOfContig);
                    nb_of_deleted_contigs++;
                    length_of_deleted_contigs += lengthOfContig;
                }
            }
            if (disconnected){
                log << "Contig " << nameOfContig << " is disconnected from the rest of the graph and contained in " << contigItMapsOn << endl;
            }
        }
        
    }
    sam.close();
    log.close();

    // cout << "Here are all the contained contigs: " << endl;
    // for (string contig : contained_contigs){
    //     cout << contig << endl;
    // }

    //rewrite the assembly graph without the contained contigs in outputFile
    std::ofstream output(outputFile);
    if (format == "gfa"){
        std::ifstream gfa(assemblyFile);
        while (std::getline(gfa, line)){
            if (line[0] == 'S'){
                string field;
                std::istringstream line2(line);
                string nameOfContig1;
                string nameOfContig2;
                int fieldNumber = 0;
                string nameOfContig;
                while (getline(line2, field, '\t')){
                    if (fieldNumber == 1){
                        nameOfContig = field;
                    }
                    fieldNumber++;
                }
                if (std::find(contained_contigs.begin(), contained_contigs.end(), nameOfContig) == contained_contigs.end()){
                    output << line << std::endl;
                }
            }
            else if (line[0] == 'L'){
                
                string field;
                std::istringstream line2(line);
                string nameOfContig1;
                string nameOfContig2;
                int fieldNumber = 0;
                while (getline(line2, field, '\t')){
                    if (fieldNumber == 1){
                        nameOfContig1 = field;
                    }
                    else if (fieldNumber == 3){
                        nameOfContig2 = field;
                    }
                    fieldNumber++;
                }
                
                if (std::find(contained_contigs.begin(), contained_contigs.end(), nameOfContig1) == contained_contigs.end() && std::find(contained_contigs.begin(), contained_contigs.end(), nameOfContig2) == contained_contigs.end()){
                    output << line << std::endl;
                }
            }
            else{
                output << line << std::endl;
            }
        }
    }
    else if (format=="fasta"){
        std::ifstream fasta(assemblyFile);
        bool write = true;
        while (std::getline(fasta, line)){
            if (line[0] == '>'){
                string nameOfContig = line.substr(1, line.find(' ')-1);
                if (std::find(contained_contigs.begin(), contained_contigs.end(), nameOfContig) == contained_contigs.end()){
                    output << line << std::endl;
                }
                else{
                    write = false;
                }
            }
            else if (write){
                output << line << std::endl;
            }
        }
    }
}