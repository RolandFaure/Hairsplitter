#include <chrono>
#include <list>
#include <set>
#include <fstream>
#include <sstream>

#include "input_output.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::string;
using std::ofstream;
using std::ifstream;
using std::array;
using std::set;
using std::stoi;
using robin_hood::unordered_map;
using namespace std::chrono;

//input : file containing all reads in fastq or fasta format
//output : all reads stored in allreads
void parse_reads(std::string fileReads, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices){

    char format = '@'; //a character to keep track of whether the input file is a fasta or a fastq
    if ((fileReads.size()>6 && fileReads.substr(fileReads.size()-6,6) == ".fasta") || fileReads.substr(fileReads.size()-3,3) == ".fa"){
        format = '>';
    }

    ifstream in(fileReads);
    if (!in){
        cout << "problem reading files in index_reads, while trying to read " << fileReads << endl;
        throw std::invalid_argument( "Input file could not be read" );
    }

    long int sequenceID = 0; //counting the number of sequences we have already seen 

    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format && buffer.size() > 0){
            //then first we append the last read we saw
            Read r(buffer[1]);
            allreads.push_back(r);

            ///compute the name of the sequence as it will appear in minimap (i.e. up to the first blank space)
            string nameOfSequence = "";
            for (unsigned int i = 1 ; i < buffer[0].size() ; i++){
                if (buffer[0][i] == ' '){
                    break;
                }
                else{
                    nameOfSequence.push_back(buffer[0][i]);
                }
            }

            ///link the minimap name to the index in allreads
            indices[nameOfSequence] = sequenceID;
            sequenceID++;

            //then we reset the buffer
            buffer = {line};
        }
        else {
            buffer.push_back(line);
        }

    }

    //now append the last read
    Read r(buffer[1]);
    allreads.push_back(r);

    ///compute the name of the sequence as it will appear in minimap (i.e. up to the first blank space)
    string nameOfSequence = "";
    for (unsigned int i = 1 ; i < buffer[0].size() ; i++){
        if (buffer[0][i] == ' '){
            break;
        }
        else{
            nameOfSequence.push_back(buffer[0][i]);
        }
    }
    ///link the minimap name to the index in allreads
    indices[nameOfSequence] = sequenceID;
    sequenceID++;

}

//input : file containing an assembly in fasta format
//output : the contigs appended to the end of allreads and marked as backbones
void parse_assembly(std::string fileAssembly, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, vector<unsigned long int> &backbone_reads){

    char format = '>';

    ifstream in(fileAssembly);
    if (!in){
        cout << "problem reading files in index_reads, while trying to read " << fileAssembly << endl;
        throw std::invalid_argument( "Input file could not be read" );
    }

    long int sequenceID = allreads.size(); //counting the number of sequences we have already seen 

    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format && buffer.size() > 0){
            //then first we append the last read we saw
            Read r(buffer[1]);
            backbone_reads.push_back(sequenceID);
            allreads.push_back(r);

            ///compute the name of the sequence as it will appear in minimap (i.e. up to the first blank space)
            string nameOfSequence = "";
            for (unsigned int i = 1 ; i < buffer[0].size() ; i++){
                if (buffer[0][i] == ' '){
                    break;
                }
                else{
                    nameOfSequence.push_back(buffer[0][i]);
                }
            }

            ///link the minimap name to the index in allreads
            indices[nameOfSequence] = sequenceID;
            sequenceID++;

            //then we reset the buffer
            buffer = {line};
        }
        else {
            buffer.push_back(line);
        }

    }

    //now append the last read
    Read r(buffer[1]);
    backbone_reads.push_back(allreads.size());
    allreads.push_back(r);


    ///compute the name of the sequence as it will appear in minimap (i.e. up to the first blank space)
    string nameOfSequence = "";
    for (unsigned int i = 1 ; i < buffer[0].size() ; i++){
        if (buffer[0][i] == ' '){
            break;
        }
        else{
            nameOfSequence.push_back(buffer[0][i]);
        }
    }
    ///link the minimap name to the index in allreads
    indices[nameOfSequence] = sequenceID;
    sequenceID++;

}

//input : a file containing overlaps
//output : the set of all overlaps, updated allreads with overlaps, and optionnaly a list of backbone reads
void parse_PAF(std::string filePAF, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, vector<unsigned long int> &backbones_reads, bool computeBackbones){

    ifstream in(filePAF);
    if (!in){
        cout << "problem reading PAF file " << filePAF << endl;
        throw std::invalid_argument( "Input file could not be read" );
    }


    vector<bool> backbonesReads (allreads.size(), true); //for now all reads can be backbones read, they will be filtered afterward

    string line;

    while(getline(in, line)){

        string field;
        std::istringstream line2(line);

        unsigned long int sequence1 = 10000;
        int pos1_1= -1;
        int pos1_2= -1;
        int length1= -1;
        unsigned long int sequence2= 100000;
        int pos2_1= -1;
        int pos2_2= -1;
        int length2= -1;
        bool positiveStrand;

        bool allgood = true;
        //now go through the fields of the line
        short fieldnumber = 0;
        while (std::getline(line2, field, '\t'))
        {
            if (fieldnumber == 0){
                try{
                    sequence1 = indices[field];
                }
                catch(...){
                    cout << "There is a sequence in the PAF I did not find in the fasta: " << field << endl;
                    allgood = false;
                }
            }
            else if (fieldnumber == 1){
                length1 = stoi(field);
            }
            else if (fieldnumber == 2){
                pos1_1 =  stoi(field);
            }
            else if (fieldnumber == 3){
                pos1_2 = stoi(field);
            }
            else if (fieldnumber == 4){
                positiveStrand = (field == "+");
            }
            else if (fieldnumber == 5){
                try{
                    sequence2 = indices[field];
                }
                catch(...){
                    cout << "There is a sequence in the PAF I did not find in the fasta/q:" << field << ":" << endl;
                    allgood = false;
                }
            }
            else if (fieldnumber == 6){
                length2 = stoi(field);
            }
            else if (fieldnumber == 7){
                pos2_1 =  stoi(field);
            }
            else if (fieldnumber == 8){
                pos2_2 = stoi(field);
            }
            //std::cout << "my field is : " << field << std::endl;
            fieldnumber++;
        }

        if (allgood && fieldnumber > 11 && sequence2 != sequence1){

            Overlap overlap;
            overlap.sequence1 = sequence1;
            overlap.sequence2 = sequence2;
            overlap.position_1_1 = pos1_1;
            overlap.position_1_2 = pos1_2;
            overlap.position_2_1 = pos2_1;
            overlap.position_2_2 = pos2_2;
            overlap.strand = positiveStrand;

            //cout << "The overlap I'm adding looks like this: " << overlap.sequence1 << " " << overlap.sequence2 << " " << overlap.strand << endl;

            allreads[sequence1].add_overlap(allOverlaps.size());
            if (sequence1 != sequence2){
                allreads[sequence2].add_overlap(allOverlaps.size());
            }
            allOverlaps.push_back(overlap);

            //now take care of backboneReads: two overlapping reads cannot both be backbone
            if (backbonesReads[sequence1] && backbonesReads[sequence2]){
                if (allreads[sequence1].sequence_.size() < allreads[sequence2].sequence_.size()){
                    backbonesReads[sequence1] = false;
                }
                else{
                    backbonesReads[sequence2] = false;
                }
            }
        }
    }


    //determine backbone_ reads if asked
    if (computeBackbones){
        for (auto i = 0 ; i<backbonesReads.size() ; i++){
            if (backbonesReads[i]){
                backbones_reads.push_back(i);
            }
        }
    }
    //backbones_reads = {0}; //just for fun now


}

//input : original file of overlaps, allreads and partitions
//output : the same file of overlaps, but with all spurious overlap filtered out
void output_filtered_PAF(std::string fileOut, std::string fileIn, std::vector <Read> &allreads, std::vector<Partition> &partitions, robin_hood::unordered_map<std::string, unsigned long int> &indices){

    ifstream in(fileIn);
    if (!in){
        cout << "problem reading PAF file " << fileIn << endl;
        throw std::invalid_argument( "Input file could not be read" );
    }

    ofstream out(fileOut);
    if (!out){
        cout << "problem opening " << fileOut << ". Do you have the right permissions ?" << endl;
        throw std::invalid_argument( "Output file could not be written" );
    }

    string line;

    while(getline(in, line)){

        std::istringstream line2(line);
        string field;
        short fieldnumber = 0;

        long int sequence1;
        long int sequence2;

        string name1;
        string name2;

        bool allgood = true;

        while (std::getline(line2, field, '\t'))
        {
            if (fieldnumber == 0){
                try{
                    sequence1 = indices[field];
                    name1 = field;
                }
                catch(...){
                    allgood = false;
                }
            }
            else if (fieldnumber == 5){
                try{
                    sequence2 = indices[field];
                    name2 = field;
                }
                catch(...){
                    allgood = false;
                }
            }
            fieldnumber++;
        }

        //now that we have the two sequences, check if they were partitionned separately
        bool goodOverlap = true;
        if (allgood){

            for (int backbone1 = 0 ; backbone1 < allreads[sequence1].backbone_seq.size() ; backbone1 ++){
                for (int backbone2 = 0 ; backbone2 < allreads[sequence2].backbone_seq.size() ; backbone2 ++){
                    int backbone = allreads[sequence1].backbone_seq[backbone1].first;
                    if (backbone == allreads[sequence2].backbone_seq[backbone2].first){ //then they lean on the same backbone read
                        if (name1[0] == name2[0]){
                            if (partitions[backbone].getPartition()[allreads[sequence2].backbone_seq[backbone2].second] 
                            != partitions[backbone].getPartition()[allreads[sequence1].backbone_seq[backbone1].second]){
                                // cout << "comparing badly " << name1 << " " << name2 << endl;
                                // partitions[backbone].print();
                                // partitions[backbone].getConfidence();
                                // for (auto i : partitions[backbone].getConfidence()){
                                //     cout << i << ",";
                                // }
                                // cout << endl;
                                // cout << partitions[backbone].getConfidence()[allreads[sequence2].backbone_seq[backbone2].second]<< ","<<
                                // partitions[backbone].getConfidence()[allreads[sequence1].backbone_seq[backbone1].second] << endl;
                            }
                        }
                        if (partitions[backbone].getPartition()[allreads[sequence2].backbone_seq[backbone2].second] 
                            != partitions[backbone].getPartition()[allreads[sequence1].backbone_seq[backbone1].second]){
                                goodOverlap = false;
                                //cout << "not validating " << name1 << " vs " << name2 << endl;
                        }

                    }
                }
            }
        }
        if (goodOverlap){
            out << line << endl;
        }

    }

}








