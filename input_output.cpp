#include <chrono>
#include <list>
#include <set>
#include <fstream>
#include <sstream>
#include <tuple>

#include "input_output.h"
#include "tools.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::tuple;
using std::make_tuple;
using std::get;
using std::make_pair;
using std::string;
using std::ofstream;
using std::ifstream;
using std::array;
using std::set;
using std::stoi;
using std::min;
using std::max;
using robin_hood::unordered_map;


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
    long int linecount = 0;
    long int lastoffset = 0;
    long int offset = 0;

    while(getline(in, line)){

        if (line[0] == format && buffer.size() > 0){
            //then first we append the last read we saw
            // Read r(buffer[1]);
            Read r(""); //append the read without the sequence to be light on memory. The sequences are only needed when they are needed
            r.name = buffer[0];
            r.set_position_in_file(lastoffset+buffer[0].size()+1);
            lastoffset = offset;
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

        linecount++;
        offset += 1+line.size();
    }

    //now append the last read
    Read r("");//append the read without the sequence to be light on memory. The sequences are only needed when they are needed
    r.set_position_in_file(lastoffset+buffer[0].size()+1);
    r.name = buffer[0];
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

/**
 * @brief Uploads the sequence of the reads that align on backbone in allreads
 * 
 * @param fileReads file containing the reads
 * @param backbone index of the backbone read in allreads
 * @param allOverlaps vector of all overlaps between backbone reads and normal reads
 * @param allreads vector of all the reads (including backbone)
 */
void parse_reads_on_contig(std::string fileReads, long int backbone, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads){
 
    ifstream in(fileReads);
    string line;

    //cout << "Loading reads that align on " << allreads[backbone].name << " " << allreads[backbone].neighbors_.size() << endl;

    for (long int n: allreads[backbone].neighbors_){
        long int read;
        if (allOverlaps[n].sequence1 != backbone){
            read = allOverlaps[n].sequence1;
        }
        else{
            read = allOverlaps[n].sequence2;
        }

        in.seekg(allreads[read].get_position_in_file());
        getline(in, line);

        allreads[read].upload_sequence(line);
        // cout << "sequence of read " << allreads[read].name << " is recovered " << line.size() << " " << allreads[read].sequence_.size() << endl;// << line << endl;
    }
    in.close();
}

//input : file containing an assembly in gfa
//output : the contigs appended to the end of allreads and marked as backbones
void parse_assembly(std::string fileAssembly, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices
    , vector<unsigned long int> &backbone_reads, vector<Link> &allLinks){

    ifstream in(fileAssembly);
    if (!in){
        cout << "problem reading files in index_reads, while trying to read " << fileAssembly << endl;
        throw std::invalid_argument( "Input file could not be read" );
    }

    long int sequenceID = allreads.size(); //counting the number of sequences we have already seen 

    string line;
    string nameOfSequence = "";
    string comments = "";


    while(getline(in, line)){

        if (line[0] == 'S'){
            
            string field;
            std::istringstream line2(line);
            int fieldNb = 0;
            unsigned int chunk = 0; //corresponds to lengthOfContig/MAX_SIZE_OF_CONTIGS
            while(getline(line2, field, '\t')){
                if (fieldNb == 1){ // name of the sequence
                    nameOfSequence = field;
                }
                else if (fieldNb == 2){ //here is the sequence

                    Read r(field);
                    r.name = nameOfSequence;
                    backbone_reads.push_back(sequenceID);
                    allreads.push_back(r);

                    //link the minimap name to the index in allreads
                    indices[r.name] = sequenceID;
                    sequenceID++;
                }
                else if (field.substr(0,2) == "dp" || field.substr(0,2) == "DP"){
                    auto depth = std::atoi(field.substr(5, field.size()-5).c_str());
                    //set the depth of all contigs that were created from this sequence
                    allreads[backbone_reads.back()].depth = depth;
                }
                else if (fieldNb > 2){
                    allreads[backbone_reads.back()].comments += "\t"+field;
                }

                fieldNb += 1;
            }
        }


        if (line[0] == 'L'){
            
            string field;
            std::istringstream line2(line);
            int fieldNb = 0;
            Link link;
            
            try{
                string name1 = "";
                string name2 = "";
                while(getline(line2, field, '\t')){
                    if (fieldNb == 1){ // name of the sequence1
                        name1 = field;
                    }
                    else if (fieldNb == 2){ //here is the sequence
                        if (field == "+"){
                            if (indices.find(name1) == indices.end()){
                                cout << "Problem in reading the link : " << line << endl;
                                throw std::invalid_argument( "Problem in file formatting" );
                            }
                            link.neighbor1 = indices[name1];
                            link.end1 = 1;
                        }
                        else if (field == "-"){
                            if (indices.find(name1) == indices.end()){
                                cout << "Problem in reading the link : " << line << endl;
                                throw std::invalid_argument( "Problem in file formatting" );
                            }
                            link.neighbor1 = indices[name1];
                            link.end1 = 0;
                        }
                        else{
                            cout << "Problem in reading the link : " << line << endl;
                            throw std::invalid_argument( "Problem in file formatting" );
                        }
                    }
                    else if (fieldNb == 3){ // name of the sequence1
                        name2 = field;
                    }
                    else if (fieldNb == 4){ //here is the sequence
                        if (field == "+"){
                            if (indices.find(name2) == indices.end()){
                                cout << "Problem in reading the link : " << line << endl;
                                throw std::invalid_argument( "Problem in file formatting" );
                            }
                            link.neighbor2 = indices[name2];
                            link.end2 = 0;
                        }
                        else if (field == "-"){
                            if (indices.find(name2) == indices.end()){
                                cout << "Problem in reading the link : " << line << endl;
                                throw std::invalid_argument( "Problem in file formatting" );
                            }
                            link.neighbor2 = indices[name2];
                            link.end2 = 1;
                        }
                        else{
                            cout << "Problem in reading the link : " << line << endl;
                            throw std::invalid_argument( "Problem in file formatting" );
                        }
                    }
                    else if (fieldNb == 5){
                        link.CIGAR = field;
                    }

                    fieldNb += 1;
                }
                //describe the link between the two contigs
                allLinks.push_back(link);
                allreads[link.neighbor1].add_link(allLinks.size()-1, link.end1);
                allreads[link.neighbor2].add_link(allLinks.size()-1, link.end2);
                // cout << line << endl;
                // cout << "Link beddtween " << name1 << " " << name2 << endl;
                // cout << "Link beddtween " << allreads[link.neighbor1].name << " and " << allreads[link.neighbor2].name << endl;
            }
            catch(...){
                cout << "Problem while reading GFA file " + fileAssembly + ". Ensure that all the contigs described in 'L' lines are present in 'S' lines." << endl;
                throw std::invalid_argument("Invalid GFA");
            }
        }

    }

    for (auto l : allLinks){
        int neighbor1 = l.neighbor1;
        int neighbor2 = l.neighbor2;
        //check if the neighbors are in the backbone
    }

}

/**
 * @brief Parses the PAF files of all the reads aligned on the assembly
 * 
 * @param fileSAM Name of SAM file
 * @param allOverlaps vector containing all the overlaps
 * @param allreads vector containing all the reads as well as the contigs. This is updated with all the overlaps
 * @param indices maps the name of the reads to their index in allreads (comes from parse_reads)
 * @param backbones_reads indicates which of the reads are actually contigs in allreads
 * @param computeBackbones in case where there is no backbone assembly, set this to true. Backbones are then selected from the reads 
 */
void parse_PAF(std::string filePAF, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, vector<unsigned long int> &backbones_reads, bool computeBackbones){

    ifstream in(filePAF);
    if (!in){
        cout << "problem reading PAF file " << filePAF << endl;
        throw std::invalid_argument( "Input file '"+filePAF +"' could not be read" );
    }

    vector<bool> backbonesReads (allreads.size(), true); //for now all reads can be backbones read, they will be filtered afterward

    string line;
    long int linenumber = 0;
    while(getline(in, line)){

        string name1;
        int pos1_1= -1;
        int pos1_2= -1;
        int length1= -1;
        string name2;
        int pos2_1= -1;
        int pos2_2= -1;
        int length2= -1;
        bool positiveStrand;
        float diff = -1; //distance between the sequences [0,1]

        int mapq_quality = 255;

        string cigar;

        bool allgood = true;
        //now go through the fields of the line
        short fieldnumber = 0;
        string field;
        std::istringstream line2(line);
        while (getline(line2, field, '\t'))
        {
            if (fieldnumber == 0){
                try{
                    if (indices.find(field)==indices.end()){
                        cout << "Read in PAF not found in FASTA: " << field << endl; // m54081_181221_163846/4391584/9445_12374 for example
                        allgood = false;
                    }
                    
                    std::istringstream line3(field);
                    string unit;
                    while(getline(line3, unit, ' ')){
                        name1 = unit;
                        break;
                    }

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
                    std::istringstream line3(field);
                    string unit;
                    while(getline(line3, unit, ' ')){
                        name2 = unit;
                        break;
                    }
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
            else if (fieldnumber == 11){
                mapq_quality = stoi(field);
            }
            else if (field.substr(0,5) == "cg:Z:"){
                cigar = field.substr(5, field.size()-5);
                cigar = convert_cigar(cigar);
            }
            else if (field.substr(0,5) == "dv:f:"){
                diff = std::atof(field.substr(5, field.size()-5).c_str());
            }
            //std::cout << "my field is : " << field << std::endl;
            fieldnumber++;
        }

        if (allgood && mapq_quality > 5 && fieldnumber > 11 && name2 != name1){

            //now let's check if this overlap extends to the end of the reads (else, it means only a small portion of the read is aligned)
            bool fullOverlap = false;
            float limit1 = std::min(float(1000), float(pos1_2-pos1_1)/7);
            float limit2 = std::min(float(1000),float(pos2_2-pos2_1)/7);
            if (positiveStrand){

                int slideLeft = min(pos1_1, pos2_1);
                int slideRight = min(length2-pos2_2, length1-pos1_2);

                if (pos1_1 < limit1 || pos2_1 < limit2){ 
                    if (length2-pos2_2 < limit2 || length1-pos1_2 < limit1){
                        fullOverlap = true;
                        pos1_2 += slideRight;
                        pos2_2 += slideRight;
                        pos1_1 -= slideLeft;
                        pos2_1 -= slideLeft;
                    }
                }
            }
            else {

                int slideLeft = min(pos1_1, length2-pos2_2);
                int slideRight = min(pos2_1, length1-pos1_2);

                if (pos1_1 < limit1 || length2-pos2_2 < limit2){
                    if (pos2_1 < limit2 || length1-pos1_2 < limit1){
                        fullOverlap = true;
                        pos1_2 += slideRight;
                        pos2_1 -= slideRight;
                        pos1_1 -= slideLeft;
                        pos2_2 += slideLeft;
                    }
                }
            }
            
            // if (name1 == "m54081_181221_163846/4194380/0_3763"){
            //     cout << "reading overlap with @2_@DRR198813.7882 : " << fullOverlap << " " << length1 << " " << pos1_1 << " " << pos1_2 << " " << name2 << endl;
            //     cout << "on practive, we mark them as overlapping... " << allreads[sequence1].name << " " << allreads[sequence2].name << endl;
            //     cout << sequence1 << " " << sequence2 << endl;
            // }
            //add the overlap if it's a full overlap

            if (fullOverlap){

                // cout << allreads[sequence1].name << endl;
                // if (allreads[sequence1].name == "@2_@DRR198813.1544 1544 length=56567"){
                //     cout << "Mapping 2_@DRR198813.1544 to " << allreads[sequence2].name << endl;
                // }

                //now add the overlap to all the concerned chunks of the contigs
                Overlap overlap;
                overlap.sequence1 = indices[name1];
                overlap.sequence2 = indices[name2];
                overlap.position_2_1 = pos2_1;
                overlap.position_2_2 = pos2_2;
                
                overlap.position_1_1 = pos1_1;
                overlap.position_1_2 = pos1_2;

                overlap.strand = positiveStrand;
                overlap.CIGAR = cigar;
                overlap.diff = diff;

                allreads[overlap.sequence1].add_overlap(allOverlaps.size());
                allreads[overlap.sequence2].add_overlap(allOverlaps.size());
                allOverlaps.push_back(overlap);

                //describe the overlap
                // cout << "Overlap between " << allreads[overlap.sequence1].name << " and " << allreads[overlap.sequence2].name << " with positions " << overlap.position_1_1 << " " <<
                //  overlap.position_1_2 << " " << overlap.position_2_1 << " " << overlap.position_2_2 << endl;

                //now take care of backboneReads: two overlapping reads cannot both be backbone
                if (backbonesReads[overlap.sequence1] && backbonesReads[overlap.sequence2] && pos1_2-pos1_1 > 0.5*min(length1, length2)){
                    if (allreads[overlap.sequence1].sequence_.size() < allreads[overlap.sequence2].sequence_.size()){
                        backbonesReads[overlap.sequence1] = false;
                    }
                    else{
                        backbonesReads[overlap.sequence2] = false;
                    }
                }
            }
            // else{ //DEBUG
            //     cout << "Not full overlap ! " << allreads[sequence1].name << "," << allreads[sequence2].name << " " << positiveStrand
            //     << " " << pos1_1 << "," << pos1_2 << "," << pos2_1 << "," << pos2_2 << " " << length1-pos1_2 << "," << length2-pos2_2 << endl;
            // }
        }
        linenumber++;
    }

    // for (auto b : backbones_reads){
    //     if (allreads[b].name == "contig_102"){
    //         cout << "backbone: " << allreads[b].name << endl;

    //         for (auto c : allreads[b].neighbors_){
    //             cout << "on bb: " << allreads[allOverlaps[c].sequence1].name << endl;
    //         }
    //     }
    // }

    //determine backbone_ reads if asked
    if (computeBackbones){
        for (auto i = 0 ; i<backbonesReads.size() ; i++){
            if (backbonesReads[i]){
                backbones_reads.push_back(i);
            }
        }
    }

}

/**
 * @brief Parses the SAM files of all the reads aligned on the assembly
 * 
 * @param fileSAM Name of SAM file
 * @param allOverlaps vector containing all the overlaps
 * @param allreads vector containing all the reads as well as the contigs
 * @param indices maps the name of the reads to their index in allreads (comes from parse_reads)
 * @param backbonesReads indicates which of the reads are actually contigs in allreads
 */
void parse_SAM(std::string fileSAM, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices){

    ifstream in(fileSAM);
    if (!in){
        cout << "problem reading SAM file " << fileSAM << endl;
        throw std::invalid_argument( "Input file '"+fileSAM +"' could not be read" );
    }

    string line;
    long int linenumber = 0;
    while(getline(in, line)){

        if (line[0] != '@'){
            string field;
            std::istringstream line2(line);

            unsigned long int sequence1 = -1;
            string name1;
            int length1 = 0;
            unsigned long int sequence2= -2;
            string name2;
            int pos2_1= -1;
            bool positiveStrand = true;

            string cigar;

            bool allgood = true;
            //now go through the fields of the line
            short fieldnumber = 0;
            while (getline(line2, field, '\t'))
            {
                if (fieldnumber == 0){
                    try{
                        if (indices.find(field)==indices.end()){
                            cout << "WARNING: read in the sam file not found in fasta file, ignoring: " << field << endl; // m54081_181221_163846/4391584/9445_12374 for example
                            allgood = false;
                        }
                        
                        sequence1 = indices[field];     

                    }
                    catch(...){
                        cout << "There is a sequence in the PAF I did not find in the fasta: " << field << endl;
                        allgood = false;
                    }
                }
                else if (fieldnumber == 1){ //this is the flag
                    int flag = stoi(field);
                    if (flag%8 >= 4 || flag%512 >= 256 || flag>=2048){ //this means that 1) the reads does not map well or 2) this is a secondary alignment 3) is a supplementary alignment
                        allgood = false;
                    }
                    if (flag%32 >= 16){
                        positiveStrand = false;
                    }
                }
                else if (fieldnumber == 2){
                    try{
                        sequence2 = indices[field];
                        std::istringstream line3(field);
                        string unit;
                        while(getline(line3, unit, ' ')){
                            name2 = unit;
                            break;
                        }
                    }
                    catch(...){
                        cout << "There is a sequence in the SAM I did not find in the fasta/q:" << field << ":" << endl;
                        allgood = false;
                    }
                }
                else if (fieldnumber == 3){
                    pos2_1 = stoi(field);
                }
                else if (fieldnumber == 5){
                    cigar = field;
                }
                else if (fieldnumber == 9){
                    length1 = field.size();
                }
                //std::cout << "my field is : " << field << std::endl;
                fieldnumber++;
            }

            if (allgood && fieldnumber > 10 && sequence2 != sequence1){

                //do not store the 'H' and 'S' at the ends of the cigar string
                int nbH_start = 0;
                int nbH_end = 0;
                string string_nbH_start = "";
                //find the first letter in the CIGAR
                int firstH = 0;
                for (int i = 0 ; i < cigar.size() ; i++){
                    if (cigar[i] > '9' || cigar[i] < '0'){
                        if (cigar[i] != 'H') {
                            string_nbH_start = "";
                        }
                        firstH = i;
                        break;
                    }
                    string_nbH_start += cigar[i];
                }
                nbH_start = 0;
                if (string_nbH_start != ""){
                    nbH_start = stoi(string_nbH_start);
                }

                string string_nbH_end = "";
                int lastH = 0;
                //find the last letter in the CIGAR
                for (int i = cigar.size()-1 ; i >= 0 ; i--){
                    if ((cigar[i] > '9' || cigar[i] < '0') && cigar[i] != 'H'){
                        lastH = i;
                        break;
                    }
                    string_nbH_end = cigar[i] + string_nbH_end;
                }
                nbH_end = 0;
                if (string_nbH_end != ""){
                    //strip last character
                    string_nbH_end = string_nbH_end.substr(0, string_nbH_end.size()-1);
                    nbH_end = stoi(string_nbH_end);
                }

                if (!positiveStrand){
                    auto tmp = nbH_start;
                    nbH_start = nbH_end;
                    nbH_end = tmp;
                }

                int nbS_start = 0;
                int nbS_end = 0;
                string string_nbS_start = "";
                //find the first letter in the CIGAR
                int firstS = 0;
                for (int i = 0 ; i < cigar.size() ; i++){
                    if (cigar[i] > '9' || cigar[i] < '0'){
                        if (cigar[i] != 'S') {
                            string_nbS_start = "";
                        }
                        firstS = i;
                        break;
                    }
                    string_nbS_start += cigar[i];
                }
                nbS_start = 0;
                if (string_nbS_start != ""){
                    nbS_start = stoi(string_nbS_start);
                }

                string string_nbS_end = "";
                int lastS = 0;
                //find the last letter in the CIGAR
                for (int i = cigar.size()-1 ; i >= 0 ; i--){
                    if ((cigar[i] > '9' || cigar[i] < '0') && cigar[i] != 'S'){
                        lastS = i;
                        break;
                    }
                    string_nbS_end = cigar[i] + string_nbS_end;
                }
                nbS_end = 0;
                if (string_nbS_end != ""){
                    //strip last character
                    string_nbS_end = string_nbS_end.substr(0, string_nbS_end.size()-1);
                    nbS_end = stoi(string_nbS_end);
                }

                if (!positiveStrand){
                    auto tmp = nbS_start;
                    nbS_start = nbS_end;
                    nbS_end = tmp;
                }

                if (nbH_start+nbH_end > 0.2*length1){
                    allgood = false;
                    //cout << "inpout outpout qkdldkj c " << line << endl;
                }
                else{

                    Overlap overlap;
                    overlap.sequence1 = sequence1;
                    overlap.sequence2 = sequence2;
                    overlap.position_1_1 = nbS_start; //the whole read is used
                    overlap.position_1_2 = nbS_start + length1;
                    overlap.position_2_1 = pos2_1-1; //-1 because the SAM file is 1-based
                    overlap.position_2_2 = pos2_1+length1;
                    overlap.strand = positiveStrand;
                    overlap.CIGAR = cigar;

                    // if (sequence1 == 8641){
                    //     cout << "heereei is the overlap from " << allreads[overlap.sequence1].name.substr(0,10) << " on " << allreads[overlap.sequence2].name.substr(0,10) 
                    //         << " " << overlap.position_1_1 << " " << overlap.position_1_2 << " " << overlap.position_2_1 << " " << overlap.position_2_2 << " " << overlap.strand << endl;
                    //     cout << cigar.substr(0,30) << endl;
                    //     cout << allreads[overlap.sequence1].sequence_.subseq(overlap.position_1_1, overlap.position_1_2-overlap.position_1_1).reverse_complement().str().substr(0,30) << endl;
                    // }

                    // cout << "The overlap I'm adding looks like this: " << overlap.sequence1 << " " << overlap.sequence2 << " " << overlap.strand << endl;

                    allreads[sequence1].add_overlap(allOverlaps.size());
                    if (sequence1 != sequence2){
                        allreads[sequence2].add_overlap(allOverlaps.size());
                    }
                    // if (sequence1 == 0){
                    //     cout << "DEBUG added edge bnbnb " << endl; 
                    // }
                    allOverlaps.push_back(overlap);
                }
            }
            linenumber++;
        }
    }

    in.close();
}

/*
//input : a VCF file containing a list of variants
//output : all variants loaded in a variant list
void parse_VCF(std::string fileVCF, robin_hood::unordered_map<std::string, std::vector <Variant>> &allVariants){
    ifstream in(fileVCF);
    if (!in){
        cout << "problem reading VCF file " << fileVCF << endl;
        throw std::invalid_argument( "Input file '"+fileVCF +"' could not be read" );
    }

    string line;
    while(getline(in, line)){

        if (line[0] != '#'){

            std::istringstream line2(line);
            string field;
            int nbfield = 0;
            Variant v;

            while (getline(line2, field, '\t')){

                if (nbfield == 0){
                    v.refSeq = field;
                }
                else if (nbfield == 1){
                    v.position = std::stoi(field);
                }
                else if (nbfield == 7){
                    std::istringstream field2(field);
                    std::string info;
                    while (getline(field2, info, ';')){
                        if (info.substr(0,6) ==  "RNAMES"){
                            
                            std::istringstream info2(info.substr(7, info.size()));
                            string read;
                            while(getline(info2, read, ',')){
                                cout << "read : " << read << endl;
                                v.readsWithVariant.emplace(read);
                                cout << "size : " << v.readsWithVariant.size() << endl;
                            }

                            break;
                        }
                    }
                }
                nbfield++;
            }
            // cout << "variant in : " << v.readsWithVariant.size() << endl;
            allVariants[v.refSeq].push_back(v);
        }

    }
    for (auto refseqs : allVariants)
    {
        std::sort(refseqs.second.begin() , refseqs.second.end());   
        // cout << "length of the variants : " << endl;
        // for (auto v : refseqs.second){
        //     cout << v.readsWithVariant.size() << endl;
        // } 
    }
}

//input : sam file of all the reads aligned to the reference
//output : updated list of variants, with for each variant also the reads that do not align on the ref but only on the variants
void parseSAM(std::string fileSAM , robin_hood::unordered_map<std::string, std::vector <Variant>> &allVariants){

    ifstream in(fileSAM);
    if (!in){
        cout << "problem reading SAM file " << fileSAM << endl;
        throw std::invalid_argument( "Input file '"+fileSAM +"' could not be read" );
    }

    string line;
    int nbreads = 0;
    while(getline(in, line)){

        if (line[0] != '@'){

            std::istringstream line2(line);
            string field;
            int fieldnb = 0;
            string seqname;
            string refseq;
            int position;
            int length;
            while (getline(line2, field, '\t')){

                if (fieldnb == 0){
                    seqname = field;
                }
                else if (fieldnb == 2){
                    refseq = field;
                }
                else if (fieldnb == 3){
                    // cout << "position : " << field << endl;
                    position = std::stoi(field);
                }
                else if (fieldnb == 9){
                    length = field.size();
                }

                fieldnb++;
            }
            //now check all the variants in this zone
            for (Variant v : allVariants[refseq]){
                if (v.position > position && v.position < position+length){
                    //cout << "the variant " << v.refSeq << " " << v.position << " is found on read " << seqname << endl;
                    if (v.readsWithVariant.find(seqname) == v.readsWithVariant.end()){
                        v.readsWithoutVariant.emplace(seqname);
                    }
                }
            }
            nbreads++;
            if (nbreads%100 == 0){
                // cout << "Read " << nbreads << " reads " << endl;
            }
        }
    }
}
*/

//input : original file of overlaps, allreads and partitions
//output : the same file of overlaps, but with all spurious overlap filtered out
void output_filtered_PAF(std::string fileOut, std::string fileIn, std::vector <Read> &allreads, std::vector<std::vector<int>> &partitions, robin_hood::unordered_map<std::string, unsigned long int> &indices){

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
    int numberOfNotBb = 0;

    while(getline(in, line)){

        std::istringstream line2(line);
        string field;
        short fieldnumber = 0;

        long int sequence1;
        long int sequence2;
        int pos1_1= -1;
        int pos1_2= -1;
        int length1= -1;
        int pos2_1= -1;
        int pos2_2= -1;
        int length2= -1;
        bool positiveStrand;

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
                    name2 = field;
                }
                catch(...){
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
            fieldnumber++;
        }

        //now let's check if this overlap extends to the end of the reads (else, it means the read did not come from the same region)
        bool fullOverlap = false;
        float limit1 = float(pos1_2-pos1_1)/2;
        float limit2 = float(pos2_2-pos2_1)/2;
        if (positiveStrand){
            if (pos1_1  < limit1){
                if (length2-pos2_2 < limit2 || length1-pos1_2 < limit1){
                    fullOverlap = true;
                }
            }
            else if (pos2_1 < limit2){
                    if (length2-pos2_2 < limit2 || length1-pos1_2 < limit1){
                    fullOverlap = true;
                }
            }
        }
        else {
            if (pos1_1 < limit1){
                if (pos2_1 < limit2 || length1-pos1_2 < limit1){
                    fullOverlap = true;
                }
            }
            else if (length2-pos2_2 < limit2){
                if (length1-pos1_2 < limit1 || pos2_1 < limit2){
                    fullOverlap = true;
                }
            }
        }

        //now that we have the two sequences, check if they were partitionned separately
        bool goodOverlap = true;
        bool commonBackbone = false;
        if (allgood){

            for (int backbone1 = 0 ; backbone1 < allreads[sequence1].backbone_seq.size() ; backbone1 ++){
                // if ("1_@SRR8184499.1.4916" == name1){
                //     cout << "backbone of 1_@SRR8184499.1.4916 " << allreads[sequence1].backbone_seq[backbone1].first << endl;
                // }

                for (int backbone2 = 0 ; backbone2 < allreads[sequence2].backbone_seq.size() ; backbone2 ++){

                    int backbone = allreads[sequence1].backbone_seq[backbone1].first;
                    if (backbone == allreads[sequence2].backbone_seq[backbone2].first){ //then they lean on the same backbone read
                        commonBackbone = true;
                        if (name1[0] == name2[0]){
                            if (partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] 
                            != partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second]){

                                if (partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] 
                                * partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second] >= 0){
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
                        }
                        // cout << "comparing " << allreads[sequence2].name << " and " << allreads[sequence1].name << endl;
                        // cout << partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] << " " << partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second] << endl;
                        if (partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] 
                            != partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second]
                                && partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] != -1
                                && partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second] != -1){ //if they're not in the same cluster
                            // if (partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] 
                            // * partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second] >= 0){ //0 is an uncertain read
                                goodOverlap = false;
                                // cout << "not validating " << name1 << " vs " << name2 << endl;
                            // }
                            // else {
                            //     goodOverlap = false;
                            // }
                        }
                        else if  (partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] 
                            == partitions[backbone][allreads[sequence1].backbone_seq[backbone1].second]
                                && partitions[backbone][allreads[sequence2].backbone_seq[backbone2].second] != -1){
                                goodOverlap = true;
                                // cout << "good overlap : " << name1 << " " << name2 << endl;
                        }

                    }
                }
            }

            // if (!commonBackbone){
            //     cout << "no common backbones : " << name1 << " and " << name2 << endl;
            // }
        }
        
        // fullOverlap = true;
        // if (!(goodOverlap && fullOverlap) && name1[0] == name2[0]){
        //     cout << name1 << " and " << name2 << endl;
        //     cout << "reason : " << goodOverlap << " " << fullOverlap << " " << commonBackbone << endl;
        // }
        // else if (goodOverlap && fullOverlap && name1[0] != name2[0]){
        //     cout << name1 << " and " << name2 << endl;
        //     cout << "reason : " << goodOverlap << " " << fullOverlap << " " << commonBackbone << endl;
        // }


        if (!commonBackbone){
            numberOfNotBb++;
        }

        if (goodOverlap && fullOverlap){
            out << line << endl;
        }
    }
    // cout << "number of non-common overlaps : " << numberOfNotBb << endl;
}

/**
 * @brief Creates the GAF corresponding to the mapping of the reads on the new GFA
 * 
 * @param allreads vector of all reads (including backbone reads which can be contigs)
 * @param backbone_reads vector of all the indices of the backbone reads in allreads
 * @param allLinks vector of all links of the GFA
 * @param allOverlaps vector of all overlaps between backbone reads and normal reads
 * @param partitions contains all the conclusions of the separate_reads algorithm
 * @param outputGAF name of the output file
 */
typedef std::tuple<int, vector<pair<string, bool>>, long int> Path; //a path is a starting position on a read, a list of contigs and their orientation relative to the read, and the index of the contig on which it aligns
void output_GAF(std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, std::vector<Link> &allLinks, std::vector <Overlap> &allOverlaps,
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions,
    std::string outputGAF){

    vector<vector<Path>> readPaths (allreads.size()); //to each read we associate a path on the graph

    int max_backbone = backbone_reads.size(); 
    for (int b = 0 ; b < max_backbone ; b++){

        long int backbone = backbone_reads[b];

        if (partitions.find(backbone) != partitions.end() && partitions[backbone].size() > 0){
            for (int n = 0 ; n < allreads[backbone].neighbors_.size() ; n++){

                auto ov = allOverlaps[allreads[backbone].neighbors_[n]];

                long int read;
                long int end;
                int start = -1;
                int nbHS_start = 0;
                int nbHS_end = 0;
                string string_nbHS_start = "";
                //find the first letter in the CIGAR
                for (int i = 0 ; i < ov.CIGAR.size() ; i++){
                    if (ov.CIGAR[i] > '9' || ov.CIGAR[i] < '0'){
                        if (ov.CIGAR[i] != 'H' && ov.CIGAR[i] != 'S') {
                            string_nbHS_start = "";
                        }
                        break;
                    }
                    string_nbHS_start += ov.CIGAR[i];
                }
                nbHS_start = 0;
                if (string_nbHS_start != ""){
                    nbHS_start = stoi(string_nbHS_start);
                }

                string string_nbHS_end = "";
                //find the last letter in the CIGAR
                for (int i = ov.CIGAR.size()-1 ; i >= 0 ; i--){
                    if ((ov.CIGAR[i] > '9' || ov.CIGAR[i] < '0') && ov.CIGAR[i] != 'H' && ov.CIGAR[i] != 'S'){
                        break;
                    }
                    string_nbHS_end = ov.CIGAR[i] + string_nbHS_end;
                }
                nbHS_end = 0;
                if (string_nbHS_end != ""){
                    //strip last character
                    string_nbHS_end = string_nbHS_end.substr(0, string_nbHS_end.size()-1);
                    nbHS_end = stoi(string_nbHS_end);
                }

                if (!ov.strand){
                    auto tmp = nbHS_start;
                    nbHS_start = nbHS_end;
                    nbHS_end = tmp;
                }

                if (ov.sequence1 != backbone){
                    read = ov.sequence1;
                    start = min(ov.position_1_1, ov.position_1_2)+nbHS_start;
                    end = max(ov.position_1_1, ov.position_1_2)-nbHS_end;
                }
                else{
                    read = ov.sequence2;
                    start = min(ov.position_2_1, ov.position_2_2)+nbHS_start;
                    end = max(ov.position_2_1, ov.position_2_2)-nbHS_end;
                }

                //go through all the intervals and see through which version this read passes
                vector<pair<string, bool>> sequence_of_traversed_contigs; //string corresponds to the name of the contig, bool is true if the contig is traversed in the same orientation as the read
                short stop = 0;
                bool firsthere = false;
                bool lasthere = false;
                int inter = 0;
                for (auto interval : partitions[backbone]){

                    if (interval.second.first[n] > -1 && stop < 2){
                        // if (allreads[read].name.substr(0,6) == "@c36dd"){
                        //     cout << "read " << allreads[read].name << " passes input_output yyu through " << interval.first.first << " " << interval.first.second << " " << interval.second.first[n] << endl;
                        //     for (auto i : interval.second.first){
                        //         if (i != -2){
                        //             cout << i << " ";
                        //         }
                        //     }
                        //     cout << endl;
                        // }
                        sequence_of_traversed_contigs.push_back(make_pair(allreads[backbone].name+"_"+std::to_string(interval.first.first)+"_"+std::to_string(interval.second.first[n])
                            , ov.strand));
                        if (inter == 0){
                            firsthere = true;
                        }
                        stop = 1;

                        // if (stop){ //you have one read that was lost just before !
                        //     cout << "WHOUOU: revival here of reads " << ov.sequence1 << endl;
                        // }
                    }
                    else if (stop == 1){
                        stop = 2;
                    }
                    inter++;
                }
                //last contig
                if (stop < 2){
                    lasthere = true;
                    int right = partitions[backbone][partitions[backbone].size()-1].first.second+1;
                    sequence_of_traversed_contigs.push_back(make_pair(allreads[backbone].name+"_"+std::to_string(right)+"_0"
                            , ov.strand));
                }

                if (!ov.strand){ //then mirror the vector
                    std::reverse(sequence_of_traversed_contigs.begin(), sequence_of_traversed_contigs.end());
                }

                if (((ov.strand && !lasthere) || (!ov.strand && !firsthere)) && ((ov.strand && !firsthere) || (!ov.strand && !lasthere))){ //mark if the read does not extend to either end
                    sequence_of_traversed_contigs.push_back(make_pair("&", ov.strand));
                }
                else if ((ov.strand && !lasthere) || (!ov.strand && !firsthere)){ //mark if the read does not extend to the end
                    sequence_of_traversed_contigs.push_back(make_pair("+", ov.strand));
                }
                else if ((ov.strand && !firsthere) || (!ov.strand && !lasthere)){ //mark if the read does not extend to the beginning
                    sequence_of_traversed_contigs.push_back(make_pair("-", ov.strand));
                }
                
                if (sequence_of_traversed_contigs.size()>0){ //this should almost always be true, but it's still safer to test
                    Path path = make_tuple(start,sequence_of_traversed_contigs, backbone);
                    // if (allreads[read].name.substr(0,6) == "@d5282"){
                    //     cout << "d5282 is a readdd, path :: "<< endl;
                    //     for (auto p : get<1>(path)){
                    //         cout << p.first << " " << p.second << endl;
                    //     }

                    // }
                    readPaths[read].push_back(path);
                }
            }
        }
        else{
            for (int n = 0 ; n < allreads[backbone].neighbors_.size() ; n++){
                auto ov = allOverlaps[allreads[backbone].neighbors_[n]];
                long int read;
                int start = -1;
                int nbHS_start = 0;
                int nbHS_end = 0;
                string string_nbHS_start = "";
                //find the first letter in the CIGAR
                for (int i = 0 ; i < ov.CIGAR.size() ; i++){
                    if (ov.CIGAR[i] > '9' || ov.CIGAR[i] < '0'){
                        if (ov.CIGAR[i] != 'H' && ov.CIGAR[i] != 'S') {
                            string_nbHS_start = "";
                        }
                        break;
                    }
                    string_nbHS_start += ov.CIGAR[i];
                }
                nbHS_start = 0;
                if (string_nbHS_start != ""){
                    nbHS_start = stoi(string_nbHS_start);
                }

                string string_nbHS_end = "";
                //find the last letter in the CIGAR
                for (int i = ov.CIGAR.size()-1 ; i >= 0 ; i--){
                    if (ov.CIGAR[i] != 'H' && ov.CIGAR[i] != 'S' && (ov.CIGAR[i] > '9' || ov.CIGAR[i] < '0')){
                        break;
                    }
                    string_nbHS_end = ov.CIGAR[i] + string_nbHS_end;
                }
                nbHS_end = 0;
                if (string_nbHS_end != ""){
                    //strip last character
                    string_nbHS_end = string_nbHS_end.substr(0, string_nbHS_end.size()-1);
                    nbHS_end = stoi(string_nbHS_end);
                }

                if (!ov.strand){
                    auto tmp = nbHS_start;
                    nbHS_start = nbHS_end;
                    nbHS_end = tmp;
                }

                if (ov.sequence1 != backbone){
                    read = ov.sequence1;
                    start = min(ov.position_1_1, ov.position_1_2)+nbHS_start;
                }
                else{
                    read = ov.sequence2;
                    start = min(ov.position_2_1, ov.position_2_2)+nbHS_start;
                }

                
                vector<pair<string, bool>> v = {make_pair(allreads[backbone].name, ov.strand)};
                Path contigpath = make_tuple(start,v, backbone);
                readPaths[read].push_back(contigpath);

                // cout << "bbb " << allreads[backbone_reads[b]].name << endl;
                // if ("edge_2@3@0" == allreads[backbone_reads[b]].name.substr(0,12)){
                //     cout << "on backbone 228 is aligned " << allreads[read].name.substr(0, 23) << endl;
                // }
            }
        }
    }

    //now merge the paths that were on different contigs
    for (auto r = 0 ; r < readPaths.size() ; r++){

        if (readPaths[r].size() > 0){

            std::sort(readPaths[r].begin(), readPaths[r].end(),[] (const auto &x, const auto &y) { return get<0>(x) < get<0>(y); }); //gets the list sorted on first element of pair, i.e. position of contig on read
            vector<Path> mergedPaths;
            Path currentPath = readPaths[r][0];
            //check if each path can be merged with next path
            for (auto p = 0 ; p<readPaths[r].size()-1 ; p++){

                long int contig = get<2> (currentPath);
                bool orientation = get<1>(currentPath)[get<1>(currentPath).size()-1].second;
                long int nextContig = get<2> (readPaths[r][p+1]);

                if (contig != nextContig){
                    vector<size_t> links;
                    if (orientation){
                        links = allreads[contig].get_links_right();
                    }
                    else{
                        links = allreads[contig].get_links_left();
                    }


                    bool merge = false;
                    for (auto li : links){
                        Link l = allLinks[li];
                        if (l.neighbor1 == nextContig || l.neighbor2 == nextContig){ //then merge
                            if ((l.end1==l.end2 && get<1>(readPaths[r][p+1])[0].second != orientation) 
                                || (l.end1!=l.end2 && get<1>(readPaths[r][p+1])[0].second == orientation)){
                                merge = true;
                            }
                        }
                    }

                    //if & in name, there is a cut there
                    char lastchar = get<1>(currentPath)[get<1>(currentPath).size()-1].first[get<1>(currentPath)[get<1>(currentPath).size()-1].first.size()-1];
                    char firstnextchar = get<1>(readPaths[r][p+1])[get<1>(readPaths[r][p+1]).size()-1].first[get<1>(readPaths[r][p+1])[get<1>(readPaths[r][p+1]).size()-1].first.size()-1];
                    if (lastchar == '&' || lastchar == '+' || firstnextchar == '-'){
                        merge = false;
                    }
                    if (lastchar == '&' || lastchar == '+' || lastchar == '-'){
                        get<1>(currentPath).erase(get<1>(currentPath).end()-1);
                    }

                    if (merge){
                        // if (get<1>(currentPath)[0].first.substr(0,12) == "edge_235@0@0" || get<1>(currentPath)[0].first.substr(0,12) == "edge_510@0@0"){
                        //     cout << "current path : " << endl;
                        //     for (auto p : get<1>(currentPath)){
                        //         cout << p.first << " " << p.second << endl;
                        //     }
                        //     cout << "merging 235 with path : in read " << r << endl;
                        //     for (auto p : get<1>(readPaths[r][p+1])){
                        //         cout << p.first << " " << p.second << endl;
                        //     }
                        // }
                        get<1>(currentPath).insert(get<1>(currentPath).end(), get<1> (readPaths[r][p+1]).begin(), get<1> (readPaths[r][p+1]).end());
                        get<2>(currentPath) = get<2> (readPaths[r][p+1]);
                    }
                    else{
                        mergedPaths.push_back(currentPath);
                        currentPath = readPaths[r][p+1];
                    }
                }
            }
            //if & in name, delete it
            char lastchar = get<1>(currentPath)[get<1>(currentPath).size()-1].first[get<1>(currentPath)[get<1>(currentPath).size()-1].first.size()-1];
            if (lastchar == '&' || lastchar == '+' || lastchar == '-'){
                get<1>(currentPath).erase(get<1>(currentPath).end()-1);
            }

            mergedPaths.push_back(currentPath);
            readPaths[r] = mergedPaths;
        }

    }

    //now the paths have been determined, output the file
    ofstream out(outputGAF);
    for (auto p = 0 ; p < readPaths.size() ; p++){
        for (Path path : readPaths[p]){
            if (get<1>(path).size() > 1){
                //output the path
                out << allreads[p].name << "\t-1\t"<< get<0>(path) <<"\t-1\t+\t";
                for (auto contig : get<1>(path)){
                    if (contig.second){
                        out << ">";
                    }
                    else{
                        out << "<";
                    }
                    out << contig.first;
                }
                out << "\t-1\t-1\t-1\t-1\t-1\t255\n";
            }
        }
    }

}

/**
 * @brief Outputs the fasta file
 * 
 * @param allreads vector of all reads (including backbone reads which can be contigs)
 * @param backbone_reads vector of all the indices of the backbone reads in allreads
 * @param fileOut output file
 */
void output_FASTA(std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, std::string fileOut){
    ofstream out(fileOut);
    for (auto r : backbone_reads){
        if (allreads[r].name != "delete_me" && allreads[r].sequence_.size() > 0){
            Read read = allreads[r];
            out << ">" << read.name << " "<< read.comments << "\n" << read.sequence_.str() << "\n" ;    
        }
    }
}

//input : the list of all reads. Among those, backbone reads are actually contigs
//output : a new gfa file with all contigs splitted
void output_GFA(vector <Read> &allreads, vector<unsigned long int> &backbone_reads, string fileOut, vector<Link> &allLinks)
{
    ofstream out(fileOut);
    for (auto r : backbone_reads){
        if (allreads[r].name != "delete_me"){
            Read read = allreads[r];
            out << "S\t"<< read.name << "\t" << read.sequence_.str();// << read.comments; // the comments may not be compatible with the new contigs
            if (read.depth != -1){
                out << "\tDP:f:" << std::to_string(read.depth);
            }
            out << "\n";
        }
    }
    for (auto l : allLinks){
        if (l.end1 != -1 && l.end2 != -1){
            string end1 = "-";
            if (l.end1 == 1) {end1 = "+";}
            string end2 = "-";
            if (l.end2 == 0) {end2 = "+";}
            out << "L\t" << allreads[l.neighbor1].name << "\t" << end1 << "\t" << allreads[l.neighbor2].name 
                << "\t" << end2 << "\t" << l.CIGAR << "\n";
        }
    }
}

/**
 * @brief Outputs the readGroups file
 * 
 * @param readGroupsFile file to write the readGroups in (.txt)
 * @param allreads 
 * @param backbone_reads 
 * @param partitions 
 */
void output_readGroups(std::string readGroupsFile, std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions, std::vector <Overlap> &allOverlaps){

    ofstream out(readGroupsFile);
    for (auto r : backbone_reads){
        Read contig = allreads[r];
        out << "@CONTIG " << contig.name << "\t" << contig.comments << "\n";
        for (auto p : partitions[r]){
            out << "@POS_ON_CONTIG " << p.first.first << " <-> " << p.first.second << "\n";
            for (auto neighbor = 0 ; neighbor < p.second.first.size() ; neighbor++){
                Overlap ov = allOverlaps[contig.neighbors_[neighbor]];
                Read read = allreads[ov.sequence1];
                out << read.name << "\t" << p.second.first[neighbor] << "\t" << ov.sequence1 << "\n";
            }
            out << "\n";
        }
        out << "\n\n";
        
    }

}

void outputTruePar(Partition truePar, std::string id){
    string fileOut = "/home/rfaure/Documents/these/overlap_filtering/species/Escherichia/triploid/trueHaps_"+id+".tsv";
    ofstream out(fileOut);

    int n = 0;
    for (auto read : truePar.getPartition()){
        if (read == 1){
            out << "1";
        }
        else{
            out << "0";
        }
        n++;
    }
    out.close();
}

/*
//input: list of suspicious column
//output: a file written on disk with 0s and 1s
void outputMatrix(std::vector<Column> &snps, std::vector<size_t> suspectPostitions, string id){
    
    string fileOut = "/home/rfaure/Documents/these/overlap_filtering/species/Escherichia/triploid/matrix_"+id+".tsv";
    ofstream out(fileOut);

    int nbReads = 0;
    for (auto c : snps){
        for (auto pp : c.readIdxs){
            if (pp > nbReads){
                nbReads = pp+1;
            }
        }
    }

    vector<vector<int>> matrix (nbReads, vector<int> (suspectPostitions.size(), 0));

    int n = 0;

    for (auto sus : suspectPostitions){
        auto c =  Partition(snps[sus], 0, snps[sus].content[0]);
        int n2 = 0;
        for (auto p : c.getReads()){
            if (c.getPartition()[n2] == 1){
                matrix[p][n] = 1;
            }
            // else{
            //     cout << "content : " << c.getPartition()[n2] << " did not ocm ei" << endl;
            // }
            n2++;
        }
        n++;
    }

    for (auto line : matrix){
        for (auto c : line){
            out << c << "\t";
        }
        out << "\n";
    }
    out.close();
}
*/
void outputGraph(std::vector<std::vector<int>> &adj,std::vector<int> &clusters, std::string fileOut){

    ofstream out(fileOut);

    out << "nodedef>name VARCHAR,label VARCHAR, cluster VARCHAR\n";
    for (auto i = 0 ; i < adj.size() ; i++){
        out << i << ", " << i << ", " <<clusters[i] << "\n";
    }
    out << "edgedef>node1 VARCHAR,node2 VARCHAR, weight DOUBLE\n";
    for (auto i = 0 ; i < adj.size() ; i++){
        // cout << "line " << i << ", length of line i : " << 
        for (auto j=0 ; j < adj[i].size(); j++){
            if (adj[i][j] > 0){
                out << i << ", " << j << ", " << adj[i][j] << "\n";
            }
            // else{
            //     // cout << "Not working on me : " << adj[i][j] << endl;
            // }
        }
        
    }
}

void outputGraph_several_clusterings(std::vector<std::vector<int>> &adj,std::vector<std::vector<int>> &clusters, std::string fileOut){

    ofstream out(fileOut);

    out << "nodedef>name VARCHAR,label VARCHAR, ";
    for (auto i = 0 ; i < clusters.size() ; i++){
        out << "cluster_" << i << " VARCHAR, ";
    }
    out << "\n";
    for (auto i = 0 ; i < adj.size() ; i++){
        out << i << ", " << i << ", ";
        for (auto j = 0 ; j < clusters.size() ; j++){
            out << clusters[j][i] << ", ";
        }
        out << "\n";
    }
    out << "edgedef>node1 VARCHAR,node2 VARCHAR, weight DOUBLE\n";
    for (auto i = 0 ; i < adj.size() ; i++){
        // cout << "line " << i << ", length of line i : " << 
        for (auto j=0 ; j < adj[i].size(); j++){
            if (adj[i][j] > 0){
                out << i << ", " << j << ", " << adj[i][j] << "\n";
            }
            // else{
            //     // cout << "Not working on me : " << adj[i][j] << endl;
            // }
        }
        
    }
}

/**
 * @brief Output graph in simple format (src dest)
 * 
 * @param adj 
 * @param fileOut 
 */
void output_simple_graph(std::vector<std::vector<int>> &adj, std::string fileOut){
    
        ofstream out(fileOut);
    
        for (auto i = 0 ; i < adj.size() ; i++){
            for (auto j=0 ; j < adj[i].size(); j++){
                if (adj[i][j] > 0){
                    out << i << " " << j << "\n";
                }
            }
        }
        out.close();
}



