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


/**
 * @brief Parses a fasta or fastq file and stores the reads in allreads
 * 
 * @param fileReads file containing all reads in fastq or fasta format
 * @param allreads vector to store the reads
 * @param indices maps the name of a read to its index in allreads
 */
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

        if (line[0] == format && (buffer.size() == 4 && format == '@' || buffer.size() == 2 && format == '>')){
            //then first we append the last read we saw

            ///parse the name of the sequence as it will appear in minimap (i.e. up to the first blank space)
            string nameOfSequence = "";
            std::stringstream line2(buffer[0].substr(1)); //exclude the first character (either @ or >)
            std::getline(line2, nameOfSequence, ' ');

            Read r("", buffer[1].size()); //append the read without the sequence to be light on memory. The sequences are only needed when they are needed
            r.name = nameOfSequence;
            r.set_position_in_file(lastoffset+buffer[0].size()+1);
            lastoffset = offset;
            allreads.push_back(r);

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
    string nameOfSequence = "";
    std::stringstream line2(buffer[0].substr(1)); //exclude the first character (either @ or >)
    std::getline(line2, nameOfSequence, ' ');

    Read r("", buffer[1].size());//append the read without the sequence to be light on memory. The sequences are only needed when they are needed
    r.set_position_in_file(lastoffset+buffer[0].size()+1);
    r.name = nameOfSequence;
    allreads.push_back(r);

    ///link the minimap name to the index in allreads
    indices[nameOfSequence] = sequenceID;
    sequenceID++;

}

/**
 * @brief Parses a gfa file and stores the contigs in allreads
 * 
 * @param fileAssembly GFA file containing the assembly
 * @param allreads Vector allreads appended with the contigs
 * @param indices Indices linking the name of the contigs to their index in allreads
 * @param backbone_reads Index of the sequences in allreads that are contigs
 * @param allLinks List of all the links between contigs
 */
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
                    //name of sequence as is in minimap file is field up to the first blank space
                    std::stringstream line3(field);
                    std::getline(line3, nameOfSequence, ' ');                    
                }
                else if (fieldNb == 2){ //here is the sequence

                    Read r(field, field.size());
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
            int flag;

            string cigar;

            bool allgood = true;
            //now go through the fields of the line
            short fieldnumber = 0;

            //print all keys of indices
            // for (auto i : indices){
            //     cout << "flkdjs : fqdsjf input output " <<  i.first << "\n";
            // }

            while (getline(line2, field, '\t'))
            {
                if (fieldnumber == 0){
                    try{
                        if (indices.find(field)==indices.end()){
                            cout << "WARNING: read in the sam file not found in reads file, ignoring: " << field << endl; // m54081_181221_163846/4391584/9445_12374 for example
                            allgood = false;
                        }
                        
                        sequence1 = indices[field];     

                    }
                    catch(...){
                        cout << "There is a sequence in the SAM I did not find in the fasta: " << field << endl;
                        allgood = false;
                    }
                }
                else if (fieldnumber == 1){ //this is the flag
                    flag = stoi(field);
                    if (flag%8 >= 4 || flag%512 >= 256){ //this means that 1) the reads does not map well or 2) this is a secondary alignment 3) flag>2048 is a supplementary alignment but can just be the follow-up of the first alignment
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

                int length_alignment = 0;
                string alignment = convert_cigar(cigar);
                for (auto a : alignment){
                    if (a == 'M' || a == 'D' || a == '=' || a == 'X'){
                        length_alignment++;
                    }
                }
                auto pos2_2 = pos2_1 + length_alignment;

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

                if (nbH_start+nbH_end > 0.2*length1 && flag < 2048){ //flag>=2048 means it is a supplementary alignment, in which case H can generally be changed in S. S alignments can be tolerated, if the read does not align elsewhere )
                    allgood = false;
                    //cout << "inpout outpout qkdldkj c " << line << endl;
                }
                else{

                    Overlap overlap;
                    overlap.sequence1 = sequence1;
                    overlap.sequence2 = sequence2;
                    overlap.position_1_1 = nbS_start + nbH_start; //the whole read is used
                    overlap.position_1_2 = nbS_start + nbH_start + length1;
                    overlap.position_2_1 = pos2_1-1; //-1 because the SAM file is 1-based
                    overlap.position_2_2 = pos2_1+length1;
                    overlap.strand = positiveStrand;
                    overlap.CIGAR = cigar;

                    // if (allreads[sequence1].name.substr(0,19) == ">SRR14289618.827569"){
                    //     cout << "heereei is the overlap from " << allreads[overlap.sequence1].name << " on " << allreads[overlap.sequence2].name
                    //         << " " << overlap.position_1_1 << " " << overlap.position_1_2 << " " << overlap.position_2_1 << " " << overlap.position_2_2 << " " << overlap.strand << endl;
                    //     cout << cigar.substr(0,30) << endl;
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
void parse_PAF(std::string filePAF, 
    std::vector <Overlap> &allOverlaps, 
    std::vector <Read> &allreads, 
    robin_hood::unordered_map<std::string, unsigned long int> &indices,
    vector<unsigned long int> &backbones_reads, 
    bool computeBackbones,
    bool filterUncompleteAlignments){

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

        // cout << "fqdio indice " << endl;
        // for (auto i : indices){
        //     cout << i.first << " " << i.second << endl;
        // }
        // exit(0);

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
                        exit(0);
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

        if (((allgood && mapq_quality > 5) || !filterUncompleteAlignments) && fieldnumber > 11 && name2 != name1){

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

            if (fullOverlap || !filterUncompleteAlignments){

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
            out << "S\t"<< read.name << "\t" << read.sequence_.str() << "\t";// << read.comments; // the comments may not be compatible with the new contigs
            if (read.depth != -1){
                out << "DP:f:" << std::to_string(read.depth);
            }
            out << " LN:i:" << std::to_string(read.sequence_.size()) << "\n";
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
                if (p.second.first[neighbor] != -2){
                    Overlap ov = allOverlaps[contig.neighbors_[neighbor]];
                    Read read = allreads[ov.sequence1];
                    out << read.name << "\t" << p.second.first[neighbor] << "\t" << ov.sequence1 << "\n";
                }
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

void outputGraph_several_clusterings(std::vector<std::vector<int>> &adj,std::vector<std::vector<int>> &clusters, std::vector<bool> &mask, std::string fileOut){

    ofstream out(fileOut);

    out << "nodedef>name VARCHAR,label VARCHAR, ";
    for (auto i = 0 ; i < clusters.size() ; i++){
        out << "cluster_" << i << " VARCHAR, ";
    }
    out << "\n";
    for (auto i = 0 ; i < adj.size() ; i++){
        if (mask[i]){
            out << i << ", " << i << ", ";
            for (auto j = 0 ; j < clusters.size() ; j++){
                out << clusters[j][i] << ", ";
            }
            out << "\n";
        }
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



