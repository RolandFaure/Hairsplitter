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
using std::min;
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

    while(getline(in, line)){

        if (line[0] == 'S'){
            
            string field;
            std::istringstream line2(line);
            int fieldNb = 0;
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
                    indices[nameOfSequence] = sequenceID;
                    sequenceID++;
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
                while(getline(line2, field, '\t')){
                    if (fieldNb == 1){ // name of the sequence1
                        link.neighbor1 = indices[field];
                    }
                    else if (fieldNb == 2){ //here is the sequence
                        if (field == "+"){
                            link.end1 = 1;
                        }
                        else if (field == "-"){
                            link.end1 = 0;
                        }
                        else{
                            cout << "Problem in reading the link : " << line << endl;
                        }
                    }
                    else if (fieldNb == 3){ // name of the sequence1
                        link.neighbor2 = indices[field];
                    }
                    else if (fieldNb == 4){ //here is the sequence
                        if (field == "+"){
                            link.end2 = 0;
                        }
                        else if (field == "-"){
                            link.end2 = 1;
                        }
                        else{
                            cout << "Problem in reading the link : " << line << endl;
                        }
                    }
                    else if (fieldNb == 5){
                        link.CIGAR = field;
                    }

                    fieldNb += 1;
                }
                allLinks.push_back(link);
                allreads[link.neighbor1].add_link(allLinks.size()-1, link.end1);
                allreads[link.neighbor2].add_link(allLinks.size()-1, link.end2);
            }
            catch(...){
                cout << "Problem while reading your GFA file. Please ensure all 'L' lines are below 'S' lines" << endl;
                throw std::invalid_argument("Invalid GFA");
            }
        }

    }
}

//input : a file containing overlaps
//output : the set of all overlaps, updated allreads with overlaps, and optionnaly a list of backbone reads
void parse_PAF(std::string filePAF, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, vector<unsigned long int> &backbones_reads, bool computeBackbones){

    ifstream in(filePAF);
    if (!in){
        cout << "problem reading PAF file " << filePAF << endl;
        throw std::invalid_argument( "Input file '"+filePAF +"' could not be read" );
    }


    vector<bool> backbonesReads (allreads.size(), true); //for now all reads can be backbones read, they will be filtered afterward

    string line;

    while(getline(in, line)){

        string field;
        std::istringstream line2(line);

        unsigned long int sequence1 = 10000;
        string name1;
        int pos1_1= -1;
        int pos1_2= -1;
        int length1= -1;
        unsigned long int sequence2= 100000;
        string name2;
        int pos2_1= -1;
        int pos2_2= -1;
        int length2= -1;
        bool positiveStrand;

        bool allgood = true;
        //now go through the fields of the line
        short fieldnumber = 0;
        while (getline(line2, field, '\t'))
        {
            if (fieldnumber == 0){
                try{
                    sequence1 = indices[field];
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
                    sequence2 = indices[field];
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
            //std::cout << "my field is : " << field << std::endl;
            fieldnumber++;
        }

        if (allgood && fieldnumber > 11 && sequence2 != sequence1){

            //now let's check if this overlap extends to the end of the reads (else, it means the read did not come from the same region)
            bool fullOverlap = false;
            float limit1 = float(pos1_2-pos1_1)/2;
            float limit2 = float(pos2_2-pos2_1)/2;
            if (positiveStrand){
                if (pos1_1  < limit1){ 
                    if (length2-pos2_2 < limit2){
                        fullOverlap = true;
                        pos2_2 = length2;
                    }
                    if (length1-pos1_2 < limit1){
                        fullOverlap = true;
                        pos1_2 = length1;
                    }
                    pos1_1 = 0;
                }
                else if (pos2_1 < limit2){
                    if (length2-pos2_2 < limit2){
                        fullOverlap = true;
                        pos2_2 = length2;
                    }
                    if (length1-pos1_2 < limit1){
                        fullOverlap = true;
                        pos1_2 = length1;
                    }
                    pos2_1 = 0;
                }
            }
            else {
                if (pos1_1 < limit1){
                    if (pos2_1 < limit2){
                        fullOverlap = true;
                        pos2_1 = 0;
                    }
                    if (length1-pos1_2 < limit1){
                        fullOverlap = true;
                        pos1_2 = length1;
                    }
                    pos1_1 = 0;
                }
                else if (length2-pos2_2 < limit2){
                    if (length1-pos1_2 < limit1){
                        fullOverlap = true;
                        pos1_2 = length1;
                    }
                    if (pos2_1 < limit2){
                        fullOverlap = true;
                        pos2_1 = 0;
                    }
                    pos2_2 = length2;
                }
            }
            
            // if (name1 == "2_@DRR198813.15575"){
            //     cout << "reading overlap with 2_@DRR198813.15575 : " << fullOverlap << " " << length1 << " " << pos1_1 << " " << pos1_2 << endl;
            // }
            //add the overlap if it's a full overlap

            if (fullOverlap){

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
                if (backbonesReads[sequence1] && backbonesReads[sequence2] && pos1_2-pos1_1 > 0.5*min(length1, length2)){
                    if (allreads[sequence1].sequence_.size() < allreads[sequence2].sequence_.size()){
                        backbonesReads[sequence1] = false;
                    }
                    else{
                        backbonesReads[sequence2] = false;
                    }
                }
            }
            // else{ //DEBUG
            //     cout << "Not full overlap ! " << allreads[sequence1].name << "," << allreads[sequence2].name << " " << positiveStrand
            //     << " " << pos1_1 << "," << pos1_2 << "," << pos2_1 << "," << pos2_2 << " " << length1-pos1_2 << "," << length2-pos2_2 << endl;
            // }
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

}

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

//input : the list of all reads. Among those, backbone reads are actually contigs

void output_GFA(vector <Read> &allreads, vector<unsigned long int> &backbone_reads, string fileOut, vector<Link> &allLinks)
{
    ofstream out(fileOut);
    for (auto r : backbone_reads){
        Read read = allreads[r];
        out << "S\t"<< read.name << "\t" << read.sequence_.str() << "\n";
    }
    for (auto l : allLinks){
        string end1 = "-";
        if (l.end1 == 1) {end1 = "+";}
        string end2 = "-";
        if (l.end2 == 0) {end2 = "+";}
        out << "L\t" << allreads[l.neighbor1].name << "\t" << end1 << "\t" << allreads[l.neighbor2].name 
            << "\t" << end2 << "\t" << l.CIGAR << "\n";
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
        auto c =  Partition(snps[sus], 0);
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

void outputGraph(std::vector<std::vector<float>> &adj,std::vector<int> &clusters, std::string fileOut){

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








