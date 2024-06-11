//author: Roland Faure
//creation date: 2024-01-17
//aim of the program: scaffold an assembly based on reads. Also break existing contigs if necessary

#include "scaffold.h"
#include "align.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <set>
#include <ctime>

#include "robin_hood.h"
#include "clipp.h" //library to build command line interfaces

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::cout;
using std::endl;
using std::set;
using robin_hood::unordered_map;
using std::min;
using std::max;

// ANSI escape codes for text color
#define RED_TEXT "\033[1;31m"
#define GREEN_TEXT "\033[1;32m"
#define RESET_TEXT "\033[0m"

string version = "0.3.2";
string last_update = "2024-06-03";

vector<string> split(string& s, string& delimiter){
    vector<string> res;
    size_t pos = 0;
    string token = "";

    for (char c : s){
        bool is_delimiter = false;
        for (char d : delimiter){
            if (c == d){
                is_delimiter = true;
                break;
            }
        }
        if (!is_delimiter){
            token += c;
        }
        else {
            res.push_back(token);
            token = "";
        }
    }

    if (token != ""){
        res.push_back(token);
    }
    return res;
}

string reverse_complement(string seq){
    string res = "";
    for (int i = seq.size()-1; i >= 0; i--){
        if (seq[i] == 'A'){
            res += 'T';
        }
        else if (seq[i] == 'T'){
            res += 'A';
        }
        else if (seq[i] == 'C'){
            res += 'G';
        }
        else if (seq[i] == 'G'){
            res += 'C';
        }
        else{
            res += 'N';
        }
    }
    return res;

}

/**
 * @brief From the gaf file, check which reads are aligned end to end, which are partially aligned and which are not aligned at all
 * 
 * @param gaf_file 
 * @param read_file 
 * @param number_of_well_aligned_reads 
 * @param number_of_partially_aligned_reads 
 * @param number_of_unaligned_reads 
 */
void count_unaligned_reads(std::string gaf_file, std::string read_file, int& number_of_unaligned_reads){
    std::ifstream gaf_stream(gaf_file);
    int too_long = 1000; //minimum length of unaligned part to consider it worth reassembling
    number_of_unaligned_reads = 0;

    //first inventoriate the unaligned parts of the reads
    string line;
    string current_read;
    int current_length;
    unordered_map<string, bool> full_alignment_detected;
    unordered_map<string, int> length_of_reads;
    unordered_map<string, int> length_of_alignments_of_reads;
    while (std::getline(gaf_stream, line)){
        string name_of_read;
        int length_of_read;
        int start_of_mapping;
        int end_of_mapping;
        string strand;
        string path;
        int path_length;
        int path_start;
        int path_end;
        string nothing;
        int quality;

        std::istringstream iss(line);
        iss >> name_of_read >> length_of_read >> start_of_mapping >> end_of_mapping >> strand >> path >> path_length >> path_start >> path_end >> nothing >> nothing >> quality;

        //check if the read aligns end to end
        if (quality > 0 && start_of_mapping <= min((double)too_long, 0.2*length_of_read) && end_of_mapping >= max(0.8*length_of_read, (double)length_of_read-too_long)){
            full_alignment_detected[name_of_read] = true;
            length_of_reads[name_of_read] = length_of_read;
            length_of_alignments_of_reads[name_of_read] = end_of_mapping-start_of_mapping;
        }
        else if (quality > 0){
            if (length_of_alignments_of_reads.find(name_of_read) == length_of_alignments_of_reads.end()){
                length_of_alignments_of_reads[name_of_read] = end_of_mapping-start_of_mapping;
                length_of_reads[name_of_read] = length_of_read;
            }
            else{
                length_of_alignments_of_reads[name_of_read] += end_of_mapping-start_of_mapping;
            }
        }

    }
    gaf_stream.close();

    //now go through the reads and count
    std::ifstream read_stream(read_file);
    string read_name;
    while (std::getline(read_stream, line)){
        if (line[0] == '>' || line[0] == '@'){
            read_name = line.substr(1, line.find_first_of(" \t\n")-1);
            if (length_of_reads.find(read_name) == length_of_reads.end()){
                number_of_unaligned_reads++;
            }
        }
    }
}

/**
 * @brief Module to inventoriate all the reads that do not map well on the assembly and re-assemble them using raven
 * 
 * @param gaf_file 
 * @param read_file 
 * @param new_assembly_file 
 * @param new_gaf_file 
 * @param path_minimap2 
 * @param path_miniasm 
 * @param path_minipolish 
 */
void reassemble_unaligned_reads(std::string gaf_file, std::string read_file, std::string old_assembly, std::string new_assembly_file, std::string new_gaf_file, std::string path_raven, std::string path_minigraph, int threads, string tmp_folder){
    //go through the gaf file and retrieve the parts of the reads that are not aligned (at least too_long bp)
    string file_of_unaligned_reads = tmp_folder + "tmp_unaligned_reads.fa";
    string file_of_full_unaligned_reads = tmp_folder + "tmp_full_unaligned_reads.fa";
    std::ifstream gaf_stream(gaf_file);
    int too_long = 1000; //minimum length of unaligned part to consider it worth reassembling

    unordered_map<string,vector<pair<int,int>>> unaligned_reads; //name of the read and the start and end of the unaligned part

    //first inventoriate the unaligned parts of the reads
    string line;
    string current_read;
    int current_length;
    vector<pair<int,int>> aligned_part_of_this_read;
    while (std::getline(gaf_stream, line)){
        string name_of_read;
        int length_of_read;
        int start_of_mapping;
        int end_of_mapping;
        string strand;
        string path;
        int path_length;
        int path_start;
        int path_end;
        string nothing;
        int quality;

        std::istringstream iss(line);
        iss >> name_of_read >> length_of_read >> start_of_mapping >> end_of_mapping >> strand >> path >> path_length >> path_start >> path_end >> nothing >> nothing >> quality;

        if (name_of_read != current_read){
            if (aligned_part_of_this_read.size() > 0){
                //sort the aligned parts of the read by start position
                std::sort(aligned_part_of_this_read.begin(), aligned_part_of_this_read.end(), [](pair<int,int>& a, pair<int,int>& b){return a.first < b.first;});
                //then go through the aligned parts and see if there are unaligned parts
                int start_of_unaligned = 0;
                for (auto part : aligned_part_of_this_read){
                    if (part.first-start_of_unaligned > too_long){
                        if (unaligned_reads.find(current_read) == unaligned_reads.end()){
                            unaligned_reads[current_read] = {};
                        }
                        unaligned_reads[current_read].push_back({start_of_unaligned, part.first});
                    }
                    start_of_unaligned = part.second;
                }
                if (current_length-start_of_unaligned > too_long){
                    if (unaligned_reads.find(current_read) == unaligned_reads.end()){
                        unaligned_reads[current_read] = {};
                    }
                    unaligned_reads[current_read].push_back({start_of_unaligned, current_length});
                }
            }
            current_length = length_of_read;
            current_read = name_of_read;
            aligned_part_of_this_read.clear();
            aligned_part_of_this_read.push_back({start_of_mapping, end_of_mapping});
        }
        else{
            if (quality == 60 && end_of_mapping-start_of_mapping > too_long){
                aligned_part_of_this_read.push_back({start_of_mapping, end_of_mapping});
            }
        }
    }
    gaf_stream.close();

    //now output the unaligned parts of the reads in a file
    std::ifstream read_stream(read_file);
    std::ofstream unaligned_parts_of_reads(file_of_unaligned_reads);
    std::ofstream full_unaligned_parts_of_reads(file_of_full_unaligned_reads);
    string current_seq;
    
    while (std::getline(read_stream, line)){
        if (line[0] == '>' || line[0] == '@'){
            if (current_read != ""){
                if (unaligned_reads.find(current_read) != unaligned_reads.end()){
                    for (auto part : unaligned_reads[current_read]){
                        unaligned_parts_of_reads << ">" << current_read << "_" << part.first << "_" << part.second << endl;
                        unaligned_parts_of_reads << current_seq.substr(part.first, part.second-part.first) << endl;
                    }
                    full_unaligned_parts_of_reads << ">" << current_read << endl;
                    full_unaligned_parts_of_reads << current_seq << endl;
                }
            }
            current_read = line.substr(1, line.find_first_of(" \t\n")-1);
            current_seq = "";
        }
        else{
            current_seq += line;
        }
    }
    if (unaligned_reads.find(current_read) != unaligned_reads.end()){
        for (auto part : unaligned_reads[current_read]){
            unaligned_parts_of_reads << ">" << current_read << "_" << part.first << "_" << part.second << endl;
            unaligned_parts_of_reads << current_seq.substr(part.first, part.second-part.first) << endl;
        }
        full_unaligned_parts_of_reads << ">" << current_read << endl;
        full_unaligned_parts_of_reads << current_seq << endl;
    }
    unaligned_parts_of_reads.close();
    full_unaligned_parts_of_reads.close();
    read_stream.close();

    //then reassemble the unaligned parts of the reads using raven
    string raven_asm = tmp_folder + "tmp_raven_asm.fa";
    string command = path_raven + " " + file_of_unaligned_reads + " -t " + std::to_string(threads) + " > " + raven_asm + " 2> "+ tmp_folder +"raven.log";
    int res = system(command.c_str());
    if (res != 0){
        // cout << RED_TEXT << "ERROR: raven failed to reassemble the unaligned parts of the reads" << RESET_TEXT << endl;
        // cout << "Command was: " << command << endl;
        // exit(1);
    }
    //read the assembly file (fasta format) and add all the new contigs to the new_assembly_file (gfa format)
    std::ifstream raven_asm_stream(raven_asm);
    std::ofstream new_assembly_stream(new_assembly_file);

    string depth = "0";
    while(std::getline(raven_asm_stream, line)){
        if (line[0] == '>'){
            new_assembly_stream << "S\traven_created_" << line.substr(1, line.find_first_of(" \t\n")-1) << "\t";
            //extract the RC:i: field to get the depth of the contig
            int start = line.find("RC:i:")+5;
            int end = line.find_first_of(" \t\n", start);
            depth = line.substr(start, end-start);
        }
        else{
            new_assembly_stream << line << "\tdp:i:" << depth << endl;
        }
    }
    raven_asm_stream.close();
    new_assembly_stream.close();

    //append the old assembly with the new assembly
    system(("cat " + old_assembly + " " + new_assembly_file + " > "+ tmp_folder +"tmp_assembly.gfa").c_str());
    system(("mv "+ tmp_folder +"tmp_assembly.gfa " + new_assembly_file).c_str());

    //map the unaligned reads on the new assembly with minigraph
    command = path_minigraph + " -c -t " + std::to_string(threads) + " " + new_assembly_file + " " + file_of_full_unaligned_reads + " > " + new_gaf_file + " 2> "+ tmp_folder +"minigraph.log";
    res = system(command.c_str());
    if (res != 0){
        cout << RED_TEXT << "ERROR: minigraph failed to map the unaligned reads on the new assembly" << RESET_TEXT << endl;
        exit(1);
    }

    //concat the old and new gaf file, but delete all the re-mapped reads from the old gaf file
    std::ifstream old_gaf_stream(gaf_file);
    string new_new_gaf_file = tmp_folder +"tmp_new_new_gaf_file.gaf";
    std::ofstream new_new_gaf_stream(new_new_gaf_file);

    string old_line;
    while (std::getline(old_gaf_stream, old_line)){
        string name_of_read;
        std::istringstream iss(old_line);
        iss >> name_of_read;
        if (unaligned_reads.find(name_of_read) == unaligned_reads.end()){
            new_new_gaf_stream << old_line << endl;
        }
    }
    old_gaf_stream.close();

    std::ifstream new_gaf_stream(new_gaf_file);
    string new_line;
    while (std::getline(new_gaf_stream, new_line)){
        new_new_gaf_stream << new_line << endl;
    }
    new_gaf_stream.close();
    new_new_gaf_stream.close();

    //mv the new new gaf file to the new gaf file
    system(("mv " + new_new_gaf_file + " " + new_gaf_file).c_str());

    //remove the temporary files
    system(("rm " + raven_asm).c_str());
    system(("rm "+ tmp_folder +"raven.cereal").c_str());
    system(("rm " + file_of_unaligned_reads).c_str());
    system(("rm " + file_of_full_unaligned_reads).c_str());
    system(("rm "+ tmp_folder +"raven.log").c_str());
    system(("rm "+ tmp_folder +"minigraph.log").c_str());

}

/**
 * @brief inventoriate the bridges in a gaf file
 * @details a bridge is defined by a read that maps first to a contig and then to another contig that is not linked to the first one in the assembly graph
 * 
 * @param gaf_file 
 * @param bridges 
 */
void inventoriate_bridges_and_piers(std::string gaf_file, std::vector<Bridge>& bridges, std::vector<Pier>& piers, std::string& assembly_file, 
    int &number_of_well_aligned_reads, int& number_of_partially_aligned_reads, int& number_of_jumping_reads, int &number_of_undetected_duplex_reads, std::set <string>& list_of_duplex_reads){
    
    number_of_well_aligned_reads = 0;
    number_of_partially_aligned_reads = 0;
    number_of_undetected_duplex_reads = 0;
    number_of_jumping_reads = 0;

    std::unordered_map<std::string, long int> length_of_contigs;
    std::ifstream assembly_stream(assembly_file);
    string line;
    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            std::string name_of_contig;
            iss >> token;
            iss >> name_of_contig;
            iss >> token;
            length_of_contigs[name_of_contig] = token.size();
        }
    }
    assembly_stream.close();
    
    std::ifstream gaf_stream(gaf_file);

    unordered_map<string, vector<Mapping>> mappings;
    string delimiter = "<>";

    while (std::getline(gaf_stream, line)){
        string name_of_read;
        int length_of_read;
        int start_of_mapping;
        int end_of_mapping;
        string strand;
        string path;
        int path_length;
        int path_start;
        int path_end;
        string nothing;
        int quality;

        std::istringstream iss(line);
        iss >> name_of_read >> length_of_read >> start_of_mapping >> end_of_mapping >> strand >> path >> path_length >> path_start >> path_end >> nothing >> nothing >> quality;

        if (strand != "+"){
            cout << "ERROR: strand in GAF is not +... panicking !! SOS !!" << endl;
            exit(1);
        }

        int min_length_for_breakpoint = max(0.2*length_of_read, 500.0); //minimum length of overhang to consider a breakpoint AND minimum length of mapping to be confident
        if (quality == 60)
        {        
            vector<string> all_contigs = split(path, delimiter);
            vector<char> orientations;
            //remove the first element
            all_contigs.erase(all_contigs.begin());
            //get the orientations
            for (char c : path){
                if (c == '<' || c == '>'){
                    orientations.push_back(c);
                }
            }

            Mapping mapping;
            mapping.read = name_of_read;
            mapping.length_of_read = length_of_read;

            mapping.contig1 = all_contigs[0];
            if (length_of_contigs.find(all_contigs[0]) == length_of_contigs.end()){
                cout << "ERROR: contig " << all_contigs[0] << " not found in the assembly file " << assembly_file << endl;
                exit(1);
            }
            mapping.position_on_contig1 = path_start;
            mapping.pos_on_read1 = start_of_mapping;
            mapping.orientation_on_contig1 = false;
            mapping.orientation_on_read1 = false;
            if (orientations[0] == '<' ){
                mapping.position_on_contig1 = length_of_contigs[all_contigs[0]] - path_start;
                mapping.orientation_on_contig1 = true;
            }

            mapping.contig2 = all_contigs[all_contigs.size()-1];
            if (length_of_contigs.find(all_contigs[0]) == length_of_contigs.end()){
                cout << "ERROR: contig " << all_contigs[0] << " not found in the assembly file " << assembly_file << endl;
                exit(1);
            }
            mapping.position_on_contig2 = length_of_contigs[all_contigs[all_contigs.size()-1]] - (path_length - path_end);
            mapping.pos_on_read2 = end_of_mapping;
            mapping.orientation_on_contig2 = true;
            mapping.orientation_on_read2 = true;
            if (orientations[orientations.size()-1] == '<'){
                mapping.position_on_contig2 = path_length - path_end;
                mapping.orientation_on_contig2 = false;
            }
            if (mapping.position_on_contig2 < 0){
                cout << line << endl;
                cout << path_length << " " << path_end << " " << length_of_contigs[all_contigs[all_contigs.size()-1]] << endl;
                exit(1);
            }

            mapping.breakpoint1 = false;
            if (start_of_mapping > min_length_for_breakpoint && end_of_mapping-start_of_mapping > min_length_for_breakpoint){ //do not flag a breakpoint if the mapping is too short
                mapping.breakpoint1 = true;
            }
            mapping.breakpoint2 = false;
            if (length_of_read-end_of_mapping > min_length_for_breakpoint && end_of_mapping-start_of_mapping > min_length_for_breakpoint){
                mapping.breakpoint2 = true;
            }

            mappings[name_of_read].push_back(mapping);

            // if ("685d55c2-401c-44ae-9dcd-cfcede23cd75" == name_of_read.substr(0,36) || "8f2b22a9-35df-409d-a939-baa12ef1862f" == name_of_read){
            //     cout << "read " << name_of_read << endl;
            //     cout << "Mapping " << mapping.contig1 << " " << mapping.position_on_contig1 << " " << mapping.orientation_on_contig1 << " " << mapping.breakpoint1 
            //         << " " << mapping.pos_on_read1 << " " << mapping.pos_on_read2 << " " << mapping.breakpoint2 << " " << mapping.contig2 << " " << mapping.position_on_contig2 
            //         << " " << mapping.orientation_on_contig2 << endl;
            // }
        }
    }
    gaf_stream.close();

    //sort the mappings by position on read
    for (auto& mapping : mappings){
        std::sort(mapping.second.begin(), mapping.second.end(), [](Mapping& a, Mapping& b){return a.pos_on_read1 < b.pos_on_read1;});

        bool pier_left_or_right = false;
        bool duplex = false;
        bool bridge_or_not = false;

        //see if there are links between the breakpoints
        for (int b = 0; b < mapping.second.size(); b++){

            //see if there is a pier left of the first breakpoint
            if (b == 0 && mapping.second[b].breakpoint1 == true){
                //check that the full length of the read align after this
                int length_of_alignment = 0;
                for (int c = 0; c < mapping.second.size(); c++){
                    length_of_alignment += mapping.second[c].pos_on_read2 - mapping.second[c].pos_on_read1;
                }
                if (length_of_alignment > (mapping.second[0].length_of_read-mapping.second[0].pos_on_read1)*0.8 && mapping.second[0].pos_on_read1 > max(500.0, 0.5*mapping.second[0].length_of_read)){
                    Pier pier;
                    pier.contig = mapping.second[b].contig1;
                    pier.position = mapping.second[b].position_on_contig1;
                    pier.strand = mapping.second[b].orientation_on_contig1;
                    pier.read_name = mapping.first;
                    pier.pos_read_on_contig = mapping.second[b].pos_on_read1;
                    pier.strand_read = false;
                    piers.push_back(pier);

                    pier_left_or_right = true;
                }
            }

            if (b < mapping.second.size()-1 && mapping.second[b].breakpoint2 == true && mapping.second[b+1].breakpoint1 == true ){ //then the two breakpoints are linked !!

                //if there is only a small overlap between the two mappings, adjust the overlap to remove it
                if (mapping.second[b+1].pos_on_read1-mapping.second[b].pos_on_read2 < 0){
                    int overlap_length = mapping.second[b].pos_on_read2 - mapping.second[b+1].pos_on_read1;
                    if (overlap_length < 0.1*(mapping.second[b].pos_on_read2-mapping.second[b].pos_on_read1) && overlap_length < 0.1*(mapping.second[b+1].pos_on_read2-mapping.second[b+1].pos_on_read1)){
                        //choose which overlap to trim based on which contig is the longest
                        if (length_of_contigs[mapping.second[b].contig2] < length_of_contigs[mapping.second[b+1].contig1]){
                            
                            if (mapping.second[b].orientation_on_contig2 == false){
                                int sliding_length = min(min(overlap_length, (int) length_of_contigs[mapping.second[b].contig2] - mapping.second[b].position_on_contig2-10), mapping.second[b].pos_on_read2);
                                mapping.second[b].pos_on_read2 -= sliding_length;
                                mapping.second[b].position_on_contig2 += sliding_length;
                            }
                            else{
                                int sliding_length = min(min(overlap_length, mapping.second[b].position_on_contig2), mapping.second[b].pos_on_read2);
                                mapping.second[b].pos_on_read2 -= sliding_length;
                                mapping.second[b].position_on_contig2 -= sliding_length;
                            }
                        }
                        else{
                            if (mapping.second[b+1].orientation_on_contig1 == false){
                                int sliding_length = min(overlap_length, (int) length_of_contigs[mapping.second[b+1].contig1]- overlap_length -10);
                                mapping.second[b+1].pos_on_read1 += sliding_length;
                                mapping.second[b+1].position_on_contig1 += sliding_length ;
                            }
                            else{
                                int sliding_length = min(overlap_length, mapping.second[b+1].position_on_contig1-10);
                                mapping.second[b+1].pos_on_read1 += sliding_length;
                                mapping.second[b+1].position_on_contig1 -= sliding_length ;
                            }
                        }
                    }
                }

                if (mapping.second[b+1].pos_on_read1-mapping.second[b].pos_on_read2 >= 0) { //if there is still an overlap, then we can't use this bridge

                    Bridge bridge;
                    bridge.contig1 = mapping.second[b].contig2;
                    bridge.pos_read_on_contig1 = mapping.second[b].pos_on_read2;
                    bridge.position1 = mapping.second[b].position_on_contig2;
                    bridge.strand1 = mapping.second[b].orientation_on_contig2;
                    if (bridge.position1 < 0){
                        cout << "whattttt " << mapping.first << endl;
                        cout << mapping.second[b].position_on_contig2 << endl;
                        exit(1);
                    }

                    bridge.contig2 = mapping.second[b+1].contig1;
                    bridge.pos_read_on_contig2 = mapping.second[b+1].pos_on_read1;
                    bridge.position2 = mapping.second[b+1].position_on_contig1;
                    bridge.strand2 = mapping.second[b+1].orientation_on_contig1;

                    bridge.read_name = mapping.first;

                    //do not consider the reads that do "half-turns" on the contigs, this is an artefact of sequencing (undetected duplex reads)
                    if (bridge.contig1 != bridge.contig2 || abs(bridge.position1-bridge.position2) > 1000 || bridge.strand1 != bridge.strand2){
                        bridges.push_back(bridge);
                        bridge_or_not = true;
                    }
                    else{
                        list_of_duplex_reads.emplace(mapping.first);
                        duplex = true;
                    }

                }
            }

            //see if there is a pier right of the last breakpoint
            if (b == mapping.second.size()-1 && mapping.second[b].breakpoint2 == true){
                //check that the full length of the read align before this
                int length_of_alignment = 0;
                for (int c = 0; c < mapping.second.size(); c++){
                    length_of_alignment += mapping.second[c].pos_on_read2 - mapping.second[c].pos_on_read1;
                }
                if (length_of_alignment > mapping.second[0].pos_on_read2*0.8 && mapping.second[b].pos_on_read2 > max(500.0, 0.5*mapping.second[0].length_of_read)){
                    Pier pier;
                    pier.contig = mapping.second[b].contig2;
                    pier.position = mapping.second[b].position_on_contig2;
                    pier.strand = mapping.second[b].orientation_on_contig2;
                    pier.read_name = mapping.first;
                    pier.pos_read_on_contig = mapping.second[b].pos_on_read2;
                    pier.strand_read = true;
                    piers.push_back(pier);

                    pier_left_or_right = true;
                }
            }
        }

        if (duplex){
            number_of_undetected_duplex_reads++;
        }
        else if (!bridge_or_not && !pier_left_or_right){
            number_of_well_aligned_reads++;
        }
        else if (pier_left_or_right){
            number_of_partially_aligned_reads++;
        }
        else if (bridge_or_not){
            number_of_jumping_reads++;
        }
    }
}

/**
 * @brief agregate bridges into solid bridges
 * 
 * @param bridges 
 * @param solid_bridges 
 * @param aggregative_distance distance in base pairs to agregate bridges
 * @param min_number_of_reads 
 */
void agregate_bridges_and_piers(std::vector<Bridge>& bridges, std::vector<Pier>& piers, std::vector<SolidBridge>& solid_bridges, std::vector<SolidPier>& solid_piers, int aggregative_distance, int min_number_of_reads){
    //agregate bridges into solid bridges
    for (Bridge& bridge : bridges){

        bool found = false;
        for (auto& solid_bridge : solid_bridges){
            if (solid_bridge.contig1 == bridge.contig1 && solid_bridge.contig2 == bridge.contig2 
                && abs(solid_bridge.position1 - bridge.position1) <= aggregative_distance && abs(solid_bridge.position2 - bridge.position2) <= aggregative_distance
                && abs(abs(solid_bridge.pos_read_on_contig1[0]-solid_bridge.pos_read_on_contig2[0]) - abs(bridge.pos_read_on_contig1-bridge.pos_read_on_contig2) ) <= aggregative_distance
                && solid_bridge.strand1 == bridge.strand1 && solid_bridge.strand2 == bridge.strand2){ // the two reads are on the same strand

                if (bridge.pos_read_on_contig2-bridge.pos_read_on_contig1 < solid_bridge.pos_read_on_contig2[0]-solid_bridge.pos_read_on_contig1[0]){
                    solid_bridge.position1 = bridge.position1;
                    solid_bridge.position2 = bridge.position2;

                    solid_bridge.read_names.insert(solid_bridge.read_names.begin(), bridge.read_name);
                    solid_bridge.pos_read_on_contig1.insert(solid_bridge.pos_read_on_contig1.begin(), bridge.pos_read_on_contig1);
                    solid_bridge.pos_read_on_contig2.insert(solid_bridge.pos_read_on_contig2.begin(), bridge.pos_read_on_contig2);
                    solid_bridge.strand.insert(solid_bridge.strand.begin(), true);
                }
                else{
                    solid_bridge.read_names.push_back(bridge.read_name);
                    solid_bridge.pos_read_on_contig1.push_back(bridge.pos_read_on_contig1);
                    solid_bridge.pos_read_on_contig2.push_back(bridge.pos_read_on_contig2);
                    solid_bridge.strand.push_back(true);
                }
                found = true;
                break;
            }
            else if (solid_bridge.contig1 == bridge.contig2 && solid_bridge.contig2 == bridge.contig1 
                && abs(solid_bridge.position1 - bridge.position2) <= aggregative_distance && abs(solid_bridge.position2 - bridge.position1) <= aggregative_distance
                && abs(abs(solid_bridge.pos_read_on_contig1[0]-solid_bridge.pos_read_on_contig2[0]) - abs(bridge.pos_read_on_contig1-bridge.pos_read_on_contig2) ) <= aggregative_distance
                && solid_bridge.strand1 == bridge.strand2 && solid_bridge.strand2 == bridge.strand1){ // the two reads are on opposite strands

                if (bridge.pos_read_on_contig2-bridge.pos_read_on_contig1 < solid_bridge.pos_read_on_contig2[0]-solid_bridge.pos_read_on_contig1[0]){
                    solid_bridge.position1 = bridge.position2;
                    solid_bridge.position2 = bridge.position1;

                    solid_bridge.read_names.insert(solid_bridge.read_names.begin(), bridge.read_name);
                    solid_bridge.pos_read_on_contig1.insert(solid_bridge.pos_read_on_contig1.begin(), bridge.pos_read_on_contig1);
                    solid_bridge.pos_read_on_contig2.insert(solid_bridge.pos_read_on_contig2.begin(), bridge.pos_read_on_contig2);
                    solid_bridge.strand.insert(solid_bridge.strand.begin(), false);
                }
                else{
                    solid_bridge.read_names.push_back(bridge.read_name);
                    solid_bridge.pos_read_on_contig1.push_back(bridge.pos_read_on_contig1);
                    solid_bridge.pos_read_on_contig2.push_back(bridge.pos_read_on_contig2);
                    solid_bridge.strand.push_back(false);
                }
                found = true;
                break;
            }
        }
        if (!found){
            SolidBridge solid_bridge;
            solid_bridge.contig1 = bridge.contig1;
            solid_bridge.contig2 = bridge.contig2;
            solid_bridge.position1 = bridge.position1;
            solid_bridge.position2 = bridge.position2;
            solid_bridge.strand1 = bridge.strand1;
            solid_bridge.strand2 = bridge.strand2;
            solid_bridge.read_names.push_back(bridge.read_name);
            solid_bridge.pos_read_on_contig1.push_back(bridge.pos_read_on_contig1);
            solid_bridge.pos_read_on_contig2.push_back(bridge.pos_read_on_contig2);
            solid_bridge.strand.push_back(true);
            solid_bridges.push_back(solid_bridge);
        }
    }

    //keep only the solid bridges that have at least min_number_of_reads reads
    std::vector<SolidBridge> solid_bridges_kept;
    for (auto solid_bridge : solid_bridges){
        if (solid_bridge.read_names.size() >= min_number_of_reads){
            solid_bridges_kept.push_back(solid_bridge);
        }
    }
    solid_bridges = solid_bridges_kept;

    //agregate piers into solid piers
    for (Pier& pier : piers){
        bool found = false;
        for (auto& solid_pier : solid_piers){
            if (solid_pier.contig == pier.contig && abs(solid_pier.position - pier.position) <= aggregative_distance
                && solid_pier.strand == pier.strand){
                solid_pier.read_names.push_back(pier.read_name);
                solid_pier.pos_read_on_contig.push_back(pier.pos_read_on_contig);
                solid_pier.strands_read.push_back(pier.strand_read);
                found = true;
                // if (abs(pier.position - 2026206)  < 1000){
                //     cout << "Adding pier " << pier.read_name << " to solid pier " << solid_pier.read_names[0] << " at position " << pier.position << " and solid pier position " << solid_pier.position << endl;
                // }
                break;
            }
        }
        if (!found){
            SolidPier solid_pier;
            solid_pier.contig = pier.contig;
            solid_pier.position = pier.position;
            solid_pier.strand = pier.strand;
            solid_pier.read_names.push_back(pier.read_name);
            solid_pier.pos_read_on_contig.push_back(pier.pos_read_on_contig);
            solid_pier.strands_read.push_back(pier.strand_read);
            solid_piers.push_back(solid_pier);
        }
    }

    //now keep only the solid piers that have at least min_number_of_reads reads AND that are not just one end of a bridge but not long enough to be a bridge
    std::vector<SolidPier> solid_piers_kept;
    for (auto solid_pier : solid_piers){
        if (solid_pier.read_names.size() >= min_number_of_reads){
            bool is_pier = true;
            for (auto solid_bridge : solid_bridges){
                if (solid_bridge.contig1 == solid_pier.contig || solid_bridge.contig2 == solid_pier.contig){
                    if ((abs(solid_bridge.position1 - solid_pier.position) <= aggregative_distance  && solid_bridge.strand1 == solid_pier.strand)
                        || (abs(solid_bridge.position2 - solid_pier.position) <= aggregative_distance && solid_bridge.strand2 == solid_pier.strand)){
                        is_pier = false; //the pier is just one end of a bridge
                        break;
                    }
                }
            }
            if (is_pier){
                solid_piers_kept.push_back(solid_pier);
            }
        }
    }
    solid_piers = solid_piers_kept;
}

/**
 * @brief transform solid bridges into links, ie with precise coordinates and if needed gap filling
 * 
 * @param solid_bridges 
 * @param read_file 
 * @param links result of the transformation
 
 */
void transform_bridges_in_links(std::vector<SolidBridge>& solid_bridges, std::string read_file, std::string assembly_file, std::vector<Link>& links,
    std::string& path_minimap2, std::string& path_racon, std::string& tmp_folder){

    //index the positions of the reads in the read file
    robin_hood::unordered_map<std::string, long int> read_positions;
    std::ifstream read_stream(read_file);
    string line;
    long int position = 0;
    auto pos = read_stream.tellg();
    while (std::getline(read_stream, line)){
        pos = read_stream.tellg();
        if (line[0] == '>' || line[0] == '@'){
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            read_positions[token.substr(1)] = pos;
        }
    }
    read_stream.close();

    //index the positions of the contigs in the assembly file
    robin_hood::unordered_map<std::string, long int> contig_positions;
    std::ifstream assembly_stream(assembly_file);
    pos = assembly_stream.tellg();
    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            iss >> token;
            contig_positions[token] = pos;
        }
        pos = assembly_stream.tellg();
    }

    //then transform the solid bridges into links
    for (auto solid_bridge : solid_bridges){
        Link link;
        link.contig1 = solid_bridge.contig1;
        link.contig2 = solid_bridge.contig2;
        link.strand1 = solid_bridge.strand1;
        link.strand2 = solid_bridge.strand2;
        link.coverage = solid_bridge.read_names.size();

        //compute the sequence left_of_the_junction + read_chunk + sequence_right_of_the_junction and see if we adjust the coordinates
        //retrieve the sequence of the read in the junction
        std::ifstream read_stream(read_file);
        read_stream.seekg(read_positions[solid_bridge.read_names[0]]);
        std::string gap_filling_seq;
        std::string seq_with_overhangs;

        std::getline(read_stream, line);
        if (solid_bridge.strand[0] == true){
            gap_filling_seq = line.substr(solid_bridge.pos_read_on_contig1[0], solid_bridge.pos_read_on_contig2[0]-solid_bridge.pos_read_on_contig1[0]);
        }
        else {
            gap_filling_seq = reverse_complement(line.substr(solid_bridge.pos_read_on_contig1[0], solid_bridge.pos_read_on_contig2[0]-solid_bridge.pos_read_on_contig1[0]));
        }
        read_stream.close();

        //retrieve the sequence of the contig on the left of the junction
        std::ifstream assembly_stream(assembly_file);
        std::string contig1_sequence, contig2_sequence;

        assembly_stream.seekg(contig_positions[solid_bridge.contig1]);
        std::getline(assembly_stream, line);
        string seq_1, seq_2, nothing;
        std::istringstream iss(line);
        iss >> nothing >> nothing >> seq_1;
        if (solid_bridge.strand1 == true){
            int overhang_1 = min(solid_bridge.position1 , 200);
            contig1_sequence = seq_1.substr(solid_bridge.position1-overhang_1, overhang_1);
        }
        else{
            int overhang_1 = min(200, (int) seq_1.size()-solid_bridge.position1);
            contig1_sequence = reverse_complement(seq_1.substr(solid_bridge.position1, overhang_1));
        }

        //retrieve the sequence of the contig on the right of the junction
        assembly_stream.seekg(contig_positions[solid_bridge.contig2]);
        std::getline(assembly_stream, line);
        std::istringstream iss2(line);
        iss2 >> nothing >> nothing >> seq_2;
        if (solid_bridge.strand2 == true){
            int overhang_2 = min(200, solid_bridge.position2);
            contig2_sequence = reverse_complement(seq_2.substr(solid_bridge.position2-overhang_2, overhang_2));
        }
        else{
            int overhang_2 = min((int) seq_2.size() - solid_bridge.position2, 200);
            contig2_sequence = seq_2.substr(solid_bridge.position2, overhang_2);
        }
        assembly_stream.close();

        seq_with_overhangs = contig1_sequence + gap_filling_seq + contig2_sequence;


        //now retrieve all the reads that make the link to polish the gap filling sequence
        std::vector<std::string> reads;
        for (auto read_num = 1 ; read_num < solid_bridge.read_names.size() ; read_num++){
            string read_name = solid_bridge.read_names[read_num];
            std::ifstream read_stream(read_file);
            read_stream.seekg(read_positions[read_name]);
            std::getline(read_stream, line);
            
            //retrieve the gap filling sequence only, with 100 bp on each side
            int start = max(solid_bridge.pos_read_on_contig1[read_num]-300, 0);
            int end = min(solid_bridge.pos_read_on_contig2[read_num]+300, (int) line.size());
            if (solid_bridge.strand[read_num] == true){
                reads.push_back(line.substr(start, end-start));
            }
            else {
                reads.push_back(reverse_complement(line.substr(start, end-start)));
            }

            read_stream.close();
        }

        //polish the gap filling sequence
        string polished_gap_filling_seq = polish(seq_with_overhangs, reads, path_minimap2, path_racon, tmp_folder);

        //now align the gap filling sequence to the contigs left and right to see where the junction is exactly and what to put in between
        string contig1_sequence2, contig2_sequence2; //sequences on the other side of the junction compared to before, i.e. towards the inside of the junction
        if (solid_bridge.strand1 == true){
            int overhang_1 = min(seq_1.size() - solid_bridge.position1 , gap_filling_seq.size());
            contig1_sequence2 = seq_1.substr(solid_bridge.position1, overhang_1);
        }
        else{
            int overhang_1 = min( (int) gap_filling_seq.size(), solid_bridge.position1);
            contig1_sequence2 = reverse_complement(seq_1.substr(solid_bridge.position1-overhang_1, overhang_1));
        }
        string cigar = align(contig1_sequence2, 0, contig1_sequence2.size(), gap_filling_seq, 0, gap_filling_seq.size());

        int end_of_match_gap_fill = 0;
        int end_of_match_contig1 = 0;
        int pos_c1 = 0;
        int pos_gf = 0;
        int consecutive_matches = 5;
        int num_indel = 0;
        int num_matches = 0;

        for (auto c = 0 ; c < cigar.size() ; c++){
            if (cigar[c] == '='){
                consecutive_matches++;
                pos_c1++;
                pos_gf++;
                num_matches++;
            }
            else{
                consecutive_matches = 0;
                if (cigar[c] == 'I'){
                    pos_gf++;
                }
                else if (cigar[c] == 'D'){
                    pos_c1++;
                }
                else if (cigar[c] == 'X'){
                    pos_c1++;
                    pos_gf++;
                }
                num_indel++;
            }
            if (consecutive_matches >= 5 && num_indel <= 0.2*num_matches){
                end_of_match_gap_fill = pos_gf;
                end_of_match_contig1 = pos_c1;
            }
        }

        //slide the gap filling sequence and the contig to the right by end_of_match_gap_fill and end_of_match_contig1
        if (solid_bridge.strand1 == true){
            link.position1 = solid_bridge.position1 + end_of_match_contig1;
        }
        else{
            link.position1 = solid_bridge.position1 - end_of_match_contig1;
        }
        end_of_match_gap_fill = min(end_of_match_gap_fill, (int) gap_filling_seq.size());
        gap_filling_seq = gap_filling_seq.substr(end_of_match_gap_fill, gap_filling_seq.size()-end_of_match_gap_fill);

        if (solid_bridge.strand2 == true){
            int overhang_2 = min(seq_2.size() - solid_bridge.position2 , gap_filling_seq.size());
            contig2_sequence2 = reverse_complement(seq_2.substr(solid_bridge.position2, overhang_2));
        }
        else{
            int overhang_2 = min( (int) gap_filling_seq.size(), solid_bridge.position2);
            contig2_sequence2 = seq_2.substr(solid_bridge.position2-overhang_2, overhang_2);
        }

        //align to gap filling seq but with reverse complement
        string rcc2 = reverse_complement(contig2_sequence2);
        string rcgf = reverse_complement(gap_filling_seq);
        string cigar2 = align(rcc2, 0, contig2_sequence2.size(), rcgf, 0, gap_filling_seq.size());
        
        int end_of_match_gap_fill2 = 0;
        int end_of_match_contig2 = 0;
        int pos_c2 = 0;
        int pos_gf2 = 0;
        consecutive_matches = 5;
        num_indel = 0;
        num_matches = 0;

        for (auto c = 0 ; c < cigar2.size() ; c++){
            if (cigar2[c] == '='){
                consecutive_matches++;
                pos_c2++;
                pos_gf2++;
                num_matches++;
            }
            else{
                consecutive_matches = 0;
                if (cigar2[c] == 'I'){
                    pos_gf2++;
                }
                else if (cigar2[c] == 'D'){
                    pos_c2++;
                }
                else if (cigar2[c] == 'X'){
                    pos_c2++;
                    pos_gf2++;
                }
                num_indel++;
            }
            if (consecutive_matches >= 5 && num_indel <= 0.2*num_matches){
                end_of_match_gap_fill2 = pos_gf2;
                end_of_match_contig2 = pos_c2;
            }
        }

        //slide the gap filling sequence and the contig to the right by end_of_match_gap_fill and end_of_match_contig2
        if (solid_bridge.strand2 == true){
            link.position2 = solid_bridge.position2 + end_of_match_contig2;
        }
        else{
            link.position2 = solid_bridge.position2 - end_of_match_contig2;
        }
        end_of_match_gap_fill2 = min(end_of_match_gap_fill2, (int) gap_filling_seq.size());
        gap_filling_seq = gap_filling_seq.substr(end_of_match_gap_fill2, gap_filling_seq.size()-end_of_match_gap_fill2);

        link.extra_sequence = gap_filling_seq;
        links.push_back(link);

        // cout << "solid bridge: " << endl;
        // cout << "contig1: " << solid_bridge.contig1 << " " << solid_bridge.position1 << " " << solid_bridge.strand1 << " " << seq_1.size() << endl;
        // cout << "CIGAR of the gap filling seq on the right of contig 1: " << cigar << endl;
        // cout << "I would slide the gap filling seq by " << end_of_match_gap_fill << " bp and the contig by " << end_of_match_contig1 << endl;
        // cout << "contig 2: " << solid_bridge.contig2 << " " << solid_bridge.position2 << " " << solid_bridge.strand2 << " " << seq_2.size() << endl;
        // cout << "CIGAR of the gap filling seq on the left of contig 2: " << cigar2 << endl;
        // cout << "I would slide the gap filling seq by " << end_of_match_gap_fill2 << " bp and the contig by " << end_of_match_contig2 << endl;
        // cout << endl;
        // exit(1);
    }

}

/**
 * @brief 
 * 
 * @param solid_piers 
 * @param read_file 
 * @param assembly_file 
 * @param end_contigs 
 * @param path_minimap2 
 * @param path_racon 
 */
void build_piers(std::vector<SolidPier>& solid_piers, std::string read_file, std::string assembly_file, std::vector<End_contig>& end_contigs,
    std::string& path_minimap2, std::string& path_racon, string& tmp_folder){
    
    //index the positions of the reads in the read file
    robin_hood::unordered_map<std::string, long int> read_positions;
    std::ifstream read_stream(read_file);
    string line;
    long int position = 0;
    auto pos = read_stream.tellg();
    while (std::getline(read_stream, line)){
        pos = read_stream.tellg();
        if (line[0] == '>' || line[0] == '@'){
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            read_positions[token.substr(1)] = pos;
        }
    }
    read_stream.close();

    //index the positions of the contigs in the assembly file
    robin_hood::unordered_map<std::string, long int> contig_positions;
    std::ifstream assembly_stream(assembly_file);
    pos = assembly_stream.tellg();
    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            iss >> token;
            contig_positions[token] = pos;
        }
        pos = assembly_stream.tellg();
    }

    //now create extra sequences for the piers and link them
    for (auto solid_pier : solid_piers){

        //retrieve the sequence of the contig on the left of the junction
        std::ifstream assembly_stream(assembly_file);
        std::string contig1_sequence, contig2_sequence;

        assembly_stream.seekg(contig_positions[solid_pier.contig]);
        std::getline(assembly_stream, line);
        string seq_1, seq_2, nothing;
        std::istringstream iss(line);
        iss >> nothing >> nothing >> seq_1;
        if (solid_pier.strand == true){
            int overhang_1 = min(solid_pier.position , 200);
            contig1_sequence = seq_1.substr(solid_pier.position-overhang_1, overhang_1);
        }
        else{
            int overhang_1 = min(200, (int) seq_1.size()-solid_pier.position);
            contig1_sequence = reverse_complement(seq_1.substr(solid_pier.position, overhang_1));
        }

        //retrieve the longest sequence of the reads that make the pier
        int longest_length = 0;
        string seq_to_polish;
        vector<string> reads;
        bool strand_of_longest_seq;
        for (auto read_num = 0 ; read_num < solid_pier.read_names.size() ; read_num++){
            string read_name = solid_pier.read_names[read_num];
            std::ifstream read_stream(read_file);
            read_stream.seekg(read_positions[read_name]);
            std::getline(read_stream, line);
            if (solid_pier.strands_read[read_num] && line.size()-solid_pier.pos_read_on_contig[read_num] > longest_length){
                longest_length = line.size()-solid_pier.pos_read_on_contig[read_num];
                seq_to_polish = line.substr(solid_pier.pos_read_on_contig[read_num], line.size()-solid_pier.pos_read_on_contig[read_num]);
                strand_of_longest_seq = solid_pier.strands_read[read_num];
            }
            else if (solid_pier.strands_read[read_num] == false && solid_pier.pos_read_on_contig[read_num] > longest_length){
                longest_length = solid_pier.pos_read_on_contig[read_num];
                seq_to_polish = reverse_complement(line.substr(0, solid_pier.pos_read_on_contig[read_num]));
                strand_of_longest_seq = solid_pier.strands_read[read_num];
            }
            reads.push_back(line);
            read_stream.close();
        }
        //polish the sequence
        string polished_seq = polish(seq_to_polish, reads, path_minimap2, path_racon, tmp_folder);

        //now align the polished sequence to the contig to see where the junction is exactly
        string contig1_sequence2; //sequences on the other side of the junction compared to before, i.e. towards the inside of the junction
        if (solid_pier.strand == true){
            int overhang_1 = min(seq_1.size() - solid_pier.position , polished_seq.size());
            contig1_sequence2 = seq_1.substr(solid_pier.position, overhang_1);
        }
        else{
            int overhang_1 = min( (int) polished_seq.size(), solid_pier.position);
            contig1_sequence2 = reverse_complement(seq_1.substr(solid_pier.position-overhang_1, overhang_1));
        }
        string cigar = align(contig1_sequence2, 0, contig1_sequence2.size(), polished_seq, 0, polished_seq.size());

        int end_of_match_polished_seq = 0;
        int end_of_match_contig1 = 0;
        int pos_c1 = 0;
        int pos_ps = 0;
        int consecutive_matches = 5;
        int num_indel = 0;
        int num_matches = 0;

        for (auto c = 0 ; c < cigar.size() ; c++){
            if (cigar[c] == '='){
                consecutive_matches++;
                pos_c1++;
                pos_ps++;
                num_matches++;
            }
            else{
                consecutive_matches = 0;
                if (cigar[c] == 'I'){
                    pos_ps++;
                }
                else if (cigar[c] == 'D'){
                    pos_c1++;
                }
                else if (cigar[c] == 'X'){
                    pos_c1++;
                    pos_ps++;
                }
                num_indel++;
            }
            if (consecutive_matches >= 5 && num_indel <= 0.2*num_matches){
                end_of_match_polished_seq = pos_ps;
                end_of_match_contig1 = pos_c1;
            }
        }

        //slide the polished sequence and the contig to the right by end_of_match_polished_seq and end_of_match_contig1
        End_contig end_contig;
        end_contig.contig = solid_pier.contig;
        if (solid_pier.strand == true){
            end_contig.position = solid_pier.position + end_of_match_contig1;
        }
        else{
            end_contig.position = solid_pier.position - end_of_match_contig1;
        }
        
        end_of_match_polished_seq = min(end_of_match_polished_seq, (int) polished_seq.size()); //should never happen but you never know
        end_contig.extra_sequence = polished_seq.substr(end_of_match_polished_seq, polished_seq.size()-end_of_match_polished_seq);
        end_contig.strand = solid_pier.strand;
        end_contig.strand_read = strand_of_longest_seq;
        
        end_contigs.push_back(end_contig);

    }

}


/**
 * @brief create a gfa file from the assembly and the links
 * 
 * @param input_assembly 
 * @param output_assembly 
 * @param links 
 */
void create_gfa(std::string& input_assembly, std::string& output_assembly, std::vector<Link>& links, std::vector<End_contig>& end_contigs){
    
    unordered_map<string, vector<int>> breakpoints_in_contigs;
    unordered_map<string, int> length_of_contigs;

    //go through the assembly file to get the length of the contigs
    std::ifstream assembly_stream(input_assembly);
    string line;
    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            std::string name_of_contig;
            iss >> token;
            iss >> name_of_contig;
            string sequence;
            iss >> sequence;

            length_of_contigs[name_of_contig] = sequence.size();
        }
    }
    assembly_stream.close();

    //go through the links and inventoriate the breakpoints
    for (auto link : links){
        string contig1 = link.contig1;
        int position1 = link.position1;

        //insert the breakpoint in the contig
        if (position1 > 0 && position1 < length_of_contigs[contig1]){
            if (breakpoints_in_contigs.find(contig1) == breakpoints_in_contigs.end()){
                breakpoints_in_contigs[contig1] = {position1};
            }
            else{ //insert the breakpoint at the right position
                int bp = 0;
                while (bp < breakpoints_in_contigs[contig1].size() && breakpoints_in_contigs[contig1][bp] < position1){
                    bp++;
                }
                if (bp == breakpoints_in_contigs[contig1].size()){
                    breakpoints_in_contigs[contig1].push_back(position1);
                }
                else if (breakpoints_in_contigs[contig1][bp] > position1){
                    breakpoints_in_contigs[contig1].insert(breakpoints_in_contigs[contig1].begin()+bp, position1);
                }
            }
        }

        string contig2 = link.contig2;
        int position2 = link.position2;
        if (position2 > 0 && position2 < length_of_contigs[contig2]){
            if (breakpoints_in_contigs.find(contig2) == breakpoints_in_contigs.end()){
                breakpoints_in_contigs[contig2] = {position2};
            }
            else{ //insert the breakpoint at the right position
                int bp = 0;
                while (bp < breakpoints_in_contigs[contig2].size() && breakpoints_in_contigs[contig2][bp] < position2){
                    bp++;
                }
                if (bp == breakpoints_in_contigs[contig2].size()){
                    breakpoints_in_contigs[contig2].push_back(position2);
                }
                else if (breakpoints_in_contigs[contig2][bp] > position2){
                    breakpoints_in_contigs[contig2].insert(breakpoints_in_contigs[contig2].begin()+bp, position2);
                }
            }
        }
    }

    //go through the end contigs and inventoriate the breakpoints
    for (auto end_contig : end_contigs){
        string contig = end_contig.contig;
        int position = end_contig.position;

        //insert the breakpoint in the contig
        if (position > 0 && position < length_of_contigs[contig]){
            if (breakpoints_in_contigs.find(contig) == breakpoints_in_contigs.end()){
                breakpoints_in_contigs[contig] = {position};
            }
            else{ //insert the breakpoint at the right position
                int bp = 0;
                while (bp < breakpoints_in_contigs[contig].size() && breakpoints_in_contigs[contig][bp] < position){
                    bp++;
                }
                if (bp == breakpoints_in_contigs[contig].size()){
                    breakpoints_in_contigs[contig].push_back(position);
                }
                else if (breakpoints_in_contigs[contig][bp] > position){
                    breakpoints_in_contigs[contig].insert(breakpoints_in_contigs[contig].begin()+bp, position);
                }
            }
        }
    }

    //now go through the assembly and output broken contigs at breakpoints
    vector<string> L_lines;
    assembly_stream = std::ifstream (input_assembly);
    std::ofstream output_stream(output_assembly);
    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            std::string name_of_contig;
            iss >> token;
            iss >> name_of_contig;
            string sequence;
            iss >> sequence;

            vector<string> tags;
            while (iss >> token){
                tags.push_back(token);
            }

            length_of_contigs[name_of_contig] = sequence.size();

            if (breakpoints_in_contigs.find(name_of_contig) != breakpoints_in_contigs.end()){
                int last_bp = 0;
                int bp_index = 0;
                for (auto bp : breakpoints_in_contigs[name_of_contig]){
                    output_stream << "S\t" << name_of_contig << "_" << last_bp << "_" << bp << "\t" << sequence.substr(last_bp, bp-last_bp);
                    for (auto tag : tags){
                        output_stream << "\t" << tag;
                    }
                    output_stream << endl;
                    int next_bp = length_of_contigs[name_of_contig];
                    if (bp_index < breakpoints_in_contigs[name_of_contig].size()-1){
                        next_bp = breakpoints_in_contigs[name_of_contig][bp_index+1];
                    }
                    //now create the L line
                    string L_line = "L\t" + name_of_contig + "_" + std::to_string(last_bp) + "_" + std::to_string(bp) + "\t+\t" + name_of_contig + "_" + std::to_string(bp) + "_" + std::to_string(next_bp) + "\t+\t0M";
                    L_lines.push_back(L_line);
                    last_bp = bp;
                    bp_index++;
                }
                output_stream << "S\t" << name_of_contig << "_" << last_bp << "_" << sequence.size() << "\t" << sequence.substr(last_bp, sequence.size()-last_bp);
                for (auto tag : tags){
                    output_stream << "\t" << tag;
                }
                output_stream << endl;
            }
            else{
                output_stream << "S\t" << name_of_contig << "_0_" << sequence.size() << "\t" << sequence;
                for (auto tag : tags){
                    output_stream << "\t" << tag;
                }
                output_stream << endl;
            }
        }
        else if (line[0] == 'L'){
            std::istringstream iss(line);
            std::string nothing, contig1, orientation1, contig2, orientation2, cigar;
            iss >> nothing >> contig1 >> orientation1 >> contig2 >> orientation2 >> cigar;

            //now find the new name of contig 1 and contig 2
            if (breakpoints_in_contigs.find(contig1) != breakpoints_in_contigs.end()){
                if (orientation1 == "+"){
                    contig1 = contig1 + "_" + std::to_string(breakpoints_in_contigs[contig1][breakpoints_in_contigs[contig1].size()-1]) 
                        + "_" + std::to_string(length_of_contigs[contig1]);
                }
                else{
                    contig1 = contig1 + "_0_" + std::to_string(breakpoints_in_contigs[contig1][0]);
                }
            }
            else{
                contig1 = contig1 + "_0_" + std::to_string(length_of_contigs[contig1]);
            }

            if (breakpoints_in_contigs.find(contig2) != breakpoints_in_contigs.end()){
                if (orientation2 == "+"){
                    contig2 = contig2 + "_0_" + std::to_string(breakpoints_in_contigs[contig2][0]);
                }
                else{
                    contig2 = contig2 + "_" + std::to_string(breakpoints_in_contigs[contig2][breakpoints_in_contigs[contig2].size()-1]) 
                        + "_" + std::to_string(length_of_contigs[contig2]);
                }
            }
            else {
                contig2 = contig2 + "_0_" + std::to_string(length_of_contigs[contig2]);
            }

            L_lines.push_back("L\t" + contig1 + "\t" + orientation1 + "\t" + contig2 + "\t" + orientation2 + "\t" + cigar);
        }
    }

    //now go through the links
    for (auto link : links){

        string contig1 = link.contig1;
        int position1 = link.position1;
        char strand1 = "-+"[link.strand1];

        string contig2 = link.contig2;
        int position2 = link.position2;
        char strand2 = "+-"[link.strand2];

        string name1;
        int pos_before = 0;
        int pos_after = length_of_contigs[contig1];
        for (int b = 0 ; b < breakpoints_in_contigs[contig1].size() ; b++){
            if (breakpoints_in_contigs[contig1][b] == position1){
                if (b < breakpoints_in_contigs[contig1].size()-1){
                    pos_after = breakpoints_in_contigs[contig1][b+1];
                }
                if (b > 0){
                    pos_before = breakpoints_in_contigs[contig1][b-1];
                }
                break;
            }
        }
        if (position1 == 0 && breakpoints_in_contigs[contig1].size() > 0){
            pos_after = breakpoints_in_contigs[contig1][0];
        }
        if (position1 == length_of_contigs[contig1] && breakpoints_in_contigs[contig1].size() > 0){
            pos_before = breakpoints_in_contigs[contig1][breakpoints_in_contigs[contig1].size()-1];
        }
        if (strand1 == '-'){
            name1 = contig1 + "_" + std::to_string(position1) + "_" + std::to_string(pos_after);
        }
        else{
            name1 = contig1 + "_" + std::to_string(pos_before) + "_" + std::to_string(position1);
        }

        string name2;
        pos_before = 0;
        pos_after = length_of_contigs[contig2];
        for (int b = 0 ; b < breakpoints_in_contigs[contig2].size() ; b++){
            if (breakpoints_in_contigs[contig2][b] == position2){
                if (b < breakpoints_in_contigs[contig2].size()-1){
                    pos_after = breakpoints_in_contigs[contig2][b+1];
                }
                if (b > 0){
                    pos_before = breakpoints_in_contigs[contig2][b-1];
                }
                break;
            }
        }
        if (position2 == 0 && breakpoints_in_contigs[contig2].size() > 0){
            pos_after = breakpoints_in_contigs[contig2][0];
        }
        if (position2 == length_of_contigs[contig2] && breakpoints_in_contigs[contig2].size() > 0){
            pos_before = breakpoints_in_contigs[contig2][breakpoints_in_contigs[contig2].size()-1];
        }
        
        if (strand2 == '+'){
            name2 = contig2 + "_" + std::to_string(position2) + "_" + std::to_string(pos_after);
        }
        else{
            name2 = contig2 + "_" + std::to_string(pos_before) + "_" + std::to_string(position2);
        }

        //if necessary, create a new contig for the link
        if (link.extra_sequence.size() > 0){
            string name_join = "join_"+ contig1 + "_" + std::to_string(position1) + strand1 + "_" + contig2 + "_" + std::to_string(position2) + strand2+"|";
            output_stream << "S\t" << name_join << "\t" << link.extra_sequence << "\tdp:i:" << link.coverage  << endl;

            //now create the two L lines
            string L1 = "L\t" + name1 + "\t" + strand1 + "\t" + name_join + "\t+\t0M";
            string L2 = "L\t" + name_join + "\t+\t" + name2 + "\t" + strand2 + "\t0M";
            L_lines.push_back(L1);
            L_lines.push_back(L2);
        }
        else{
            string L1 = "L\t" + name1 + "\t" + strand1 + "\t" + name2 + "\t" + strand2 + "\t0M";
            L_lines.push_back(L1);
        }
    }

    //now go through the end contigs
    for (auto end_contig: end_contigs){
        string contig = end_contig.contig;
        int position = end_contig.position;
        char strand = "-+"[end_contig.strand];

        string name;
        int pos_before = 0;
        int pos_after = length_of_contigs[contig];
        for (int b = 0 ; b < breakpoints_in_contigs[contig].size() ; b++){
            if (breakpoints_in_contigs[contig][b] == position){
                if (b < breakpoints_in_contigs[contig].size()-1){
                    pos_after = breakpoints_in_contigs[contig][b+1];
                }
                if (b > 0){
                    pos_before = breakpoints_in_contigs[contig][b-1];
                }
                break;
            }
        }
        if (position == 0 && breakpoints_in_contigs[contig].size() > 0){
            pos_after = breakpoints_in_contigs[contig][0];
        }
        if (position == length_of_contigs[contig] && breakpoints_in_contigs[contig].size() > 0){
            pos_before = breakpoints_in_contigs[contig][breakpoints_in_contigs[contig].size()-1];
        }
        if (strand == '-'){
            name = contig + "_" + std::to_string(position) + "_" + std::to_string(pos_after);
        }
        else{
            name = contig + "_" + std::to_string(pos_before) + "_" + std::to_string(position);
        }

        //create the extra contig
        string name_pier = "extend_"+ contig + "_" + std::to_string(position) + strand + "|";
        output_stream << "S\t" << name_pier << "\t" << end_contig.extra_sequence << "\tdp:i:" << end_contig.coverage  << endl;

        //now create the L line from the contig to the extra contig
        char orientationContig = "-+"[end_contig.strand];
        char orientationPier = '+'; //that's because we took the reverse complement if needed when polishing the pier

        string L1 = "L\t" + name + "\t" + orientationContig + "\t" + name_pier + "\t" + orientationPier + "\t0M";
        L_lines.push_back(L1);
    }

    //print all the L lines
    for (auto L_line : L_lines){
        output_stream << L_line << endl;
    }
}

/**
 * @brief Shave the graph of small dead ends and pop small bubbles, that probably result from polishing errors
 * 
 * @param input_file 
 * @param output_file 
 * @param max_length 
 */
void shave_and_pop(std::string input_file, std::string output_file, int max_length_dead_end, int max_length_bubble){
    std::ifstream input(input_file);
    if (!input.is_open())
    {
        std::cout << "Could not open file iicy " << input_file << std::endl;
        exit(1);
    }

    std::string line;
    set<string> bad_contigs;
    unordered_map<string, pair<vector<pair<string, bool>>, vector<pair<string,bool>>>> links_of_contigs;
    unordered_map<string,int> length_of_contigs;

    std::ofstream output(output_file);

    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            links_of_contigs[name] = make_pair(vector<pair<string, bool>>(), vector<pair<string, bool>>());
            length_of_contigs[name] = sequence.size();
        }
        else if (line[0] == 'L'){
            string name1;
            string name2;
            string orientation1, orientation2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2;

            if (orientation1 == "+"){
                links_of_contigs[name1].second.push_back(make_pair(name2, orientation2 == "-"));
            }
            else{
                links_of_contigs[name1].first.push_back(make_pair(name2, orientation2 == "-"));
            }
            if (orientation2 == "+"){
                links_of_contigs[name2].first.push_back(make_pair(name1, orientation1 == "+"));
            }
            else{
                links_of_contigs[name2].second.push_back(make_pair(name1, orientation1 == "+"));
            }
        }
    }

    input.close();
    input.open(input_file);

    for (auto contig : links_of_contigs){
        //first check it it is a small dead end
        if ((contig.second.first.size() == 0 || contig.second.second.size() == 0) && length_of_contigs[contig.first] < max_length_dead_end){
            bad_contigs.insert(contig.first);
        }

        //now check if there is a bubble at either of its ends
        //first left end of the contig
        for (pair<string,bool> neighbor1 : contig.second.first){
            for (pair<string,bool> neighbor2 : contig.second.first){
                if (neighbor1.first != neighbor2.first && bad_contigs.find(neighbor1.first) == bad_contigs.end() && bad_contigs.find(neighbor2.first) == bad_contigs.end()){
                    if (length_of_contigs[neighbor1.first] < max_length_bubble && length_of_contigs[neighbor2.first] < max_length_bubble){
                        //now check that it is a bubble
                        if (links_of_contigs[neighbor1.first].first.size() == 1 && links_of_contigs[neighbor2.first].first.size() == 1
                            && links_of_contigs[neighbor1.first].second.size() == 1 && links_of_contigs[neighbor2.first].second.size() == 1){
                                
                            if (neighbor1.second == neighbor2.second 
                                && links_of_contigs[neighbor1.first].first[0].first == links_of_contigs[neighbor2.first].first[0].first
                                && links_of_contigs[neighbor1.first].second[0].first == links_of_contigs[neighbor2.first].second[0].first){
                                    bad_contigs.insert(neighbor1.first);
                            }
                            else if (neighbor1.second != neighbor2.second 
                                && links_of_contigs[neighbor1.first].first[0].first == links_of_contigs[neighbor2.first].second[0].first
                                && links_of_contigs[neighbor1.first].second[0].first == links_of_contigs[neighbor2.first].first[0].first){
                                    bad_contigs.insert(neighbor1.first);
                            }
                        }
                    }
                }
            }
        }

    }


    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            if (sequence.size() > max_length_dead_end || bad_contigs.find(name) == bad_contigs.end()){
                output << line << "\n";
            }
        }
        else if (line[0] == 'L'){
            string name1;
            string name2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> dont_care >> name2;
            // std::cerr << name1 << " " << name2 << "\n";
            // for (auto i: good_contigs){
            //     std::cerr << i << ",";
            // }
            // std::cerr << "\n";
            if (bad_contigs.find(name1) == bad_contigs.end() && bad_contigs.find(name2) == bad_contigs.end()){
                output << line << "\n";
            }
        }
        else{
            output << line << "\n";
        }
    }
}

/**
 * @brief 
 * 
 * @param input_assembly 
 * @param input_reads 
 * @param gaf_file 
 * @param path_minigraph 
 * @param num_threads 
 */
void realign_reads_on_assembly(vector<SolidBridge>& solid_bridges, vector<SolidPier>& solid_piers, std::string& input_assembly, std::string& input_reads, std::string& gaf_file, std::string& path_minigraph, int num_threads, std::string tmp_folder){
    
    //index the positions of the reads in the read file
    robin_hood::unordered_map<std::string, long int> read_positions;
    std::ifstream read_stream(input_reads);
    string line;
    long int position = 0;
    auto pos = read_stream.tellg();
    while (std::getline(read_stream, line)){
        pos = read_stream.tellg();
        if (line[0] == '>' || line[0] == '@'){
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            read_positions[token.substr(1)] = pos;
        }
    }

    //blacklist all the reads that were used to build the bridges and the piers
    std::set<string> blacklist;
    for (auto solid_bridge : solid_bridges){
        for (auto read_name : solid_bridge.read_names){
            blacklist.insert(read_name);
        }
    }
    for (auto end_contig : solid_piers){
        for (auto read_name : end_contig.read_names){
            blacklist.insert(read_name);
        }
    }

    string tmp_unaligned_reads = tmp_folder +"tmp_unaligned_reads.fasta";
    //go through the gaf file: output the well-aligned alignments in tmp_gaf_file and all the other in tmp_unaligned_reads
    std::ifstream gaf_stream(gaf_file);
    std::ofstream unaligned_reads(tmp_unaligned_reads);
    std::set<string> set_unaligned_reads;

    while (std::getline(gaf_stream, line)){
        string name_of_read;
        int length_of_read;
        int start_of_mapping;
        int end_of_mapping;
        string strand;
        string path;
        int path_length;
        int path_start;
        int path_end;
        string nothing;
        int quality;

        std::istringstream iss(line);
        iss >> name_of_read >> length_of_read >> start_of_mapping >> end_of_mapping >> strand >> path >> path_length >> path_start >> path_end >> nothing >> nothing >> quality;

        int min_length_for_breakpoint = min(0.2*length_of_read, 100.0); //minimum length of overhang to consider a breakpoint
        if ((start_of_mapping > min_length_for_breakpoint || length_of_read-end_of_mapping > min_length_for_breakpoint) && quality == 60 && set_unaligned_reads.find(name_of_read) == set_unaligned_reads.end() 
            && blacklist.find(name_of_read) == blacklist.end()) //read is not mapped end to end ! (or is blacklisted because it should be aligned end-to-end -> to ensure the program finishes)
        {        
            set_unaligned_reads.insert(name_of_read);
            unaligned_reads << ">" << name_of_read << endl;
            std::ifstream read_stream(input_reads);
            read_stream.seekg(read_positions[name_of_read]);
            std::getline(read_stream, line); 
            unaligned_reads << line << endl;
            read_stream.close();
        }
    }
    gaf_stream.close();
    unaligned_reads.close();

    //now run minigraph on the assembly and the unaligned reads
    auto minigraph_run = system((path_minigraph + " -c --secondary=no -t " + std::to_string(num_threads) + " " + input_assembly + " " + tmp_unaligned_reads + " >" + gaf_file + " 2> "+ tmp_folder +"logminigraph.gt.txt").c_str());
    if (minigraph_run != 0){
        std::cerr << "Error running minigraph" << std::endl;
        std::cerr << "Command: " << path_minigraph + " -c --secondary=no -t " + std::to_string(num_threads) + " " + input_assembly + " " + tmp_unaligned_reads + " >" + gaf_file + " 2> "+ tmp_folder +"logminigraph.gt.txt" << std::endl;
        exit(1);
    }

    //remove the temporary files
    remove(tmp_unaligned_reads.c_str());
}

/**
 * @brief Align the reads on the new graph and deduce the coverage of all contigs
 * 
 * @param assembly 
 * @param output_assembly 
 * @param reads 
 * @param num_threads 
 */
void last_cleanup(std::string& input_assembly, std::string& output_assembly, std::string& reads, std::string& gaf_file, std::string& path_minigraph, int num_threads, string tmp_folder){
    

    auto minimap_run = system((path_minigraph + " --secondary=no -t " + std::to_string(num_threads) + " " + input_assembly + " " + reads + " >" + gaf_file + " 2> "+ tmp_folder +"logminigraph.gt.txt").c_str());
    if (minimap_run != 0){
        std::cerr << "Error running minimap: GGY" << std::endl;
        std::cerr << "Command: " << path_minigraph + " --secondary=no -t " + std::to_string(num_threads) + " " + input_assembly + " " + reads + " >" + gaf_file + " 2> "+ tmp_folder +"logminigraph.gt.txt" << std::endl;
        exit(1);
    }

    //index the name and length of the contigs
    robin_hood::unordered_map<std::string, int> length_of_contigs;
    std::ifstream assembly_stream(input_assembly);
    string line;
    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            std::string name_of_contig;
            iss >> token;
            iss >> name_of_contig;
            string sequence;
            iss >> sequence;

            length_of_contigs[name_of_contig] = sequence.size();
        }
    }

    //now go through the gaf file and compute the coverage of all the contigs
    robin_hood::unordered_map<std::string, double> coverage_of_contigs;
    robin_hood::unordered_map<std::string, bool> well_alignness_of_contigs;
    std::ifstream paf_stream(gaf_file);
    while (std::getline(paf_stream, line)){
        string name_of_read;
        int length_of_read;
        int start_of_mapping;
        int end_of_mapping;
        string strand;
        string path;
        int path_length;
        int path_start;
        int path_end;
        string nothing;
        int quality;

        std::istringstream iss(line);
        iss >> name_of_read >> length_of_read >> start_of_mapping >> end_of_mapping >> strand >> path >> path_length >> path_start >> path_end >> nothing >> nothing >> quality;

        if (quality == 60){

            //split the path on > and <
            vector<string> contigs;
            string current_contig;
            for (auto c : path){
                if (c == '<' || c == '>'){
                    if (current_contig.size() > 0){
                        contigs.push_back(current_contig);
                        current_contig = "";
                    }
                }
                else{
                    current_contig += c;
                }
            }
            contigs.push_back(current_contig);

            if (contigs.size() > 1){
                //now add the coverage of the read to the coverage of the contigs
                int total_length = 0;
                //begin by the first contig
                if (coverage_of_contigs.find(contigs[0]) == coverage_of_contigs.end()){
                    coverage_of_contigs[contigs[0]] = 0;
                }
                coverage_of_contigs[contigs[0]] += 1-(path_start/(double)length_of_contigs[contigs[0]]);
                total_length += length_of_contigs[contigs[0]];

                //now the middle contigs
                for (int c = 1 ; c < contigs.size()-1 ; c++){
                    if (coverage_of_contigs.find(contigs[c]) == coverage_of_contigs.end()){
                        coverage_of_contigs[contigs[c]] = 0;
                    }
                    coverage_of_contigs[contigs[c]] += 1;
                    total_length += length_of_contigs[contigs[c]];
                }

                //now the last contig
            
                if (coverage_of_contigs.find(contigs[contigs.size()-1]) == coverage_of_contigs.end()){
                    coverage_of_contigs[contigs[contigs.size()-1]] = 0;
                }
                coverage_of_contigs[contigs[contigs.size()-1]] += (path_end-total_length)/(double)length_of_contigs[contigs[contigs.size()-1]];
            }
            else{
                if (coverage_of_contigs.find(contigs[0]) == coverage_of_contigs.end()){
                    coverage_of_contigs[contigs[0]] = 0;
                }
                coverage_of_contigs[contigs[0]] += (path_end-path_start)/(double)length_of_contigs[contigs[0]];
            }

        }
    }

    //now output the assembly with dp and ln tags, leaving out the contigs with a coverage of 0
    std::ofstream output_stream(output_assembly);
    assembly_stream = std::ifstream (input_assembly);
    std::set<string> contigs_with_enough_coverage;

    while (std::getline(assembly_stream, line)){
        if (line[0] == 'S'){
            std::istringstream iss(line);
            std::string token;
            std::string name_of_contig;
            iss >> token;
            iss >> name_of_contig;
            string sequence;
            iss >> sequence;

            if (coverage_of_contigs.find(name_of_contig) != coverage_of_contigs.end() && coverage_of_contigs[name_of_contig] > 1){
                output_stream << "S\t" << name_of_contig << "\t" << sequence << "\tDP:f:" << coverage_of_contigs[name_of_contig] << "\tLN:i:" << length_of_contigs[name_of_contig] << endl;
                contigs_with_enough_coverage.insert(name_of_contig);
            }
        }
        else if (line[0] == 'L'){
            std::istringstream iss(line);
            std::string token;
            std::string name1, name2, orientation1, orientation2, cigar;
            iss >> token >> name1 >> orientation1 >> name2 >> orientation2 >> cigar;
            if (contigs_with_enough_coverage.find(name1) != contigs_with_enough_coverage.end() && contigs_with_enough_coverage.find(name2) != contigs_with_enough_coverage.end()){
                output_stream << line << endl;
            }
        }
        else{
            output_stream << line << endl;
        }
    }
    output_stream.close();
    assembly_stream.close();
}


/**
 * @brief Output the solid piers and the solid bridges into a file
 * 
 * @param error_file 
 * @param solid_bridges 
 * @param solid_piers 
 */
void output_errors(string error_file, std::vector<SolidBridge>& solid_bridges, std::vector<SolidPier>& solid_piers){
    
    std::ofstream output_stream(error_file);
    output_stream << "Total number of detected structural misassemblies: " << solid_bridges.size() << std::endl;
    output_stream << "Total number of additional detected breakpoints: " << solid_piers.size() << std::endl;
    output_stream << "Misassemblies:" << std::endl;
    for (auto solid_bridge : solid_bridges){
        output_stream << "\tReads support a connection between these two contigs:\n\t\t-" << solid_bridge.contig1 << " " << solid_bridge.position1 << " " << solid_bridge.strand1 << "\n\t\t-" << solid_bridge.contig2 << " " << solid_bridge.position2 << " " << solid_bridge.strand2 << std::endl;
        output_stream << "\tReads supporting the connection:\n";
        for (auto read_name : solid_bridge.read_names){
            output_stream << "\t\t-" << read_name << std::endl;
        }
    }
    output_stream << "\nAdditionnal breakpoints:" << std::endl;
    for (auto solid_pier : solid_piers){
        output_stream << "\tBreakpoint: " << solid_pier.contig << " " << solid_pier.position << " " << solid_pier.strand << std::endl;
        output_stream << "\tReads supporting the breakpoint:\n";
        for (auto read_name : solid_pier.read_names){
            output_stream << "\t\t-" << read_name << std::endl;
        }
    }
    
    output_stream.close();
}

int main(int argc, char *argv[])
{
    //just check which reads are well aligned in "tmp_last_cleanup.gaf"
    // int num_well_aligned_readsz = 0;
    // int num_partially_aligned_readsz = 0;
    // int num_unaligned_readsz = 0;
    // string gagaf = "tmp_last_cleanup.gaf";
    // gagaf = "reads_aligned_on_assembly.gaf";
    // string reads = "hs/tmp/reads.fasta";
    // check_which_reads_are_well_aligned(gagaf, reads, num_well_aligned_readsz, num_partially_aligned_readsz, num_unaligned_readsz);
    // cout << "ECICIC " << num_well_aligned_readsz << " " << num_partially_aligned_readsz << " " << num_unaligned_readsz << endl;
    // exit(1);

    //build a clipp.h command line parser
    bool help = false;
    bool print_version = false;
    bool correct = false;
    string input_assembly, input_reads, gaf_file, output_scaffold, error_file, correct_string;
    int num_threads = 1;
    int min_num_reads_for_link = 5;
    string path_minigraph = "minigraph";
    string path_minimap2 = "minimap2";
    string path_racon = "racon";
    string path_to_raven = "raven";
    string path_tmp_folder = "./";
    string path_new_reads = "";

    auto cli = clipp::group(
        clipp::required("-i", "--input_assembly").doc("input assembly in gfa format") & clipp::value("input_assembly", input_assembly),
        clipp::required("-r", "--input_reads").doc("input reads in fasta/q format") & clipp::value("input_reads", input_reads),
        clipp::required("-m", "--mode").doc("mode: correct or detect") & clipp::value("mode", correct_string),
        clipp::required("-e", "--output_errors").doc("output file describing the errors found in the assembly") & clipp::value("output_errors", error_file),
        clipp::option("-o", "--output_assembly").doc("output assembly in gfa format (required if correct mode)") & clipp::value("output_assembly", output_scaffold),
        clipp::option("-d", "--output_non_duplexed_reads").doc("file to output a file of non-duplexed reads") & clipp::value("output_non_duplexed_reads", path_new_reads),
        clipp::option("-p", "--path-to-tmp-folder").doc("path to a temporary folder where the intermediate files will be stored [./]") & clipp::value("path-to-tmp-folder", path_tmp_folder),
        clipp::option("-b", "--minimum-number-of-reads").doc("minimum number of reads to support a breakpoint [5]") & clipp::value("minimum-number-of-reads", min_num_reads_for_link),
        clipp::option("-t", "--threads").doc("number of threads to use for minigraph [1]") & clipp::value("threads", num_threads),
        clipp::option("-g", "--gaf_file").doc("gaf file (NO SECONDARY ALIGNMENTS). Will be generated with minigraph if not provided") & clipp::value("gaf_file", gaf_file),
        clipp::option("--minigraph").doc("path to minigraph") & clipp::value("minigraph", path_minigraph),
        clipp::option("--minimap2").doc("path to minimap2") & clipp::value("minimap2", path_minimap2),
        clipp::option("--racon").doc("path to racon") & clipp::value("racon", path_racon),
        clipp::option("--path-to-raven").doc("path to raven") & clipp::value("path-to-raven", path_to_raven),
        clipp::option("-h", "--help").set(help).doc("print this help message and exit"),
        clipp::option("-v", "--version").set(print_version).doc("print version information and exit")
    );

    //parse the command line
    //if the command line is invalid or the user asked for help, print the usage and exit
    if(!parse(argc, argv, cli) || help) {
        cout << make_man_page(cli, argv[0]);
        return 0;
    }

    //add a slash at the end of the path to the temporary folder if it is not there
    if (path_tmp_folder.back() != '/'){
        path_tmp_folder += "/";
    }
    //check that the temporary folder exists, if not create it
    system(("mkdir -p " + path_tmp_folder).c_str());

    if (correct_string == "correct"){
        correct = true;
    }
    else if (correct_string == "detect"){
        correct = false;
    }
    else{
        std::cerr << "Error: mode should be either 'correct' or 'detect'" << std::endl;
        return 1;
    }

    if (print_version){
        cout << "version: " << version << endl;
        cout << "last update: " << last_update << endl;
        cout << "author: Roland Faure" << endl;
        return 0;
    }

    //welcome message
    cout << endl;
    cout << "   ****************************************" << endl;
    cout << "   *                                      *" << endl;
    cout << "   *                                      *" << endl;
    cout << "   *            Welcome to the            *" << endl;
    cout << "   *             Genome Tailor            *" << endl;
    cout << "   *   Tailor-made assemblies since 2024  *" << endl;
    cout << "   *                                      *" << endl;
    cout << "   *                                      *" << endl;
    cout << "   ****************************************" << endl;

    cout << endl;
    cout << "Running Genome Tailor version " << version << endl;
    cout << "Last update: " << last_update << endl;
    cout << "Author: Roland Faure" << endl;
    cout << "Command line: ";
    for (auto i = 0 ; i < argc ; i++){
        cout << argv[i] << " ";
    }
    cout << endl << endl;

    //check the dependencies
    cout << "Checking dependencies..." << endl;
    bool minimap_ok, minigraph_ok, racon_ok, raven_ok;
    minimap_ok = (system((path_minimap2 + " -h >trash.tmp 2>trash.tmp").c_str()) == 0);
    minigraph_ok = (system((path_minigraph + " --version >trash.tmp 2>trash.tmp").c_str()) == 0);
    racon_ok = (system((path_racon + " -h >trash.tmp 2>trash.tmp").c_str()) == 0);
    raven_ok = (path_to_raven != "" && system((path_to_raven + " --version >trash.tmp 2>trash.tmp").c_str()) == 0);
    system("rm trash.tmp");

    // Print the table of dependencies
    std::cout << "_______________________________" << std::endl;
    std::cout << "|    Dependency     |  Found  |" << std::endl;
    std::cout << "|-------------------|---------|" << std::endl;
    std::cout << "|    minigraph      |   " << (minigraph_ok ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    std::cout << "|    minimap2       |   " << (minimap_ok ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    std::cout << "|    racon          |   " << (racon_ok ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    std::cout << "|    raven          |   " << (raven_ok ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    std::cout << "-------------------------------" << std::endl;

    if (!minimap_ok || !minigraph_ok || !racon_ok || !raven_ok){
        std::cout << "Error: some dependencies are missing. Please install them or provide a valid path with the options." << std::endl;
        return 1;
    }


    std::string input_assembly_format = input_assembly.substr(input_assembly.find_last_of(".") + 1);
    std::string input_reads_format = input_reads.substr(input_reads.find_last_of(".") + 1);
    std::string output_scaffold_format = output_scaffold.substr(output_scaffold.find_last_of(".") + 1);

    bool input_assembly_is_gfa = (input_assembly_format == "gfa");
    bool input_reads_is_fasta = (input_reads_format == "fasta" || input_reads_format == "fa");
    bool input_reads_is_fastq = (input_reads_format == "fastq" || input_reads_format == "fq");
    bool input_assembly_exists = (std::ifstream(input_assembly).good());
    bool input_reads_exists = (std::ifstream(input_reads).good()); 
    bool output_scaffold_is_gfa = (output_scaffold_format == "gfa" || !correct);
    bool gaf_file_is_gaf = (gaf_file.substr(gaf_file.find_last_of(".") + 1) == "gaf");


    cout << "Checking file formats..." << endl;
    //print a table of the input and output files
    std::cout << "_________________________________" << std::endl;
    std::cout << "|       File        | Format ok |" << std::endl;
    std::cout << "|-------------------|-----------|" << std::endl;
    if (input_assembly_exists){
        std::cout << "|  input_assembly   |  " << (input_assembly_is_gfa ? GREEN_TEXT "  gfa  " : RED_TEXT "not gfa") << RESET_TEXT "  |" << std::endl;
    }
    else{
        std::cout << "|  input_assembly   | " << RED_TEXT "not found" << RESET_TEXT " |" << std::endl;
    }
    if (input_reads_exists){
        std::cout << "|  input_reads      |" << (input_reads_is_fasta || input_reads_is_fastq ? GREEN_TEXT "  fasta/q  " : RED_TEXT "not fasta/q") << RESET_TEXT "|" << std::endl;
    }
    else{
        std::cout << "|  input_reads      | " << RED_TEXT "not found" << RESET_TEXT " |" << std::endl;
    }
    if (correct){
        std::cout << "|  output_scaffold  |  " << (output_scaffold_is_gfa ? GREEN_TEXT "  gfa  " : RED_TEXT "not gfa") << RESET_TEXT "  |" << std::endl;
    }
    if (gaf_file != ""){
        std::cout << "|  gaf_file         |  " << (gaf_file_is_gaf ? GREEN_TEXT "  gaf  " : RED_TEXT "not gaf") << RESET_TEXT "  |" << std::endl;
    }
    std::cout << "---------------------------------" << std::endl << std::endl;

    if (!input_assembly_is_gfa || (!input_reads_is_fasta && !input_reads_is_fastq) || !output_scaffold_is_gfa || (gaf_file != "" && !gaf_file_is_gaf)){
        std::cout << "Error: some input or output files are not in the right format. Please check the format and try again." << std::endl;
        return 1;
    }
    if (!input_assembly_exists || !input_reads_exists){
        std::cout << "Error: some input files are missing. Please check the file paths and try again." << std::endl;
        return 1;
    }
    if (correct && output_scaffold == ""){
        std::cout << "Error: the output assembly file (-o) is required in correct mode. Please provide a valid path and try again." << std::endl;
        return 1;
    }

    //if the gaf file is not provided, create it with minigraph
    if (gaf_file == ""){
        gaf_file = path_tmp_folder + "reads_aligned_on_assembly.gaf";
        cout << "0) Creating a gaf file by aligning reads on assembly with minigraph, creating the file " << gaf_file << endl;
        auto minigraph_run = system((path_minigraph + " -c --secondary=no -t " + std::to_string(num_threads) + " " + input_assembly + " " + input_reads + " >" + gaf_file + " 2> "+path_tmp_folder+"logminigraph.gt.txt").c_str());
        if (minigraph_run != 0){
            cout << "Error: minigraph failed. Please check the log file logminigraph.gt.txt and try again. The command line tried was:" << endl;
            cout << path_minigraph + " -c --secondary=no -t " + std::to_string(num_threads) + " " + input_assembly + " " + input_reads + " >" + gaf_file + " 2> "+path_tmp_folder+"logminigraph.gt.txt" << endl;
            return 1;
        }
        // cout << "NOT RUNNING MINIGRPAH, TO BE REMOVED" << endl;
        //delete the log file
        system("rm logminigraph.gt.txt");
    }    

    //count the number of reads that align end-to-end, partially and not at all
    int num_well_aligned_reads_start = 0;
    int num_partially_aligned_reads_start = 0;
    int num_unaligned_reads_start = 0;
    int num_jumpings_reads_start = 0;
    int num_undetected_duplex_reads_start = 0;
    count_unaligned_reads(gaf_file, input_reads, num_unaligned_reads_start);
    int num_well_aligned_reads_end = 0;
    int num_partially_aligned_reads_end = 0;
    int num_unaligned_reads_end = 0;
    int num_jumpings_reads_end = 0;
    int num_undetected_duplex_reads_end = 0;

    //reassemble unaligned reads with raven
    string assembly_completed = path_tmp_folder+"tmp_assembly_completed.gfa";
    string gaf_completed = path_tmp_folder+"tmp_gaf_completed.gaf";
    if (correct){
        cout << "Reassembling the unaligned reads with raven..." << endl;
        reassemble_unaligned_reads(gaf_file, input_reads, input_assembly, assembly_completed, gaf_completed, path_to_raven, path_minigraph, num_threads, path_tmp_folder);
    }
    else{
        system(("cp " + input_assembly + " " + assembly_completed).c_str());
        system(("cp " + gaf_file + " " + gaf_completed).c_str());
    }
    // system(("cp " + input_assembly + " " + assembly_completed).c_str());
    // system(("cp " + gaf_file + " " + gaf_completed).c_str());
    // cout << "NOT RUNNING RAVEN, TO BE REMOVED" << endl;

    if (correct){
        cout << "\n==== Now looping and iteratively modify the GFA until all reads align end-to-end on the assembly graph ====" << endl;
        bool some_reads_still_unaligned = true;
        int iteration = 0;
        while (some_reads_still_unaligned){

            cout << endl << " Loop iteration " << iteration << "... " << endl;

            //inventoriate the bridges
            cout << "  - Going through the gaf file and listing the reads that do not align end-to-end on the assembly graph..." << endl;
            std::vector<Bridge> bridges;
            std::vector<Pier> piers;
            std::set<string> list_of_duplex_reads;
            if (iteration == 0){
                inventoriate_bridges_and_piers(gaf_completed, bridges, piers, assembly_completed, num_well_aligned_reads_start, num_partially_aligned_reads_start, num_jumpings_reads_start, num_undetected_duplex_reads_start, list_of_duplex_reads);
                //discard all the duplex reads both in a new rreads file if given as input
                if (path_new_reads != ""){
                    cout << "Outputting non-duplexed reads in " << path_new_reads << ", discarding in total " << list_of_duplex_reads.size() << " reads" << endl;
                    std::ofstream new_reads_file(path_new_reads);
                    std::ifstream read_stream(input_reads);
                    string line;
                    char format = 'f';
                    if (input_reads.substr(input_reads.find_last_of(".") + 1) == "fastq" || input_reads.substr(input_reads.find_last_of(".") + 1) == "fq"){
                        format = 'q';
                    }
                    int line_number = 0;
                    while (std::getline(read_stream, line)){
                        if ((format == 'f' && line_number % 2 == 0) || (format == 'q' && line_number % 4 == 0)){
                            std::istringstream iss(line);
                            std::string token;
                            iss >> token;
                            if (list_of_duplex_reads.find(token.substr(1)) == list_of_duplex_reads.end()){
                                new_reads_file << ">"<< line.substr(1) << std::endl;
                                std::getline(read_stream, line);
                                new_reads_file << line << std::endl;
                            }
                        }
                        line_number++;
                    }
                    new_reads_file.close();
                    read_stream.close();
                }
            }
            else{
                int nothing;
                inventoriate_bridges_and_piers(gaf_completed, bridges, piers, assembly_completed, num_well_aligned_reads_end, num_partially_aligned_reads_end, num_jumpings_reads_end, num_undetected_duplex_reads_end, list_of_duplex_reads);
            }

            // cout << "Found " << bridges.size() << " bridges and " << piers.size() << " piers." << endl;

            //agregate the bridges
            cout << "  - Pooling the reads that display similar behaviour on the assembly graph..." << endl;
            std::vector<SolidBridge> solid_bridges;
            std::vector<SolidPier> solid_piers;
            agregate_bridges_and_piers(bridges, piers, solid_bridges, solid_piers, 1000, min_num_reads_for_link);
            if (iteration == 0){
                output_errors(error_file, solid_bridges, solid_piers);
            }


            if ((/*solid_piers.size() == 0 &&*/ solid_bridges.size() == 0)){
                some_reads_still_unaligned = false;
                cout << "\n==== Graph cannot be improved further. Moving on to the last step, computing the coverage ====\n" << endl;
                break;
            }
            else{
                cout << "  - Found " << solid_bridges.size() << " solid bridges and " << solid_piers.size() << " solid piers." << endl;
            }

            //transform the solid bridges into links in the assembly
            cout << "  - Computing the exact location of new links in the GFA and gap-filling if necessary..." << endl;
            std::vector<Link> links;
            std::vector<End_contig> end_contigs;
            transform_bridges_in_links(solid_bridges, input_reads, assembly_completed, links, path_minimap2, path_racon, path_tmp_folder);
            //piers are tricky: you just need one read to align on one contig: this is a very weak signal
            // build_piers(solid_piers, input_reads, assembly_completed, end_contigs, path_minimap2, path_racon);

            //create a gfa file from the assembly and the links
            cout << "  - Outputting the new assembly..." << endl;
            string tmp_gfa = path_tmp_folder+"tmp_drt_gfa.gfa";
            create_gfa(assembly_completed, tmp_gfa, links, end_contigs);

            //shave the graph of small dead ends, that probably result from polishing errors
            cout << "  - Shaving the graph of small dead ends and popping small bubbles that may have appeared in previous steps..." << endl;
            shave_and_pop(tmp_gfa, output_scaffold, 60, 20);
            system(("cp " + output_scaffold + " " + assembly_completed).c_str()); //so that the next iteration starts from the new assembly

            string gfa_iteration = path_tmp_folder+"tmp_gfa_iteration_"+ std::to_string(iteration) +"_.gfa";
            system(("cp " + output_scaffold + " " + gfa_iteration).c_str());

            //realign the reads on the new assembly
            cout << "  - Realigning the reads on the new assembly..." << endl;
            realign_reads_on_assembly(solid_bridges, solid_piers, assembly_completed, input_reads, gaf_completed, path_minigraph, num_threads, path_tmp_folder); //the new result is stored in gaf_completed

            //remove the temporary gfa file
            system(("rm " + tmp_gfa).c_str());

            iteration++;
        }

        //remove the temporary files
        // system(("rm " + assembly_completed).c_str());
        // system(("rm " + gaf_completed).c_str());

        cout << "Re-aligning the reads on the final assembly to compute coverage and remove all non-covered contigs..." << endl;
        int num_well_aligned_reads = 0;
        int num_partially_aligned_reads = 0;
        int num_unaligned_reads = 0;
        string final_gaf = path_tmp_folder+"tmp_last_cleanup.gaf";
        last_cleanup(assembly_completed, output_scaffold, input_reads, final_gaf, path_minigraph, num_threads, path_tmp_folder);

        //count the number of reads that align end-to-end, partially and not at all
        count_unaligned_reads(final_gaf, input_reads, num_unaligned_reads);

        //print the number of reads that aligned end-to-end, partially and not at all before and after 
        cout << "____________________________________________________________________________________________________" << endl;
        cout << "|                                           |    Before GenomeTailor    |    After GenomeTailor    |" << endl;
        cout << "|-------------------------------------------|---------------------------|--------------------------|" << endl;
        string num_well_aligned_reads_start_str = std::to_string(num_well_aligned_reads_start);
        cout << "|    Number of end-to-end aligned reads     | " << num_well_aligned_reads_start_str;
        for (int i = 0 ; i < 25-num_well_aligned_reads_start_str.size() ; i++){
            cout << " ";
        }
        cout << " | " << num_well_aligned_reads;
        for (int i = 0 ; i < 25-std::to_string(num_well_aligned_reads).size() ; i++){
            cout << " ";
        }
        cout << "|" << endl;

        string num_undetected_duplex_reads_end_str = std::to_string(num_undetected_duplex_reads_end);
        cout << "|    Number of undetected duplex reads      | " << num_undetected_duplex_reads_end_str;
        for (int i = 0 ; i < 25-num_undetected_duplex_reads_end_str.size() ; i++){
            cout << " ";
        }
        cout << " | " << num_undetected_duplex_reads_start;
        for (int i = 0 ; i < 25-std::to_string(num_undetected_duplex_reads_start).size() ; i++){
            cout << " ";
        }
        cout << "|" << endl;

        string num_jumpings_reads_end_str = std::to_string(num_jumpings_reads_end);
        cout << "|    Number of jumping reads                 | " << num_jumpings_reads_end_str;
        for (int i = 0 ; i < 25-num_jumpings_reads_end_str.size() ; i++){
            cout << " ";
        }
        cout << " | " << num_jumpings_reads_start;
        for (int i = 0 ; i < 25-std::to_string(num_jumpings_reads_start).size() ; i++){
            cout << " ";
        }
        cout << "|" << endl;

        string num_partially_aligned_reads_start_str = std::to_string(num_partially_aligned_reads_start);
        cout << "|    Number of partially aligned reads       | " << num_partially_aligned_reads_start_str;
        for (int i = 0 ; i < 25-num_partially_aligned_reads_start_str.size() ; i++){
            cout << " ";
        }
        cout << " | " << num_partially_aligned_reads;
        for (int i = 0 ; i < 25-std::to_string(num_partially_aligned_reads).size() ; i++){
            cout << " ";
        }
        cout << "|" << endl;

        string num_unaligned_reads_start_str = std::to_string(num_unaligned_reads_start);
        cout << "|    Number of unaligned reads              | " << num_unaligned_reads_start_str;
        for (int i = 0 ; i < 25-num_unaligned_reads_start_str.size() ; i++){
            cout << " ";
        }
        cout << " | " << num_unaligned_reads;
        for (int i = 0 ; i < 25-std::to_string(num_unaligned_reads).size() ; i++){
            cout << " ";
        }
        cout << "|" << endl;
        cout << "----------------------------------------------------------------------------------------------------" << endl;
    }
    else{ //just detect the misassemblies
        //inventoriate the bridges
        cout << "  - Going through the gaf file and listing the reads that do not align end-to-end on the assembly graph..." << endl;
        std::vector<Bridge> bridges;
        std::vector<Pier> piers;
        int num_well_aligned_reads = 0;
        int num_partially_aligned_reads = 0;
        int num_unaligned_reads = 0;
        int num_jumpings_reads = 0;
        int num_undetected_duplex_reads = 0;
        std::set<string> list_of_duplex_reads;
        inventoriate_bridges_and_piers(gaf_completed, bridges, piers, assembly_completed, num_well_aligned_reads, num_partially_aligned_reads, num_jumpings_reads, num_undetected_duplex_reads, list_of_duplex_reads);
        
        //discard all the duplex reads both in a new rreads file if given as input
        if (path_new_reads != ""){
            cout << "Outputting non-duplexed reads in " << path_new_reads << ", discarding in total " << list_of_duplex_reads.size() << " reads" << endl;
            std::ofstream new_reads_file(path_new_reads);
            std::ifstream read_stream(input_reads);
            string line;
            char format = 'f';
            if (input_reads.substr(input_reads.find_last_of(".") + 1) == "fastq" || input_reads.substr(input_reads.find_last_of(".") + 1) == "fq"){
                format = 'q';
            }
            int line_number = 0;
            while (std::getline(read_stream, line)){
                if ((format == 'f' && line_number % 2 == 0) || (format == 'q' && line_number % 4 == 0)){
                    std::istringstream iss(line);
                    std::string token;
                    iss >> token;
                    if (list_of_duplex_reads.find(token.substr(1)) == list_of_duplex_reads.end()){
                        new_reads_file << ">"<< line.substr(1) << std::endl;
                        std::getline(read_stream, line);
                        new_reads_file << line << std::endl;
                    }
                }
                line_number++;
            }
            new_reads_file.close();
            read_stream.close();
        }
        
        count_unaligned_reads(gaf_completed, input_reads, num_unaligned_reads);
        cout << "Found " << bridges.size() << " bridges and " << piers.size() << " piers." << endl;

        //agregate the bridges
        cout << "  - Pooling the reads that display similar behaviour on the assembly graph..." << endl;
        std::vector<SolidBridge> solid_bridges;
        std::vector<SolidPier> solid_piers;
        agregate_bridges_and_piers(bridges, piers, solid_bridges, solid_piers, 1000, min_num_reads_for_link);

        //output the errors
        output_errors(error_file , solid_bridges, solid_piers);

        //output the number of reads that align end-to-end, partially and not at all
        cout << "_________________________________________________________________________" << endl;
        cout << "|    Number of end-to-end aligned reads     | " << num_well_aligned_reads << " |" << endl;
        cout << "|    Number of partially aligned reads       | " << num_partially_aligned_reads << " |" << endl;
        cout << "|    Number of unaligned reads              | " << num_unaligned_reads << " |" << endl;
        cout << "|    Number of jumping reads                 | " << num_jumpings_reads << " |" << endl;
        cout << "|    Number of undetected duplex reads      | " << num_undetected_duplex_reads << " |" << endl;
        cout << "-------------------------------------------------------------------------" << endl;
    
    }

    cout << endl << "Done! Customer service at github.com/RolandFaure/GenomeTailor" << endl;
}