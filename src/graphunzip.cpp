#include "graphunzip.h"

#include <iostream>
#include <fstream>
#include <sstream>

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::mutex;
using std::pair;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::stringstream;
using robin_hood::unordered_map;
using std::set;

string reverse_complement(string &seq){
    string rev_comp = "";
    for (int c = seq.size() - 1 ; c >= 0 ; c--){
        if (seq[c] == 'A'){
            rev_comp += 'T';
        }
        else if (seq[c] == 'T'){
            rev_comp += 'A';
        }
        else if (seq[c] == 'C'){
            rev_comp += 'G';
        }
        else if (seq[c] == 'G'){
            rev_comp += 'C';
        }
        else{
            rev_comp += 'N';
        }
    }
    return rev_comp;
}

class Segment{
    public:

        string name;
        int ID;
        vector<pair<vector<pair<int,int>>, vector<string>>>  links; //first vector is for the links to the left, second vector is for the links to the right. Each link is the ID of the neighbor and its end (0 for left, 1 for right) and the CIGAR

        Segment(){};

        Segment(string name, int ID, long int pos_in_file, double coverage){
            this->name = name;
            this->ID = ID;
            this->pos_in_file = pos_in_file;
            this->coverage = coverage;
            this->haploid = false;
            this->links = vector<pair<vector<pair<int,int>>, vector<string>>>(2);
        }

        Segment(string name, int ID, vector<pair<vector<pair<int,int>>, vector<string>>> links, long int pos_in_file, double coverage){
            this->name = name;
            this->ID = ID;
            this->links = links;
            this->pos_in_file = pos_in_file;
            this->coverage = coverage;
            this->haploid = false;
            this->links = vector<pair<vector<pair<int,int>>, vector<string>>>(2);
        }

        bool is_haploid(){return this->haploid;}
        vector<pair<int,bool>> get_consensus_left(){return this->consensus_left;}
        vector<pair<int,bool>> get_consensus_right(){return this->consensus_right;}
        long int get_pos_in_file(){return this->pos_in_file;}
        double get_coverage(){return this->coverage;}

        string get_seq(string& gfa_file){
            
            if (seq != ""){
                // cout << "the seq is already loaded" << endl;
                return seq;
            }

            ifstream gfa(gfa_file);
            gfa.seekg(pos_in_file);
            string line;
            getline(gfa, line);
            stringstream ss(line);
            string nothing, name, seq;
            ss >> nothing >> name >> seq;
            gfa.close();
            return seq;
        }

        void add_neighbor(vector<pair<int,bool>> new_neighbor, bool left);

        void decrease_coverage(double coverage_out){
            coverage -= coverage_out;
            if (coverage < 1){
                coverage = 1;
            }
        }

        void compute_consensuses();

        vector<vector<pair<int,bool>>> get_strong_neighbors_left();
        vector<vector<pair<int,bool>>> get_strong_neighbors_right();

        vector<pair<vector<pair<int,bool>>, int>> get_neighbors_left(){return neighbors_left;}
        vector<pair<vector<pair<int,bool>>, int>> get_neighbors_right(){return neighbors_right;}

        //the hash of the segment is the hash of the name
        size_t hash() const{
            return std::hash<string>{}(name);
        }

        string seq;

    private:
        //consensus sequences of contigs left and right
        vector<pair<int,bool>> consensus_left;
        vector<pair<int,bool>> consensus_right;

        vector<vector<pair<int,bool>>> representative_neighbors_left;
        vector<vector<pair<int,bool>>> representative_neighbors_right;

        vector<pair<vector<pair<int,bool>>, int>> neighbors_left; //each path is a vector of pairs (ID, orientation) of the contigs in the path and a strength (number of reads)
        vector<pair<vector<pair<int,bool>>, int>> neighbors_right; //each path is a vector of pairs (ID, orientation) of the contigs in the path and a strength (number of reads)

        long int pos_in_file;
        double coverage;

        bool haploid;

};

void Segment::add_neighbor(vector<pair<int,bool>> new_neighbor, bool left){


    if (new_neighbor.size() == 0){
        return;
    }
    if (left){
        //only add the neighbor if it is not contained in another neighbor. If another neighobr is contained in the new neighbor, remove it
        bool discard_new_neighbor = false;
        int old_neighbor_to_discard = -1;
        int idx = 0;
        for (pair<vector<pair<int,bool>>,int> neighbor : neighbors_left){
            for (int c = 0 ; c < new_neighbor.size() ; c++){
                if (c >= neighbor.first.size()){
                    old_neighbor_to_discard = idx;
                    break;
                }
                else if (neighbor.first[c] != new_neighbor[c]){ //a difference is found
                    break;
                }
                else if (c == new_neighbor.size()-1){
                    discard_new_neighbor = true;
                    neighbors_left[idx].second += 1;
                    return;
                }
            }
            idx++;
        }

        //remove the old neighbor if needed
        int new_strength = 1; 
        if (old_neighbor_to_discard != -1){
            new_strength += neighbors_left[old_neighbor_to_discard].second;
            neighbors_left.erase(neighbors_left.begin() + old_neighbor_to_discard);
        }
        neighbors_left.push_back({new_neighbor, new_strength});

        // if (neighbors_left.size() > 10){
        //     cout << "*******\nnew neighobr: " << endl;
        //     for (pair<int,bool> contig : new_neighbor){
        //         cout << contig.first << " ";
        //     }
        //     cout << endl;

        //     cout << "existing neighbors: " << endl;
        //     for (vector<pair<int,bool>> neighbor : neighbors_left){
        //         for (pair<int,bool> contig : neighbor){
        //             cout << contig.first << " ";
        //         }
        //         cout << endl;
        //     }
        // }
    }
    else{
        //only add the neighbor if it is not contained in another neighbor. If another neighobr is contained in the new neighbor, remove it
        bool discard_new_neighbor = false;
        int old_neighbor_to_discard = -1;
        int idx = 0;
        for (pair<vector<pair<int,bool>>,int> neighbor : neighbors_right){
            for (int c = 0 ; c < new_neighbor.size() ; c++){
                if (c >= neighbor.first.size()){
                    old_neighbor_to_discard = idx;
                    break;
                }
                else if (neighbor.first[c] != new_neighbor[c]){
                    break;
                }
                else if (c == new_neighbor.size()-1){
                    discard_new_neighbor = true;
                    neighbors_right[idx].second += 1;
                    return;
                }
            }
            idx++;
        }

        //remove the old neighbor if needed
        int new_strength = 1;
        if (old_neighbor_to_discard != -1){
            new_strength += neighbors_right[old_neighbor_to_discard].second;
            neighbors_right.erase(neighbors_right.begin() + old_neighbor_to_discard);
        }
        neighbors_right.push_back({new_neighbor, new_strength});
    }
}


void load_GFA(string gfa_file, vector<Segment> &segments, unordered_map<string, int> &segment_IDs){
    //load the segments from the GFA file
    
    //in a first pass index all the segments by their name
    ifstream gfa(gfa_file);
    string line;
    while (getline(gfa, line)){
        if (line[0] == 'S'){
            long int pos_in_file = gfa.tellg() - line.size() - 1;
            stringstream ss(line);
            string nothing, name;
            ss >> nothing >> name;

            double coverage = 0;
            //try to find a DP: tag
            string tag;
            while (ss >> tag){
                if (tag.substr(0, 3) == "DP:"){
                    coverage = std::stof(tag.substr(5, tag.size() - 5));
                }
            }

            Segment s(name, segments.size(), vector<pair<vector<pair<int,int>>, vector<string>>>(2), pos_in_file, coverage);

            segment_IDs[name] = s.ID;
            segments.push_back(s);
        }
    }
    gfa.close();

    //in a second pass, load the links
    gfa.open(gfa_file);
    while (getline(gfa, line)){
        if (line[0] == 'L'){
            stringstream ss(line);
            string nothing, name1, name2;
            string orientation1, orientation2;
            string cigar;

            ss >> nothing >> name1 >> orientation1 >> name2 >> orientation2 >> cigar;

            int end1 = 1;
            int end2 = 0;

            if (orientation1 == "-"){
                end1 = 0;
            }
            if (orientation2 == "-"){
                end2 = 1;
            }

            int ID1 = segment_IDs[name1];
            int ID2 = segment_IDs[name2];

            //check that the link did not already exist
            bool already_exists = false;
            for (pair<int,int> link : segments[ID1].links[end1].first){
                if (link.first == ID2 && link.second == end2){
                    already_exists = true;
                }
            }
            if (!already_exists){

                segments[ID1].links[end1].first.push_back({ID2, end2});
                segments[ID1].links[end1].second.push_back(cigar);

                segments[ID2].links[end2].first.push_back({ID1, end1});
                segments[ID2].links[end2].second.push_back(cigar);
            }
        }
    }
    gfa.close();
}

void load_GAF(string gaf_file, vector<Segment> &segments, unordered_map<string, int> &segments_IDs){
    //load the paths from the GAF file
    ifstream gaf(gaf_file);
    string line;
    int nb_lines = 0;
    while (getline(gaf, line)){

        //get the name of the read as the first space-delimited word, and the path as the 8th tab-delimited column
        string name_of_read;
        string path;
        string nothing;
        int field = 0;
        string field_content = "";
        stringstream line_stream(line);
        while (std::getline(line_stream, field_content, '\t')){
            if (field == 0){
                name_of_read = field_content.substr(0, field_content.find_first_of(' '));
            }
            else if (field == 5){
                path = field_content;
            }
            field++;
        }

        if (nb_lines % 100000 == 0){
            cout << nb_lines << " lines read in load_GAF\r" << std::flush;
        }
        // if (nb_lines > 1000000){
        //     cout << "NOT READING EVERYTHING FOR DEBGUGGING PUROSPS" << endl;
        //     break;
        // }
        nb_lines += 1;
        
        bool orientation_now;
        string str_now;
        vector<pair<int,bool>> segments_now;
        vector<vector<pair<int,bool>>> paths_in_this_read;

        for (int c = 0 ; c < path.size(); c++){
            if (path[c] == '>' || path[c] == '<' || c == path.size()-1){

                if (c == path.size()-1){
                    str_now += path[c];
                }

                if (c != 0){
                    if (segments_IDs.find(str_now) == segments_IDs.end()){
                        cout << "Error: the segment " << str_now << " is in the gaf but cannot be found in the GFA" << endl;
                        cout << line << endl;
                        exit(1);
                    }
                    int new_contig_ID = segments_IDs[str_now];

                    //check that the segment can indeed be found next to the previous segment
                    bool found = false;
                    if (segments_now.size() > 0){
                        for (pair<int,int> link : segments[segments_now.back().first].links[segments_now.back().second].first){
                            if (link.first == new_contig_ID && link.second == !orientation_now){
                                found = true;
                            }
                        }
                    }
                    if (segments_now.size() == 0 || found){
                        segments_now.push_back({new_contig_ID, orientation_now});
                    }
                    else{ //then the neighbor is not the neighbor of a previous segment, cut
                        paths_in_this_read.push_back(segments_now);
                        segments_now = {{new_contig_ID, orientation_now}};  
                    }
                }

                if (path[c] == '>'){
                    orientation_now = true;
                }
                else{
                    orientation_now = false;
                }
                str_now = "";
            }
            else{
                str_now += path[c];
            }
        }
        paths_in_this_read.push_back(segments_now);

        for (vector<pair<int,bool>> pair_path : paths_in_this_read){

            //compute the reverse path too
            vector<pair<int,bool>> reverse_path_contigs;
            for (int contig = pair_path.size() - 1 ; contig >= 0 ; contig--){
                reverse_path_contigs.push_back({pair_path[contig].first, !pair_path[contig].second});
            }


            for (int contig = 0 ; contig < pair_path.size()  ; contig++){
                if (pair_path[contig].second){
                    //push back the vector beginning at contig, excluding contig
                    if (contig != pair_path.size() - 1){
                        segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(pair_path.begin() + contig + 1, pair_path.end()), false);
                    }
                    segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(reverse_path_contigs.begin() + reverse_path_contigs.size() - contig, reverse_path_contigs.end()), true);
                } 
                else{
                    //push back the vector beginning at contig
                    if (contig != pair_path.size() - 1){
                        segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(pair_path.begin() + contig + 1, pair_path.end()), true);
                    }
                    segments[pair_path[contig].first].add_neighbor(vector<pair<int,bool>>(reverse_path_contigs.begin() + reverse_path_contigs.size() - contig, reverse_path_contigs.end()), false);
                }

            }
        }
    }
    gaf.close();
}

/**
 * @brief Compute the consensual path of contigs left and right
 * 
 */
void Segment::compute_consensuses(){

    this->haploid = true;

    //go through the nighbors left and discard all reads that are subreads of another read //already done in add_neighbor
    // vector<bool> to_discard_left = vector<bool>(neighbors_left.size(), false);
    // for (int p = 0 ; p < neighbors_left.size() ; p++){
    //     for(int p2 = 0 ; p2 < p ; p2++){
    //         if (!to_discard_left[p2]){
    //             bool is_subread = true;
    //             for (int c = 0 ; c < neighbors_left[p].first.size() ; c++){
    //                 if (c >= neighbors_left[p2].first.size()){
    //                     to_discard_left[p2] = true;
    //                     is_subread = false;
    //                     break;
    //                 }
    //                 else if (neighbors_left[p2].first[c] != neighbors_left[p].first[c]){
    //                     is_subread = false;
    //                     break;
    //                 }
    //             }
    //             if (is_subread){
    //                 to_discard_left[p] = true;
    //                 break;
    //             }
    //         }
    //     }
    // }

    //go through the remaining reads: if there are more than one strong then not haploid, else set the only strong read as the consensus
    vector<int> strengths_left;
    for (auto neighbor : neighbors_left){
        strengths_left.push_back(neighbor.second);
    }
    //find the two strongest neighbors
    int max_strength = 0;
    int idx_max = -1;
    int second_max_strength = 0;
    int idx_second_max = -1;
    for (int s = 0 ; s < strengths_left.size() ; s++){
        if (strengths_left[s] > max_strength){
            second_max_strength = max_strength;
            idx_second_max = idx_max;
            max_strength = strengths_left[s];
            idx_max = s;
        }
        else if (strengths_left[s] > second_max_strength){
            second_max_strength = strengths_left[s];
            idx_second_max = s;
        }
    }

    if (second_max_strength-1 > 0.2*(max_strength-1) || second_max_strength > 4){
        this->haploid = false;
    }
    else{
        if (idx_max != -1){
            this->consensus_left = neighbors_left[idx_max].first;
        }
        else{
            this->consensus_left = {};
        }
    }

    //go through the nighbors right and discard all reads that are subreads of another read //already done in add_neighbor
    // vector<bool> to_discard_right = vector<bool>(neighbors_right.size(), false);
    // for (int p = 0 ; p < neighbors_right.size() ; p++){
    //     for(int p2 = 0 ; p2 < p ; p2++){
    //         if (!to_discard_right[p2]){
    //             bool is_subread = true;
    //             for (int c = 0 ; c < neighbors_right[p].size() ; c++){
    //                 if (c >= neighbors_right[p2].size()){
    //                     to_discard_right[p2] = true;
    //                     is_subread = false;
    //                     break;
    //                 }
    //                 else if (neighbors_right[p2][c] != neighbors_right[p][c]){
    //                     is_subread = false;
    //                     break;
    //                 }
    //             }
    //             if (is_subread){
    //                 to_discard_right[p] = true;
    //                 break;
    //             }
    //         }
    //     }
    // }

    //go through the remaining reads: if there are more than one strong then not haploid, else set the only strong read as the consensus
    vector<int> strengths_right;
    for (auto neighbor : neighbors_right){
        strengths_right.push_back(neighbor.second);
    }
    //find the two strongest neighbors
    max_strength = 0;
    idx_max = -1;
    second_max_strength = 0;
    idx_second_max = -1;
    for (int s = 0 ; s < strengths_right.size() ; s++){
        if (strengths_right[s] > max_strength){
            second_max_strength = max_strength;
            idx_second_max = idx_max;
            max_strength = strengths_right[s];
            idx_max = s;
        }
        else if (strengths_right[s] > second_max_strength){
            second_max_strength = strengths_right[s];
            idx_second_max = s;
        }
    }

    if ((second_max_strength-1) > 0.2*(max_strength-1) || second_max_strength > 4){
        this->haploid = false;
    }
    else{
        if (idx_max != -1){
            this->consensus_right = neighbors_right[idx_max].first;
        }
        else{
            this->consensus_right = {};
        }
    }
}

vector<vector<pair<int,bool>>> Segment::get_strong_neighbors_left(){
    vector<vector<pair<int,bool>>> strong_neighbors;
    int max_strength = 0;
    for (auto neighbor : neighbors_left){
        if (neighbor.second > max_strength){
            max_strength = neighbor.second;
        }
    }
    for (auto neighbor : neighbors_left){
        if (neighbor.second > 0.2*max_strength || neighbor.second > 4){
            strong_neighbors.push_back(neighbor.first);
        }
    }
    return strong_neighbors;
}

vector<vector<pair<int,bool>>> Segment::get_strong_neighbors_right(){
    vector<vector<pair<int,bool>>> strong_neighbors;
    int max_strength = 0;
    for (auto neighbor : neighbors_right){
        if (neighbor.second > max_strength){
            max_strength = neighbor.second;
        }
    }
    for (auto neighbor : neighbors_right){
        if (neighbor.second > 0.2*max_strength || neighbor.second > 4){
            strong_neighbors.push_back(neighbor.first);
        }
    }
    return strong_neighbors;
}

/**
 * @brief Determine which contigs are haploid.
 * 
 * @param segments 
 * @param min_coverage 
 */
void determine_haploid_contigs(vector<Segment> &segments, int min_coverage){
    //determine which contigs are haploid. For this, a simple metric: the paths should be consensual left and right
    for (auto s = 0 ; s < segments.size() ; s++){
        segments[s].compute_consensuses();
        // if (segments[s].name == "4953--1"){
        //     cout << "neighbors left: " << endl;
        //     for (pair<vector<pair<int,bool>>,int> neighbor : segments[s].get_neighbors_left()){
        //         for (pair<int,bool> contig : neighbor.first){
        //             cout << segments[contig.first].name << " ";
        //         }
        //         cout << " strength: " << neighbor.second << endl;
        //     }
        //     cout << "neighbors right: " << endl;
        //     for (pair<vector<pair<int,bool>>,int> neighbor : segments[s].get_neighbors_right()){
        //         for (pair<int,bool> contig : neighbor.first){
        //             cout << segments[contig.first].name << " ";
        //         }
        //         cout << " strength: " << neighbor.second << endl;
        //     }
        //     cout << "consensus left: " << endl;
        //     for (pair<int,bool> contig : segments[s].get_consensus_left()){
        //         cout << segments[contig.first].name << " ";
        //     }
        //     cout << endl;
        //     cout << "consensus right: " << endl;
        //     for (pair<int,bool> contig : segments[s].get_consensus_right()){
        //         cout << segments[contig.first].name << " ";
        //     }
        //     cout << endl;
        //     cout << "haploid: " << segments[s].is_haploid() << endl;
        // }
    }    
}

/**
 * @brief Unzip the graph and create a new list of segments. 
 * 
 * @param old_segments 
 * @param new_segments 
 * @param min_coverage 
 */
void create_haploid_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, int min_coverage){

    new_segments = {};

     //associates all old contigs ids to their new IDs (can be multiple if contig is duplicated)
    int new_ID = 0;

    //begin by building bridges between haploid contigs
    for (auto old_segment = 0 ; old_segment < old_segments.size() ; old_segment++){

        // if (old_segments[old_segment].name != "2" && old_segments[old_segment].name != "666603"){
        //     continue;
        // }
        
        if (old_segments[old_segment].is_haploid()){

            // cout << "in haploid contig " << old_segments[old_segment].name << endl;

            //first create the contig if needed
            if (old_ID_to_new_IDs.find(old_segment) == old_ID_to_new_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_segment].name + "_0" , new_ID, old_segments[old_segment].get_pos_in_file(), old_segments[old_segment].get_coverage()));
                old_ID_to_new_IDs[old_segment] = {new_ID};
                new_ID++;
            }
            int ID_of_new_contig = old_ID_to_new_IDs[old_segment][0];

            //now see if it already has neighbors left and right and create them until the next haploid contig if not

            //left
            vector<pair<int,bool>> cons_left = old_segments[old_segment].get_consensus_left();
            bool there_is_a_bridge = false;
            double coverage_bridge = old_segments[old_segment].get_coverage();
            for (auto contig_and_orientation : cons_left){
                if (old_segments[contig_and_orientation.first].is_haploid()){
                    there_is_a_bridge = true;
                    coverage_bridge = (old_segments[contig_and_orientation.first].get_coverage() + old_segments[old_segment].get_coverage())/2;
                }
            }
            if (there_is_a_bridge){
                int contig_in_cons = 0;
                int previous_ID = ID_of_new_contig;
                int previous_old_ID = old_segment;
                int previous_end = 0;
                int ID_of_new_contig_left;
                while (true){

                    int end_of_contig_left = 0;
                    if (!cons_left[contig_in_cons].second){
                        end_of_contig_left = 1;
                    }

                    //create the contig if needed
                    if (old_ID_to_new_IDs.find(cons_left[contig_in_cons].first) == old_ID_to_new_IDs.end()){
                        new_segments.push_back(Segment(old_segments[cons_left[contig_in_cons].first].name + "_0" , new_ID, old_segments[cons_left[contig_in_cons].first].get_pos_in_file(), coverage_bridge));
                        old_ID_to_new_IDs[cons_left[contig_in_cons].first] = {new_ID};
                        new_ID++;
                    }
                    else{
                        //check if the link already exists
                        bool link_exists = false;
                        int next_ID;
                        for (int potential_valid_new_ids : old_ID_to_new_IDs[cons_left[contig_in_cons].first]){
                            for (pair<int,int> link : new_segments[potential_valid_new_ids].links[end_of_contig_left].first){
                                if (link.first == previous_ID && link.second == previous_end){
                                    link_exists = true;
                                    next_ID = potential_valid_new_ids;
                                }
                            }
                        }

                        if (link_exists){
                            //go to the next contig
                            if (old_segments[cons_left[contig_in_cons].first].is_haploid()){
                                break;
                            }
                            previous_end = 1-end_of_contig_left;
                            previous_ID = next_ID;
                            previous_old_ID = cons_left[contig_in_cons].first;
                            contig_in_cons++;
                            continue;
                        }
                        else if (!old_segments[cons_left[contig_in_cons].first].is_haploid()){
                            new_segments.push_back(Segment(old_segments[cons_left[contig_in_cons].first].name + "_"+std::to_string(old_ID_to_new_IDs[cons_left[contig_in_cons].first].size()) , new_ID, old_segments[cons_left[contig_in_cons].first].get_pos_in_file(), coverage_bridge));
                            old_ID_to_new_IDs[cons_left[contig_in_cons].first].push_back(new_ID);
                            new_ID++;
                        }
                    }
                    
                    old_segments[cons_left[contig_in_cons].first].decrease_coverage(coverage_bridge);

                    int ID_of_new_contig_left = old_ID_to_new_IDs[cons_left[contig_in_cons].first][old_ID_to_new_IDs[cons_left[contig_in_cons].first].size() - 1];

                    //add the link
                    //find the CIGAR in the old segment
                    string cigar = "*";
                    int idx = 0;
                    for (pair<int,int> link : old_segments[cons_left[contig_in_cons].first].links[end_of_contig_left].first){
                        if (link.first == previous_old_ID && link.second == previous_end){
                            cigar = old_segments[cons_left[contig_in_cons].first].links[end_of_contig_left].second[idx];
                        }
                        idx++;
                    }
                    new_segments[previous_ID].links[previous_end].first.push_back({ID_of_new_contig_left, end_of_contig_left});
                    new_segments[previous_ID].links[previous_end].second.push_back(cigar);
                    new_segments[ID_of_new_contig_left].links[end_of_contig_left].first.push_back({previous_ID, previous_end});
                    new_segments[ID_of_new_contig_left].links[end_of_contig_left].second.push_back(cigar);

                    //go to the next contig
                    if (old_segments[cons_left[contig_in_cons].first].is_haploid()){
                        break;
                    }
                    previous_end = 1-end_of_contig_left;
                    previous_ID = ID_of_new_contig_left;
                    previous_old_ID = cons_left[contig_in_cons].first;
                    contig_in_cons++;

                }
            }

            //right
            vector<pair<int,bool>> cons_right = old_segments[old_segment].get_consensus_right();
            there_is_a_bridge = false;
            coverage_bridge = old_segments[old_segment].get_coverage();
            for (auto contig_and_orientation : cons_right){
                if (old_segments[contig_and_orientation.first].is_haploid()){
                    there_is_a_bridge = true;
                    coverage_bridge = (old_segments[contig_and_orientation.first].get_coverage() + old_segments[old_segment].get_coverage())/2;
                }
            }
            if(there_is_a_bridge){

                int contig_in_cons = 0;
                int previous_ID = ID_of_new_contig;
                int previous_old_ID = old_segment;
                int previous_end = 1;
                while (true){

                    int end_of_contig_right = 0;
                    if (!cons_right[contig_in_cons].second){
                        end_of_contig_right = 1;
                    }

                    //create the contig if needed
                    if (old_ID_to_new_IDs.find(cons_right[contig_in_cons].first) == old_ID_to_new_IDs.end()){
                        new_segments.push_back(Segment(old_segments[cons_right[contig_in_cons].first].name + "_0" , new_ID, old_segments[cons_right[contig_in_cons].first].get_pos_in_file(), coverage_bridge));
                        old_ID_to_new_IDs[cons_right[contig_in_cons].first] = {new_ID};
                        new_ID++;
                    }
                    else{
                        //check if the link already exists
                        bool link_exists = false;
                        int next_ID;
                        for (int potential_valid_new_ids : old_ID_to_new_IDs[cons_right[contig_in_cons].first]){
                            for (pair<int,int> link : new_segments[potential_valid_new_ids].links[end_of_contig_right].first){
                                if (link.first == previous_ID && link.second == previous_end){
                                    link_exists = true;
                                    next_ID = potential_valid_new_ids;
                                }
                            }
                        }

                        if (link_exists){
                            //go to the next contig
                            if (old_segments[cons_right[contig_in_cons].first].is_haploid()){
                                break;
                            }
                            previous_end = 1-end_of_contig_right;
                            previous_ID = next_ID;
                            previous_old_ID = cons_right[contig_in_cons].first;
                            contig_in_cons++;
                            continue;
                        }
                        else if (!old_segments[cons_right[contig_in_cons].first].is_haploid()){
                            new_segments.push_back(Segment(old_segments[cons_right[contig_in_cons].first].name + "_"+std::to_string(old_ID_to_new_IDs[cons_right[contig_in_cons].first].size()) , new_ID, old_segments[cons_right[contig_in_cons].first].get_pos_in_file(), coverage_bridge));
                            old_ID_to_new_IDs[cons_right[contig_in_cons].first].push_back(new_ID);
                            new_ID++;
                        }
                    }

                    old_segments[cons_right[contig_in_cons].first].decrease_coverage(coverage_bridge);

                    int ID_of_new_contig_right = old_ID_to_new_IDs[cons_right[contig_in_cons].first][old_ID_to_new_IDs[cons_right[contig_in_cons].first].size() - 1];

                    //add the link
                    //find the CIGAR in the old segment
                    string cigar = "*";
                    int idx = 0;
                    for (auto link : old_segments[cons_right[contig_in_cons].first].links[end_of_contig_right].first){
                        if (link.first == previous_old_ID && link.second == previous_end){
                            cigar = old_segments[cons_right[contig_in_cons].first].links[end_of_contig_right].second[idx];
                        }
                        idx++;
                    }
                    new_segments[previous_ID].links[previous_end].first.push_back({ID_of_new_contig_right, end_of_contig_right});
                    new_segments[previous_ID].links[previous_end].second.push_back(cigar);
                    new_segments[ID_of_new_contig_right].links[end_of_contig_right].first.push_back({previous_ID, previous_end});
                    new_segments[ID_of_new_contig_right].links[end_of_contig_right].second.push_back(cigar);

                    //go to the next contig
                    if (old_segments[cons_right[contig_in_cons].first].is_haploid()){
                        break;
                    }
                    previous_end = 1-end_of_contig_right;
                    previous_ID = ID_of_new_contig_right;
                    previous_old_ID = cons_right[contig_in_cons].first;
                    contig_in_cons++;

                }
            }

        }

    }
}

/**
 * @brief Check what is not seen yet in new_segments
 * 
 * @param old_segments 
 * @param new_segments 
 * @param old_ID_to_new_IDs 
 * @param min_coverage 
 * @param all_paths
 * @return vector<vector<pair<int,bool>>> all the non represented paths
 */
vector<vector<pair<int,bool>>> list_non_represented_paths(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, int min_coverage){

    //first index represented paths
    vector<vector<pair<int,bool>>> represented_paths;
    unordered_map<int, vector<pair<int,int>>> where_is_this_contig_represented; //associates to an ID (indices of the path, position in the path)

    //go through all the haploids segments and list the represented paths left and right in new_segments
    for (auto old_segment = 0 ; old_segment < old_segments.size() ; old_segment++){
        if (old_segments[old_segment].is_haploid()){
            
            std::vector<std::pair<int,bool>> consensus_left = {{old_segments[old_segment].ID, false}};
            std::vector<std::pair<int,bool>> new_elements = old_segments[old_segment].get_consensus_left();
            consensus_left.insert(consensus_left.end(), new_elements.begin(), new_elements.end());

            //check if it goes until another haploid contig
            bool there_is_a_bridge = false;
            int idx_of_last_haploid_contig = 0;
            for (auto contig_and_orientation : consensus_left){
                if (old_segments[contig_and_orientation.first].is_haploid() && idx_of_last_haploid_contig != 0){
                    there_is_a_bridge = true;
                    break;
                }
                idx_of_last_haploid_contig++;
            }
            if (there_is_a_bridge){
                represented_paths.push_back(vector<pair<int,bool>>(consensus_left.begin(), consensus_left.begin() + idx_of_last_haploid_contig + 1));

                //fill where_is_this_contig_represented
                int number_of_haploid_contigs = 0;
                int idx = 0;
                for (auto contig_and_orientation : consensus_left){

                    if (where_is_this_contig_represented.find(contig_and_orientation.first) == where_is_this_contig_represented.end()){
                        where_is_this_contig_represented[contig_and_orientation.first] = {};
                    }
                    where_is_this_contig_represented[contig_and_orientation.first].push_back({represented_paths.size() - 1, idx});

                    if (old_segments[contig_and_orientation.first].is_haploid()){
                        number_of_haploid_contigs++;
                        if (number_of_haploid_contigs > 1){
                            break;
                        }
                    }
                    idx++;
                }
            }

            vector<pair<int,bool>> consensus_right = {{old_segments[old_segment].ID, true}};
            new_elements = old_segments[old_segment].get_consensus_right();
            consensus_right.insert(consensus_right.end(), new_elements.begin(), new_elements.end());

            //check if it goes until another haploid contig
            there_is_a_bridge = false;
            idx_of_last_haploid_contig = 0;
            for (auto contig_and_orientation : consensus_right){
                if (old_segments[contig_and_orientation.first].is_haploid() && idx_of_last_haploid_contig != 0){
                    there_is_a_bridge = true;
                    break;
                }
                idx_of_last_haploid_contig++;
            }

            if (there_is_a_bridge){

                represented_paths.push_back(vector<pair<int,bool>>(consensus_right.begin(), consensus_right.begin() + idx_of_last_haploid_contig + 1));

                //fill where_is_this_contig_represented
                int number_of_haploid_contigs = 0;
                int idx = 0;
                for (auto contig_and_orientation : consensus_right){

                    if (where_is_this_contig_represented.find(contig_and_orientation.first) == where_is_this_contig_represented.end()){
                        where_is_this_contig_represented[contig_and_orientation.first] = {};
                    }
                    where_is_this_contig_represented[contig_and_orientation.first].push_back({represented_paths.size() - 1, idx});

                    if (old_segments[contig_and_orientation.first].is_haploid()){
                        number_of_haploid_contigs++;
                        if (number_of_haploid_contigs > 1){
                            break;
                        }
                    }
                    idx++;
                }
            }
        }
    }

    // cout << "finished indexing represented paths" << endl;
    // cout << "for example, here are all the paths where I find 1: " << endl;
    // for (pair<int,int> path_and_pos : where_is_this_contig_represented[1]){
    //     cout << "path " << path_and_pos.first << " at position " << path_and_pos.second << " : ";
    //     for (pair<int,bool> contig : represented_paths[path_and_pos.first]){
    //         cout << contig.first << " ";
    //     }
    //     cout << endl;
    // }

    //now go through all the paths and see if they are represented
    vector<vector<pair<int,bool>>> unrepresented_paths;
    int idx = 0;
    for (Segment s : old_segments){
        vector<vector<pair<int,bool>>> paths_segment = s.get_strong_neighbors_left();

        int number_of_left_paths = paths_segment.size();
        auto neighbors_right = s.get_strong_neighbors_right();

        // if (s.name == "4953--1"){
        //     cout << "strong neighbors left: " << endl;
        //     for (vector<pair<int,bool>> neighbor : paths_segment){
        //         for (pair<int,bool> contig : neighbor){
        //             cout << old_segments[contig.first].name << " ";
        //         }
        //         cout << endl;
        //     }
        //     cout << "strong neighbors right: " << endl;
        //     for (vector<pair<int,bool>> neighbor : neighbors_right){
        //         for (pair<int,bool> contig : neighbor){
        //             cout << old_segments[contig.first].name << " ";
        //         }
        //         cout << endl;
        //     }
        // }

        paths_segment.insert(paths_segment.end(), neighbors_right.begin(), neighbors_right.end());

        int idx_here = 0;
        for (vector<pair<int,bool>> path : paths_segment){

            if (idx % 100000 == 0){
                cout << idx << " paths checked\r" << std::flush;
            }

            //check if the path is represented
            vector<pair<int,bool>> path_until_haploid_contig;
            if (idx_here < number_of_left_paths){
                path_until_haploid_contig = {{s.ID, false}};
            }
            else{
                path_until_haploid_contig = {{s.ID, true}};
            }
            for (pair<int,bool> contig : path){

                if (!old_segments[contig.first].is_haploid()){
                    path_until_haploid_contig.push_back(contig);
                }
                else{
                    path_until_haploid_contig.push_back(contig);
                    bool found = false;
                    // if (s.name == "4953--1"){
                    //     cout << "checking if path " << idx_here << " ";
                    //     for (pair<int,bool> contig : path_until_haploid_contig){
                    //         cout << old_segments[contig.first].name << " ";
                    //     }
                    //     cout << " is represented" << endl;
                    // }
                    if (where_is_this_contig_represented.find(contig.first) != where_is_this_contig_represented.end()){
                        for (pair<int,int> path_and_pos : where_is_this_contig_represented[contig.first]){
                            vector<pair<int,bool>> path_to_check;
                            if (path_and_pos.second == 0){
                                path_to_check = vector<pair<int,bool>>(represented_paths[path_and_pos.first].begin(), represented_paths[path_and_pos.first].begin()+std::min(path_until_haploid_contig.size(),represented_paths[path_and_pos.first].size()));
                            }
                            else {
                                path_to_check = vector<pair<int,bool>>(represented_paths[path_and_pos.first].end()-std::min(path_until_haploid_contig.size(),represented_paths[path_and_pos.first].size()), represented_paths[path_and_pos.first].end());
                            }
                            // if (s.name == "4953--1"){
                            //     cout << "comparing to ";
                            //     for (pair<int,bool> contig : path_to_check){
                            //         cout << old_segments[contig.first].name << " " << contig.second << " ; ";
                            //     }
                            //     cout << endl;
                            // }
                            if (path_until_haploid_contig == path_to_check){
                                found = true;
                            }
                        }
                    }
                    if (!found){
                        // if (s.name == "4953--1"){
                        //     cout << "unrepresented path: ";
                        //     for (pair<int,bool> contig : path_until_haploid_contig){
                        //         cout << old_segments[contig.first].name << " " << contig.second << " ; ";
                        //     }
                        //     cout << endl;
                        // }
                        unrepresented_paths.push_back(path_until_haploid_contig);
                    }
                    // else{
                    //     cout << "path found: ";
                    //     for (pair<int,bool> contig : path_until_haploid_contig){
                    //         cout << contig.first << " ";
                    //     }
                    //     cout << endl;
                    // }
                    path_until_haploid_contig = {contig};
                }
            }
            //check the last path
            if (path_until_haploid_contig.size() > 0){
                if (where_is_this_contig_represented.find(path_until_haploid_contig.back().first) == where_is_this_contig_represented.end()){
                    // is_represented = false;
                }
                bool found = false;
                // if (s.name == "4953--1"){
                //     cout << "checking if paeeth " << idx_here << " ";
                //     for (pair<int,bool> contig : path_until_haploid_contig){
                //         cout << old_segments[contig.first].name << " ";
                //     }
                //     cout << " is represented" << endl;
                // }
                for (pair<int,int> path_and_pos : where_is_this_contig_represented[path_until_haploid_contig.back().first]){
                    auto path_to_check = vector<pair<int,bool>>(represented_paths[path_and_pos.first].begin(), represented_paths[path_and_pos.first].begin()+std::min(path_until_haploid_contig.size(),represented_paths[path_and_pos.first].size()));
                    // cout << "comparing to ";
                    // for (pair<int,bool> contig : path_to_check){
                    //     cout << contig.first << " ";
                    // }
                    // cout << endl;
                    if (path_until_haploid_contig == path_to_check){
                        found = true;
                    }
                }
                if (!found){
                    unrepresented_paths.push_back(path_until_haploid_contig);
                    // if (s.name == "4953--1"){
                    //     cout << "unrepresented path: ";
                    //     for (pair<int,bool> contig : path_until_haploid_contig){
                    //         cout << old_segments[contig.first].name << " " << contig.second << " ; ";
                    //     }
                    //     cout << endl;
                    // }
                }
            }
            idx_here++;
            idx++;
        }
    }

    return unrepresented_paths;

}

void add_unrepresented_paths(vector<Segment> &old_segments, vector<Segment> &new_segments, unordered_map<int, std::vector<int>>& old_ID_to_new_IDs, int min_coverage, vector<vector<pair<int,bool>>>& unrepresented_paths){

    unordered_map<int, int> old_IDs_to_new_non_haploid_IDs; //associates old IDs to new IDs for the contig we are going to create
    for (Segment s: old_segments){ //we are not going to create new versions of haploid contigs
        if (s.is_haploid()){
            old_IDs_to_new_non_haploid_IDs[s.ID] = old_ID_to_new_IDs[s.ID][0];
        }
    }

    //convert unrepresented path in a list of links that must be there in the final graph
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add;
    int index = 0;

    for (vector<pair<int,bool>> path : unrepresented_paths){

        // cout << "adding the contigs of path " << endl;
        // for (pair<int,bool> contig : path){
        //     cout << contig.first << " ";
        // }
        // cout << endl;

        for (int contig = 0 ; contig < path.size() - 1 ; contig++){

            int old_ID1 = path[contig].first;
            int old_ID2 = path[contig+1].first;

            int end1 = 1;
            int end2 = 0;
            if (!path[contig].second){
                end1 = 0;
            }
            if (!path[contig+1].second){
                end2 = 1;
            }

            //create the contigs if not already done
            if (old_IDs_to_new_non_haploid_IDs.find(old_ID1) == old_IDs_to_new_non_haploid_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_ID1].name + "_" + std::to_string(old_ID_to_new_IDs[old_ID1].size()) , new_segments.size(), old_segments[old_ID1].get_pos_in_file(), old_segments[old_ID1].get_coverage()));
                old_IDs_to_new_non_haploid_IDs[old_ID1] = new_segments.size() - 1;
                if (old_ID_to_new_IDs.find(old_ID1) == old_ID_to_new_IDs.end()){
                    old_ID_to_new_IDs[old_ID1] = {(int) new_segments.size() - 1};
                }
                else{
                    old_ID_to_new_IDs[old_ID1].push_back(new_segments.size() - 1);
                }
            }

            if (old_IDs_to_new_non_haploid_IDs.find(old_ID2) == old_IDs_to_new_non_haploid_IDs.end()){
                new_segments.push_back(Segment(old_segments[old_ID2].name + "_" + std::to_string(old_ID_to_new_IDs[old_ID2].size()) , new_segments.size(), old_segments[old_ID2].get_pos_in_file(), old_segments[old_ID2].get_coverage()));
                old_IDs_to_new_non_haploid_IDs[old_ID2] = new_segments.size() - 1;
                if (old_ID_to_new_IDs.find(old_ID2) == old_ID_to_new_IDs.end()){
                    old_ID_to_new_IDs[old_ID2] = {(int) new_segments.size() - 1};
                }
                else{
                    old_ID_to_new_IDs[old_ID2].push_back(new_segments.size() - 1);
                }
            }

            int new_ID1 = old_ID_to_new_IDs[old_ID1][old_ID_to_new_IDs[old_ID1].size()-1];
            int new_ID2 = old_ID_to_new_IDs[old_ID2][old_ID_to_new_IDs[old_ID2].size()-1];
            string cigar = "*";
            int idx = 0;
            for (pair<int,int> link : old_segments[old_ID2].links[end2].first){
                if (link.first == old_ID1 && link.second == end1){
                    cigar = old_segments[old_ID2].links[end2].second[idx];
                }
                idx ++;
            }

            if (links_to_add.find({{pair<int,int>(new_ID2, end2), pair<int,int>(new_ID1, end1)}, cigar}) == links_to_add.end()){
                links_to_add.insert({{pair<int,int>(new_ID1, end1), pair<int,int>(new_ID2, end2)}, cigar});
            }
        }
        // if (index > 3){
        //     cout << "breaking" << endl;
        //     break;
        // }
        index++;
    }

    //add the links
    for (pair<pair<pair<int,int>, pair<int,int>>, string> link : links_to_add){
        new_segments[link.first.first.first].links[link.first.first.second].first.push_back({link.first.second.first, link.first.second.second});
        new_segments[link.first.first.first].links[link.first.first.second].second.push_back(link.second);
        new_segments[link.first.second.first].links[link.first.second.second].first.push_back({link.first.first.first, link.first.first.second});
        new_segments[link.first.second.first].links[link.first.second.second].second.push_back(link.second);
    }
}
    
void merge_adjacent_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, string original_gfa_file){

    set<int> already_looked_at_segments;
    int seg_idx = 0;
    unordered_map<pair<int,int>,pair<int,int>> old_ID_to_new_ID; //associates (old_id, old end) with (new_id_new_end)
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add; //list of links to add, all in old IDs
    for (Segment old_seg : old_segments){

        // if (old_seg.name != "1665420_0"){
        //     seg_idx++;
        //     continue;
        // }

        if (already_looked_at_segments.find(old_seg.ID) != already_looked_at_segments.end()){
            seg_idx++;
            continue;
        }
        //check if it has either at least two neighbors left or that its neighbor left has at least two neighbors right
        cout << "in merge, looking at segment " << seg_idx << " out of " << old_segments.size() << "\r" << std::flush;
        bool dead_end_left = false;
        if (old_seg.links[0].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].links[old_seg.links[0].first[0].second].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].ID == old_seg.ID){
            dead_end_left = true;
        }

        bool dead_end_right = false;
        if (old_seg.links[1].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].links[old_seg.links[1].first[0].second].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].ID == old_seg.ID){
            dead_end_right = true;
        }

        // cout << "dead end left: " << dead_end_left << " dead end right: " << dead_end_right << endl;

        if (!dead_end_left && !dead_end_right){ //means this contig is in the middle of a long haploid contig, no need to merge
            seg_idx++;
            continue;
        }

        //create a new contig
        if (dead_end_left && dead_end_right){
            new_segments.push_back(Segment(old_seg.name, new_segments.size(), old_seg.get_pos_in_file(), old_seg.get_coverage()));
            old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 1};
            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[0].first){
                links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_seg.links[1].first){
                links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                idx_link++;
            }
            // cout << "yupee" << endl;
        }
        else if (dead_end_left){
            
            //let's see how far we can go right
            vector<string> all_names = {old_seg.name};
            vector<string> all_seqs = {old_seg.get_seq(original_gfa_file)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            int current_ID = old_seg.ID;
            int current_end = 1;

            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                already_looked_at_segments.insert(current_ID);
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 1
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
            }

            //create the new contig
            string new_name = "";
            for (string name : all_names){
                new_name += name + "_";
            }
            new_name = new_name.substr(0, new_name.size()-1);
            string new_seq = "";
            for (string seq : all_seqs){
                new_seq += seq;
            }
            double new_coverage = 0;
            for (double coverage : all_coverages){
                new_coverage += coverage;
            }
            new_coverage = new_coverage/all_coverages.size();
            new_segments.push_back(Segment(new_name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
            new_segments[new_segments.size()-1].seq = new_seq;

            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[0].first){
                links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                links_to_add.insert({{{current_ID, current_end}, link}, old_segments[current_ID].links[current_end].second[idx_link]});
                idx_link++;
            }

            old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

            // cout << "hurray" << endl;
        }
        else{
            
            //let's see how far we can go left
            vector<string> all_names = {old_seg.name};
            string seq = old_seg.get_seq(original_gfa_file);
            vector<string> all_seqs = {reverse_complement(seq)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            int current_ID = old_seg.ID;
            int current_end = 0;

            // cout << "exploring all the contigs left" << endl;
            // cout << "first exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name << endl;
            
            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                already_looked_at_segments.insert(current_ID);
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 0
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());

                // cout << "exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name  << " " << old_segments[current_ID].links[current_end].first.size() << endl;
                // cout << "here are all the links : ";
                // for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                //     cout << link.first << " " << old_segments[link.first].name << " " << link.second << " ; ";
                // }
                // cout << endl;
            }

            //create the new contig
            string new_name = "";
            for (string name : all_names){
                new_name += name + "_";
            }
            new_name = new_name.substr(0, new_name.size()-1);
            string new_seq = "";
            for (string seq : all_seqs){
                new_seq += seq;
            }
            double new_coverage = 0;
            for (double coverage : all_coverages){
                new_coverage += coverage;
            }
            new_coverage = new_coverage/all_coverages.size();
            new_segments.push_back(Segment(new_name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
            new_segments[new_segments.size()-1].seq = new_seq;

            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[1].first){
                links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                links_to_add.insert({{{current_ID, current_end}, link}, old_segments[current_ID].links[current_end].second[idx_link]});
                idx_link++;
            }

            old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

            // cout << "yay" << endl;
        }

        seg_idx++;
    }

    //now add the links in the new segments
    for (pair<pair<pair<int,int>, pair<int,int>>, string> link : links_to_add){
        new_segments[old_ID_to_new_ID[link.first.first].first].links[old_ID_to_new_ID[link.first.first].second].first.push_back({old_ID_to_new_ID[link.first.second].first, link.first.second.second});
        new_segments[old_ID_to_new_ID[link.first.first].first].links[old_ID_to_new_ID[link.first.first].second].second.push_back(link.second);
    }
}

void output_graph(string gfa_output, string gfa_input, vector<Segment> &segments){
    ofstream gfa(gfa_output);
    for (Segment s : segments){
        gfa << "S\t" << s.name << "\t" << s.get_seq(gfa_input) << "\tDP:f:" << s.get_coverage() <<  "\n";
    }
    for (Segment s : segments){
        for (int end = 0 ; end < 2 ; end++){
            for (int neigh = 0 ; neigh < s.links[end].first.size() ; neigh++){

                //to make sure the link is not outputted twice
                if (s.ID > s.links[end].first[neigh].first || (s.ID == s.links[end].first[neigh].first && end > s.links[end].first[neigh].second) ){
                    continue;
                }

                string orientation = "+";
                if (end == 0){
                    orientation = "-";
                }
                gfa << "L\t" << s.name << "\t" << orientation << "\t" << segments[s.links[end].first[neigh].first].name << "\t";
                if (s.links[end].first[neigh].second == 0){
                    gfa << "+\t";
                }
                else{
                    gfa << "-\t";
                }
                gfa << s.links[end].second[neigh] << "\n";
            }
        }
    }
    gfa.close();
}

int main(int argc, char *argv[])
{
    //HS_GraphUnzip <gfa_input> <gaf_file> <threads> <gfa_output> <exhaustive>
    if (argc != 8){
        std::cout << "Usage: HS_GraphUnzip <gfa_input> <gaf_file> <min_coverage> <threads> <gfa_output> <exhaustive> <logfile>" << std::endl;
        return 1;
    }

    std::string gfa_input = argv[1];
    std::string gaf_file = argv[2];
    int min_coverage = std::stoi(argv[3]);
    int threads = std::stoi(argv[4]);
    std::string gfa_output = argv[5];
    bool exhaustive = std::stoi(argv[6]);
    std::string logfile = argv[7];

    ofstream log(logfile);

    //load the segments from the GFA file
    unordered_map<string, int> segment_IDs;
    vector<Segment> segments; //segments is a dict of pairs
    load_GFA(gfa_input, segments, segment_IDs);
    cout << "Segments loaded" << endl;
    log << "Segments loaded" << endl;

    //load the paths from the GAF file
    load_GAF(gaf_file, segments, segment_IDs);
    cout << "Paths loaded" << endl;
    log << "Paths loaded" << endl;

    //unzip the graph
    cout << "Now unzipping the graph" << endl;
    vector<Segment> unzipped_segments;
    unordered_map<int, std::vector<int>> old_ID_to_new_IDs;
    determine_haploid_contigs(segments, min_coverage);
    cout << "Haploid contigs determined" << endl;
    cout << "Creating haploid contigs" << endl;
    create_haploid_contigs(segments, unzipped_segments, old_ID_to_new_IDs, min_coverage);
    output_graph("haploid.gfa", gfa_input, unzipped_segments);
    cout << "Haploid contigs created" << endl;
    cout << "Listing non represented paths" << endl;
    vector<vector<pair<int,bool>>> non_represented_paths = list_non_represented_paths(segments, unzipped_segments, old_ID_to_new_IDs, min_coverage);
    cout << "Non represented paths listed" << endl;
    cout << "Adding non represented paths" << endl;
    add_unrepresented_paths(segments, unzipped_segments, old_ID_to_new_IDs, min_coverage, non_represented_paths);
    cout << "Non represented paths added" << endl;

    cout << "Graph unzipped" << endl;

    cout << "Merging adjacent contigs" << endl;
    vector<Segment> merged_segments;
    merge_adjacent_contigs(unzipped_segments, merged_segments, gfa_input);
    // merged_segments = unzipped_segments;
    cout << "Adjacent contigs merged" << endl;

    //output the graph
    cout << "Outputting the graph" << endl;
    output_graph(gfa_output, gfa_input, merged_segments);
    cout << "Graph outputted" << endl;

    //TODO: consensus tolerating one error
}