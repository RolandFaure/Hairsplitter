#include "graphunzip.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <Eigen/Sparse>

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

struct Segment{
    string name;
    int ID;
    vector<set<pair<int,int>> > links = vector<set<pair<int,int>>>(2); //first vector is for the links to the left, second vector is for the links to the right. Each link is the ID of the neighbor and its end (0 for left, 1 for right)

    //the hash of the segment is the hash of the name
    size_t hash() const{
        return std::hash<string>{}(name);
    }
};


void load_GFA(string gfa_file, unordered_map<int, Segment> &segments, unordered_map<string, int> &segment_IDs){
    //load the segments from the GFA file
    
    //in a first pass index all the segments by their name
    ifstream gfa(gfa_file);
    string line;
    int ID = 0;
    while (getline(gfa, line)){
        if (line[0] == 'S'){
            stringstream ss(line);
            string nothing, name;
            ss >> nothing >> name;
            Segment s;
            s.name = name;
            s.ID = ID;
            ID += 1;
            segment_IDs[name] = s.ID;
            segments[s.ID] = s;
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

            ss >> nothing >> name1 >> orientation1 >> name2 >> orientation2;

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

            segments[ID1].links[end1].insert({ID2, end2});
            segments[ID1].links[end2].insert({ID2, end1});
        }
    }
    gfa.close();
}

void load_GAF(string gaf_file, unordered_map<int, Segment> &segments, unordered_map<string, int> &segments_IDs, Eigen::SparseMatrix<int>& contact_matrix){
    //load the paths from the GAF file
    ifstream gaf(gaf_file);
    string line;
    int nb_lines = 0;
    vector<Eigen::Triplet<int>> triplet_list;
    while (getline(gaf, line)){
        stringstream ss(line);
        string name_of_read;
        string path;
        string nothing;
        ss >> name_of_read >> nothing >> nothing >> nothing >> nothing >> nothing >> nothing >> path;

        if (nb_lines % 100000 == 0){
            cout << nb_lines << " lines read in load_GAF\r" << std::flush;
        }
        if (nb_lines > 3000000){
            cout << "NOT READING EVERYTHING FOR DEBGUGGING PUROSPS" << endl;
            break;
        }
        nb_lines += 1;
        
        bool orientation_now;
        string str_now;
        vector<int> segments_now;
        vector<bool> orientations_now;
        vector<pair<vector<int>, vector<bool>>> paths_in_this_read;
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
                    if (segments_now.size() == 0 || segments[segments_now.back()].links[orientations_now.back()].find({new_contig_ID, !orientation_now}) != segments[segments_now.back()].links[orientations_now.back()].end()){
                        segments_now.push_back(new_contig_ID);
                        orientations_now.push_back(orientation_now);
                    }
                    else{
                        //then the neighbor is not the neighbor of a previous segment, cut
                        paths_in_this_read.push_back({segments_now, orientations_now});
                        segments_now = {new_contig_ID};
                        orientations_now = {orientation_now};     
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
        paths_in_this_read.push_back({segments_now, orientations_now});

        //go through paths_in_this_read and fill the contact matrix 
        for (pair<vector<int>, vector<bool>> pair_path : paths_in_this_read){
            for (int contig1 = 0 ; contig1 < pair_path.first.size() ; contig1++){
                for (int contig2 = contig1 + 1 ; contig2 < pair_path.first.size() ; contig2++){
                    triplet_list.push_back(Eigen::Triplet<int>(pair_path.first[contig1],pair_path.first[contig2],1));
                    triplet_list.push_back(Eigen::Triplet<int>(pair_path.first[contig2],pair_path.first[contig1],1));
                }
            }
        }
    }
    gaf.close();

    contact_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());

    cout <<"here are all the contacs for contig 1951404: " << endl;
    for (int i = 0 ; i < contact_matrix.outerSize(); i++){
        for (Eigen::SparseMatrix<int>::InnerIterator it(contact_matrix,i); it; ++it){
            if (it.row() == 1951404){
                cout << it.row() << " " << it.col() << " " << it.value() << endl;
            }
        }
    }
}

void remove_unsupported_links(unordered_map<int, Segment> &segments, Eigen::SparseMatrix<int>& contact_matrix, bool careful){

    //now remove all the links that are not confirmed
    vector<pair<pair<int, int>, pair<int, int>>> links_to_remove;
    for (auto p : segments){

        
        //first look at all the links to the left
        for (auto l : p.second.links[0]){
            if (contact_matrix.coeff(p.second.ID, l.first) == 0){
                if (!careful || (p.second.links[0].size() > 1 && segments[l.first].links[l.second].size() > 1)){ //we can delete without creating dead ends
                    links_to_remove.push_back({{p.first, 0}, l});
                    // cout << "Removing link that goes from " << p.first << " to " << l.first << " at the left" << endl;
                }
            }
        }

        //then look at all the links to the right
        links_to_remove.clear();
        for (auto l : p.second.links[1]){
            if (contact_matrix.coeff(p.second.ID, l.first)){
                if (!careful || (p.second.links[1].size() > 1 && segments[l.first].links[l.second].size() > 1)){ //we can delete without creating dead ends
                    links_to_remove.push_back({{p.first, 1}, l});
                    // cout << "Removing link that goes from " << p.first << " to " << l.first << " at the right" << endl;
                }
            }
        }
    }

    //now remove the links
    for (auto l : links_to_remove){
        if (segments[l.first.first].links[l.first.second].find(l.second) != segments[l.first.first].links[l.first.second].end()){
            segments[l.first.first].links[l.first.second].erase(l.second);
        }
        if (segments[l.second.first].links[l.second.second].find(l.first) != segments[l.second.first].links[l.second.second].end()){
            segments[l.second.first].links[l.second.second].erase(l.first);
        }
    }
}

void unzip_graph(unordered_map<int, Segment> &segments, Eigen::SparseMatrix<int>& contact_matrix){

    set<int> potentially_interesting_contigs; 
    for (auto s : segments){
        potentially_interesting_contigs.emplace(s.second.ID);
    }

    //go through the contigs and duplicate those that can be duplicated
    for (auto contig_and_contig_ID : segments){
        if (potentially_interesting_contigs.find(contig_and_contig_ID.first) == potentially_interesting_contigs.end()){
            continue;
        }

        if (contig_and_contig_ID.first != 1951404){
            continue;
        }

        //find left dilemma (the end of contig where path separate left) and right dilemma
        pair<int,int> left_dilemma = {contig_and_contig_ID.first , 0}; //the ID of the contig and its end
        cout << "startign with " << left_dilemma.first << endl;
        bool go_on = true;
        while(go_on){
            //list all serious contigs left
            vector<pair<int, int>> serious_neighbors;
            for (pair<int,int> neighbor_left : segments[left_dilemma.first].links[left_dilemma.second]){

                //check if there is serious contact with the contig
                if(contact_matrix.coeff(neighbor_left.first, contig_and_contig_ID.first) > 0){
                    serious_neighbors.push_back(neighbor_left);
                }

            }

            cout << "here is the list of serious neighbors" << endl;
            for (auto s : serious_neighbors){
                cout << s.first << " " << s.second << endl;
            }

            if (serious_neighbors.size() != 1){
                go_on = false;
            }
            else{
                left_dilemma = {serious_neighbors[0].first , 1-serious_neighbors[0].second};
            }
        }

        pair<int,int> right_dilemma = {contig_and_contig_ID.first , 1}; //the ID of the contig and its end
        cout << "startign with " << right_dilemma.first << endl;
        go_on = true;
        while(go_on){
            //list all serious contigs right
            vector<pair<int, int>> serious_neighbors;
            for (pair<int,int> neighbor_right : segments[right_dilemma.first].links[right_dilemma.second]){

                //check if there is serious contact with the contig
                if(contact_matrix.coeff(neighbor_right.first, contig_and_contig_ID.first) > 0){
                    serious_neighbors.push_back(neighbor_right);
                }

            }
            cout << "here is the list of serious neighbors" << endl;
            for (auto s : serious_neighbors){
                cout << s.first << " " << s.second << endl;
            }
            if (serious_neighbors.size() != 1){
                go_on = false;
            }
            else{
                right_dilemma = {serious_neighbors[0].first , 1-serious_neighbors[0].second};
            }
        }

        cout << "right and left dilemmas of contig " << contig_and_contig_ID.first << " are " << left_dilemma.first << " and " << right_dilemma.first << endl;

    }

}


int main(int argc, char *argv[])
{
    //HS_GraphUnzip <gfa_input> <gaf_file> <threads> <gfa_output> <exhaustive>
    if (argc != 7){
        std::cout << "Usage: HS_GraphUnzip <gfa_input> <gaf_file> <threads> <gfa_output> <exhaustive> <logfile>" << std::endl;
        return 1;
    }

    std::string gfa_input = argv[1];
    std::string gaf_file = argv[2];
    int threads = std::stoi(argv[3]);
    std::string gfa_output = argv[4];
    bool exhaustive = std::stoi(argv[5]);
    std::string logfile = argv[6];

    ofstream log(logfile);

    //load the segments from the GFA file
    unordered_map<string, int> segment_IDs;
    unordered_map<int, Segment> segments; //segments is a dict of pairs (ID, segment)
    load_GFA(gfa_input, segments, segment_IDs);
    cout << "Segments loaded" << endl;
    log << "Segments loaded" << endl;

    //load the paths from the GAF file
    //create a sparse contact matrix of contact between the contigs
    Eigen::SparseMatrix<int> contact_matrix (segments.size() , segments.size());
    load_GAF(gaf_file, segments, segment_IDs, contact_matrix);
    cout << "Paths loaded" << endl;
    log << "Paths loaded" << endl;

    //remove the links that are not supported by the paths
    if (exhaustive){
        cout << "Now removing unsupported links" << endl;
        remove_unsupported_links(segments, contact_matrix, true);
        cout << "Unsupported links deleted" << endl;
        log << "Unsupported links deleted" << endl;
    }

    //unzip the graph
    cout << "Now unzipping the graph" << endl;
    unzip_graph(segments, contact_matrix);
    cout << "Graph unzipped" << endl;
}