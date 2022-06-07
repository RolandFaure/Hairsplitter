#include "modify_gfa.h"
#include<algorithm> //for "sort"
#include <fstream>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::list;
using std::vector;
using std::min;
using std::max;
using std::begin;
using std::end;
using std::unordered_map;
using std::pair;
using std::make_pair;
using std::set;
using std::to_string;

void modify_GFA(std::string refFile, std::vector <Read> &allreads, vector<unsigned long int> &backbones_reads, std::vector <Overlap> &allOverlaps,
            std::unordered_map<unsigned long int, std::vector< std::pair<std::pair<int,int>, std::vector<int>> >> &partitions, string outputFile, vector<Link> &allLinks,
            std::unordered_map <int, std::pair<int,int>> &clusterLimits){

    int max_backbone = backbones_reads.size(); //fix that because backbones will be added to the list but not separated 
    for (int b = 0 ; b < max_backbone ; b++){

        int backbone = backbones_reads[b];

        if (partitions.find(backbone) != partitions.end() && partitions[backbone].size() > 1){

            //stitch all intervals of each backbone read
            vector<unordered_map <int,set<int>>> stitches(partitions[backbone].size()); //aggregating the information from stitchesLeft and stitchesRight to know what link to keep

            for (int n = 0 ; n < partitions[backbone].size() ; n++){
                //for each interval, go through the different parts and see with what part before and after they fit best
                if (n > 0){
                    auto stitchLeft = stitch(partitions[backbone][n].second, partitions[backbone][n-1].second);

                    auto stitchRight = stitch(partitions[backbone][n-1].second, partitions[backbone][n].second);

                    for (auto s : stitchLeft){
                        stitches[n][s.first] = s.second;
                    }

                    for (auto s : stitchRight){
                        for (int neighbor : s.second){
                            stitches[n][neighbor].emplace(s.first);
                        }
                    }

                    //make sure all contigs are stitched
                    set <int> stitchedContigs;
                    for (auto s : stitches[n]){
                        for (int neighbor : s.second){
                            stitchedContigs.emplace(neighbor);
                        }
                    }
                    // cout << "stitched contigs : "; for(auto i : stitchedContigs) {cout <<i << ";";} cout << endl;
                    for (auto s : stitchRight){
                        if (stitchedContigs.find(s.first) == stitchedContigs.end()){
                            // cout << "gluing back" << endl;
                            for (auto s2 : stitchLeft){
                                stitches[n][s2.first].emplace(s.first);
                            }
                        }
                    }

                    // cout << "stitching : ";
                    // for (auto s: stitches[n]){
                    //     cout << s.first << "<->" ;
                    //     for (auto s2 : s.second) {cout << s2 << ",";}
                    //     cout << "  ;  ";
                    // } 
                    // cout << endl;
                }
            }


            //create hangingLinks, a list of links that are not yet connected but will soon be
            vector<int> hangingLinks;
            for (int linkIdx : allreads[backbone].get_links_left()){
                if (allLinks[linkIdx].neighbor1 == backbone){
                    allLinks[linkIdx].end1 = -1;
                }
                else {
                    allLinks[linkIdx].end2 = -1;
                }
                allLinks[linkIdx].group = 0;
                hangingLinks.push_back(linkIdx);
            }

            int n = 0;
            for (auto interval : partitions[backbone]){

                if (interval.first.second == interval.first.first){ //to avoid stupid segfaults
                    interval.first.second+= 1;
                }

                unordered_map<int, vector<string>> readsPerPart; //list of all reads of each part
                cout << "in interval " << interval.first.first << " <-> " << interval.first.second << endl;

                for (int r = 0 ; r < interval.second.size(); r++){
                    if (interval.second[r] != -1){
                        int clust = interval.second[r];
                        int limitLeft = allOverlaps[allreads[backbone].neighbors_[r]].position_1_1; //the limit of the read that we should use
                        int limitRight = allOverlaps[allreads[backbone].neighbors_[r]].position_1_2;
                        auto idxRead = allOverlaps[allreads[backbone].neighbors_[r]].sequence1;

                        if (readsPerPart.find(clust) == readsPerPart.end()){
                            readsPerPart[clust] = {allreads[idxRead].sequence_.subseq(limitLeft, limitRight-limitLeft).str()};
                        }
                        else {
                            readsPerPart[clust].push_back(allreads[idxRead].sequence_.subseq(limitLeft, limitRight-limitLeft).str());
                            // cout << "read aligning is : " << allreads[idxRead].name << endl;
                        }
                    }
                }

                string toPolish = allreads[backbone].sequence_.str().substr(interval.first.first, interval.first.second-interval.first.first);
                vector<int> futureHangingLinks;

                for (auto group : readsPerPart){
                    
                    string newcontig = "";
                    if (readsPerPart.size() > 1){
                        // newcontig = local_assembly(group.second);
                        if (newcontig == ""){//if the assembly was not successful for one reason or another
                            newcontig = consensus_reads(toPolish, group.second);
                        }
                        EdlibAlignResult result = edlibAlign(toPolish.c_str(), toPolish.size(),
                                    newcontig.c_str(), newcontig.size(),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

                        newcontig = newcontig.substr(result.startLocations[0], result.endLocations[0]-result.startLocations[0]);

                        edlibFreeAlignResult(result);
                        
                    }
                    else {
                        newcontig = toPolish;
                    }

                    Read r(newcontig);
                    r.name = allreads[backbone].name + "_"+ to_string(interval.first.first)+ "_" + to_string(group.first);

                    //now create all the links IF they are compatible with "stitches"  
                    set<int> linksToKeep;
        
                    if (n == 0 || stitches[n][group.first].size() == 0){
                        for (int h : hangingLinks){
                            linksToKeep.emplace(allLinks[h].group);
                        }
                    }
                    else{
                        for (int l : stitches[n][group.first]){
                            linksToKeep.emplace(l);
                        }
                    }
                    
                    //create the links
                    for (int h : hangingLinks){
                        if (linksToKeep.find(allLinks[h].group) != linksToKeep.end()){
                            Link leftLink;
                            leftLink.CIGAR = allLinks[h].CIGAR;
                            if (allLinks[h].end2 == -1){
                                leftLink.end2 = 0;
                                leftLink.neighbor2 = allreads.size();
                                leftLink.end1 = allLinks[h].end1;
                                leftLink.neighbor1 = allLinks[h].neighbor1;
                                allreads[leftLink.neighbor1].add_link(allLinks.size(), allLinks[h].end1);
                            }
                            else if (allLinks[h].end1 == -1) {
                                leftLink.end1 = 0;
                                leftLink.neighbor1 = allreads.size();
                                leftLink.end2 = allLinks[h].end2;
                                leftLink.neighbor2 = allLinks[h].neighbor2;
                                allreads[leftLink.neighbor2].add_link(allLinks.size(), allLinks[h].end2);
                            }
                            else {
                                // cout << "WHAAAT" << endl;
                            }
                            r.add_link(allLinks.size(), 0);
                            allLinks.push_back(leftLink);
                        }
                    }
                    
                    Link rightLink;
                    rightLink.CIGAR = "0M";
                    rightLink.end1 = 1;
                    rightLink.neighbor1 = allreads.size();
                    rightLink.end2 = -1;
                    rightLink.group = group.first;
                    allLinks.push_back(rightLink);
                    r.add_link(allLinks.size()-1, 1);
                    futureHangingLinks.push_back(allLinks.size()-1);

                    allreads.push_back(r);
                    backbones_reads.push_back(allreads.size()-1);
                    cout << "now creating the different contigs : " << r.name << endl;

                }
                hangingLinks = futureHangingLinks;
                n += 1;
            }

            //now wrap up the right of the contig
            int left = partitions[backbone][partitions[backbone].size()-1].first.second; //rightmost interval
            string right = allreads[backbone].sequence_.str().substr(left, allreads[backbone].sequence_.str().size()-left);
            Read r (right);
            r.name = allreads[backbone].name + "_"+ to_string(left)+ "_" + to_string(0);
            for (int h : hangingLinks){
                Link leftLink;
                leftLink.CIGAR = allLinks[h].CIGAR;
                if (allLinks[h].end2 == -1){
                    leftLink.end2 = 0;
                    leftLink.neighbor2 = allreads.size();
                    leftLink.end1 = allLinks[h].end1;
                    leftLink.neighbor1 = allLinks[h].neighbor1;
                    allreads[leftLink.neighbor1].add_link(allLinks.size(), allLinks[h].end1);
                }
                else if (allLinks[h].end1 == -1) {
                    leftLink.end1 = 0;
                    leftLink.neighbor1 = allreads.size();
                    leftLink.end2 = allLinks[h].end2;
                    leftLink.neighbor2 = allLinks[h].neighbor2;
                    allreads[leftLink.neighbor2].add_link(allLinks.size(), allLinks[h].end2);
                }
                else {
                    // cout << "WHAAAT" << endl;
                }
                r.add_link(allLinks.size(), 0);
                allLinks.push_back(leftLink);
                
            }

            for (int linkIdx : allreads[backbone].get_links_right()){
                if (allLinks[linkIdx].neighbor1 == backbone){
                    allLinks[linkIdx].end1 = 1;
                    allLinks[linkIdx].neighbor1 = allreads.size() ;
                }
                else {
                    allLinks[linkIdx].end2 = 1;
                    allLinks[linkIdx].neighbor2 = allreads.size() ;
                }
                r.add_link(linkIdx, 1);
            }

            allreads.push_back(r);
            backbones_reads.push_back(allreads.size()-1);
            cout << "now creating the different contigs : " << r.name << endl;

            allreads[backbone].name = "delete_me"; //output_gfa will understand that and delete the contig

        }
    }

    
}


unordered_map<int, set<int>> stitch(vector<int> &par, vector<int> &neighbor){

    unordered_map<int, unordered_map<int,int>> fit_left; //each parts maps to what left part ?
    unordered_map<int, int> cluster_size; 
    unordered_map<int,set<int>> stitch;

    for (auto r = 0 ; r < par.size() ; r++){
        if (par[r] != -1 && neighbor[r] != -1){
            if (fit_left.find(par[r]) != fit_left.end()){
                if (fit_left[par[r]].find(neighbor[r]) != fit_left[par[r]].end()){
                    fit_left[par[r]][neighbor[r]] += 1;
                }
                else{
                    fit_left[par[r]][neighbor[r]] = 1;
                }
                cluster_size[par[r]] += 1;
            }
            else{
                fit_left[par[r]][neighbor[r]] = 1;
                cluster_size[par[r]] = 1;
                stitch[par[r]] = {};
            }
        }
    }

    //now give all associations
    for (auto fit : fit_left){
        for (auto candidate : fit.second){
            if (candidate.second > 0.3*cluster_size[fit.first]){ //good compatibility
                stitch[fit.first].emplace(candidate.first);
            }
        }
    }

    return stitch;
}


