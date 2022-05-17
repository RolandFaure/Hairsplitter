#include "modify_gfa.h"
#include <set>
#include<algorithm> //for "sort"

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
            std::unordered_map<unsigned long int ,std::vector<int>> &partitions, string outputFile, vector<Link> &allLinks){

    for (auto backbone : backbones_reads){

        unordered_map<int, pair<int, int>> intervals; //maps each set of reads to its corresponding interval
        unordered_map<int, vector<string>> readsPerPart; //list of all reads of each part

        if (partitions.find(backbone) != partitions.end()){ //means there is a partition there

            auto partition = partitions[backbone];

            for (auto n = 0 ; n < allreads[backbone].neighbors_.size() ; n++){
                
                if (partition[n] !=  0){

                    pair <int,int> interval;
                    Overlap overlap = allOverlaps[allreads[backbone].neighbors_[n]];
                    Read r;

                    if (overlap.sequence2 == backbone){
                        interval.first = overlap.position_2_1;
                        interval.second = overlap.position_2_2;
                        r = allreads[overlap.sequence1];
                    }
                    else if (overlap.sequence1 == backbone){
                        interval.first = overlap.position_1_1;
                        interval.second = overlap.position_1_2;
                        r = allreads[overlap.sequence2];
                    }

                    if (intervals.find(partition[n]) == intervals.end()){
                        intervals[partition[n]] = interval;
                        readsPerPart[partition[n]] = {r.sequence_.str()};
                    }
                    else {
                        readsPerPart[partition[n]].push_back(r.sequence_.str());
                        if (interval.first < intervals[partition[n]].first){intervals[partition[n]].first = interval.first;}
                        if (interval.second > intervals[partition[n]].second){intervals[partition[n]].second = interval.second;}
                    }
                }

            }

            //now modify the gfa if there are more than two parts to the partition
            if (intervals.size() > 1){

                //first define the intervals where the contig needs splitting
                vector <pair<pair<int,int>, set<int>>> splitIntervals;//associate all corresponding parts to the splitting interval

                for (auto interval : intervals){
                    bool found = false;
                    auto limits = interval.second;
                    int s = 0;
                    for (auto splitInterval : splitIntervals){
                        if (limits.first <= splitInterval.first.first && limits.second > splitInterval.first.first){
                            found = true;
                            splitIntervals[s].first.first = limits.first;
                            splitIntervals[s].first.second = max(limits.second, splitInterval.first.second);
                            splitIntervals[s].second.emplace(interval.first);
                            break;
                        }
                        else if (limits.first <= splitInterval.first.second && limits.second > splitInterval.first.second){
                            found = true;
                            splitIntervals[s].first.second = limits.second;
                            splitIntervals[s].second.emplace(interval.first);
                            break;
                        }
                        s++;
                    }
                    if (!found){
                        set<int> s = {interval.first};
                        splitIntervals.push_back(make_pair(limits, s));
                    }
                }

                //sort split intervals per leftmost position
                std::sort(splitIntervals.begin(), splitIntervals.end(), [](pair<pair<int,int>, set<int>>& lhs, pair<pair<int,int>, set<int>>& rhs) {
                    return lhs.first.first < rhs.first.first;
                });

                //now all intervals where contig need splitting are defined, let's create the new contigs ! 
                //first, create hangingLinks, a list of links left of what we'll be building
                vector<int> hangingLinks;
                for (int linkIdx : allreads[backbone].get_links_left()){
                    if (allLinks[linkIdx].neighbor1 == backbone){
                        allLinks[linkIdx].end1 = -1;
                    }
                    else {
                        allLinks[linkIdx].end2 = -1;
                    }
                    hangingLinks.push_back(linkIdx);
                }

                int lastRightmostPosition = 0;
                for (auto spl : splitIntervals){

                    if (spl.second.size() > 1) { //if there is only 1 group on this stretch, no need to re-polish it
                        string toPolish = allreads[backbone].sequence_.str().substr(spl.first.first, spl.first.second-spl.first.first);

                        //build the sequence left of here and attach it to hanging links
                        Read left(allreads[backbone].sequence_.str().substr(lastRightmostPosition, spl.first.first));
                        left.name = allreads[backbone].name + "_" +  to_string(lastRightmostPosition);

                        for(auto l : hangingLinks){
                            left.add_link(l, 0);
                            if (allLinks[l].end1 == -1){
                                allLinks[l].end1 = 0;
                                allLinks[l].neighbor1 = allreads.size();
                            }
                            else{
                                allLinks[l].end2 = 0;
                                allLinks[l].neighbor2 = allreads.size();
                            }
                        }
                        allreads.push_back(left);
                        auto leftIdx = allreads.size()-1;
                        backbones_reads.push_back(allreads.size()-1);
                        hangingLinks = {};

                        //now create all the new sequences
                        int g = 0;
                        for (auto group : spl.second){
                            string newcontig = consensus_reads(toPolish, readsPerPart[group]);
                            Read r(newcontig);
                            r.name = allreads[backbone].name + "_"+ to_string(spl.first.first)+ "_" + to_string(g);
                            g++;

                            Link leftLink;
                            leftLink.CIGAR = "0M";
                            leftLink.end1 = 0;
                            leftLink.neighbor1 = allreads.size();
                            leftLink.end2 = 1;
                            leftLink.neighbor2 = leftIdx;
                            allLinks.push_back(leftLink);
                            r.add_link(allLinks.size()-1, 0);
                            allreads[leftIdx].add_link(allLinks.size()-1, 1);
                            
                            Link rightLink;
                            rightLink.CIGAR = "0M";
                            rightLink.end1 = 1;
                            rightLink.neighbor1 = allreads.size();
                            rightLink.end2 = -1;
                            allLinks.push_back(rightLink);
                            r.add_link(allLinks.size()-1, 1);
                            hangingLinks.push_back(allLinks.size()-1);

                            allreads.push_back(r);
                            backbones_reads.push_back(allreads.size()-1);
                            cout << "now creating the different contigs : " << r.name << endl;

                        }
                        lastRightmostPosition = spl.first.second;
                    }
                }

                //now wrap up the contig, link it with contig at the right
                Read last(allreads[backbone].sequence_.str().substr(lastRightmostPosition, allreads[backbone].size()-lastRightmostPosition));
                last.name = allreads[backbone].name + "_" + to_string(lastRightmostPosition);

                //links left
                for(auto l : hangingLinks){
                    last.add_link(l, 0);
                    if (allLinks[l].end1 == -1){
                        allLinks[l].end1 = 0;
                        allLinks[l].neighbor1 = allreads.size();
                    }
                    else{
                        allLinks[l].end2 = 0;
                        allLinks[l].neighbor2 = allreads.size();
                    }
                }
                allreads.push_back(last);
                backbones_reads.push_back(allreads.size()-1);

                //links right
                for(auto l : allreads[backbone].get_links_right()){
                    if (allLinks[l].neighbor1 == backbone){
                        allLinks[l].neighbor1 = allreads.size()-1;
                        allLinks[l].end1 = 1;
                    }
                    else {
                        allLinks[l].neighbor2 = allreads.size()-1;
                        allLinks[l].end2 = 1;
                    }
                }                
            }
        }

    }
    
}
