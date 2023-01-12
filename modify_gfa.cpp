#include "modify_gfa.h"
#include<algorithm> //for "sort"
#include <fstream>
#include <omp.h>
#include "input_output.h"

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

extern bool DEBUG;

/**
 * @brief Modify the input GFA according to the way the reads have been split.
 * 
 * @param readsFile File containing all the reads
 * @param allreads vector containing all the reads (without their actual sequence)
 * @param backbones_reads vector listing all the backbone reads
 * @param allOverlaps vector containing all the overlaps
 * @param partitions for each backbone, contains a vector of the intervals, and for each interval, the partition of the reads
 * @param allLinks vector containing all the links of the GFA file (new one will be added)
 * @param readLimits for each backbone, contains the limits (in term of coordinates) of all its neighbors on the backbone
 * @param num_threads number of threads to use
 */
void modify_GFA(
    std::string readsFile, 
    vector <Read> &allreads, 
    vector<unsigned long int> &backbones_reads, 
    vector <Overlap> &allOverlaps,
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions,
    vector<Link> &allLinks,
    unordered_map <int, vector<pair<int,int>>> &readLimits, 
    int num_threads)
    {

    int max_backbone = backbones_reads.size(); //fix that because backbones will be added to the list but not separated 
    string log_text = ""; //text that will be printed out in the output.txt

    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int b = 0 ; b < max_backbone ; b++){

        //first load all the reads
        parse_reads_on_contig(readsFile, backbones_reads[b], allOverlaps, allreads);

        //then separate the contigs

        string local_log_text = ""; //text that will be printed out in the output.txt generated in this thread
        if (DEBUG){
            #pragma omp critical
            {
                cout << "Thread " << omp_get_thread_num() << " looking at " << allreads[backbones_reads[b]].name << ", can it be found ? " <<
                     (partitions.find(backbones_reads[b]) != partitions.end()) << endl;
                if (partitions.find(backbones_reads[b]) != partitions.end()){
                    cout << "If ueie so, how many intervals ? " << partitions[backbones_reads[b]].size() << endl;
                }
            }
        }

        local_log_text += "---- contig: " + allreads[backbones_reads[b]].name + " ----\n\n";

        string thread_id = std::to_string(omp_get_thread_num());
        int backbone = backbones_reads[b];

        //if partitions[backbone].size() == 0, see if we need to repolish or not, depending on wether the coverage is coherent or not
        bool dont_recompute_contig = false;
        if (partitions[backbone].size() == 0 && allreads[backbones_reads[b]].depth > 1){
            double total_depth = 0;
            for (auto n : allreads[backbones_reads[b]].neighbors_){
                total_depth += allOverlaps[n].position_1_2 - allOverlaps[n].position_1_1;
            }
            double new_depth = total_depth / allreads[backbones_reads[b]].sequence_.size();

            if (new_depth / allreads[backbones_reads[b]].depth > 0.7){
                dont_recompute_contig = true;
            }
        }
        if (partitions.find(backbone) != partitions.end() && partitions[backbone].size() > 0 && !dont_recompute_contig){
            //stitch all intervals of each backbone read
            vector<unordered_map <int,set<int>>> stitches(partitions[backbone].size()); //aggregating the information from stitchesLeft and stitchesRight to know what link to keep

            for (int n = 0 ; n < partitions[backbone].size() ; n++){
                //for each interval, go through the different parts and see with what part before and after they fit best
                if (n > 0){
                    auto stitchLeft = stitch(partitions[backbone][n].second.first, partitions[backbone][n-1].second.first);
                    auto stitchRight = stitch(partitions[backbone][n-1].second.first, partitions[backbone][n].second.first);

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
                if (allLinks[linkIdx].neighbor1 == backbone && allLinks[linkIdx].end1 == 0){
                    allLinks[linkIdx].end1 = -1;
                }
                else {
                    allLinks[linkIdx].end2 = -1;
                }
                allLinks[linkIdx].group = 0;
                hangingLinks.push_back(linkIdx);
            }

            //construct singlepolish, the sequence polished by all the reads.
            //it may be different from the polished version we already have because it has been polished using reads that maybe align elsewhere: compare the expected depth with number of reads aligning
            //compute the depth from the number of aligning reads
            std::pair<int,int> limitsAll= std::make_pair(0, allreads[backbone].sequence_.size()-1);
            vector<int> partition1 (allreads[backbone].neighbors_.size(), 1);
            unordered_map <int, double> newdepths = recompute_depths(limitsAll , partition1, readLimits[backbone], allreads[backbone].depth);

            string singlepolish = allreads[backbone].sequence_.str();
            if (newdepths[1]<0.8*allreads[backbone].depth || newdepths[1]>1.3*allreads[backbone].depth || allreads[backbone].depth == 0){ //if the depths do not match
                vector<string> allneighbors;
                for (int n = 0 ; n < allreads[backbone].neighbors_.size() ; n++){
                    auto idxRead = allOverlaps[allreads[backbone].neighbors_[n]].sequence1;
                    if (allOverlaps[allreads[backbone].neighbors_[n]].strand){
                        allneighbors.push_back(allreads[idxRead].sequence_.str());
                    }
                    else{
                        allneighbors.push_back(allreads[idxRead].sequence_.reverse_complement().str());
                    }
                }
                string seqbackbone = allreads[backbone].sequence_.str();
                singlepolish = consensus_reads(seqbackbone, allneighbors, thread_id);
            }

            int n = 0;
            for (auto interval : partitions[backbone]){

                unordered_map<int, vector<string>> readsPerPart; //list of all reads of each part
                if (omp_get_thread_num() == 0 && DEBUG){
                    cout << "in interval " << interval.first.first << " <-> " << interval.first.second << endl;
                }
                local_log_text += " - Between positions " + to_string(interval.first.first) + " and " + to_string(interval.first.second) + " of the contig, I've created these contigs:\n";

                //taking exactly the right portion of read we need
                for (int r = 0 ; r < interval.second.first.size(); r++){
                    if (interval.second.first[r] != -1){
                        int clust = interval.second.first[r];
                        int limitLeft = allOverlaps[allreads[backbone].neighbors_[r]].position_1_1; //the limit of the read that we should use
                        int limitRight = allOverlaps[allreads[backbone].neighbors_[r]].position_1_2;
                        auto idxRead = allOverlaps[allreads[backbone].neighbors_[r]].sequence1;

                        string clippedRead; //the read we're aligning with good orientation and only the part we're interested in

                        if (allOverlaps[allreads[backbone].neighbors_[r]].strand){
                            clippedRead = allreads[idxRead].sequence_.subseq(limitLeft, limitRight-limitLeft+1).str();
                        }
                        else{
                            clippedRead = allreads[idxRead].sequence_.subseq(limitLeft, limitRight-limitLeft+1).reverse_complement().str();
                        }

                        if (readsPerPart.find(clust) == readsPerPart.end()){
                            readsPerPart[clust] = {clippedRead};
                        }
                        else {
                            readsPerPart[clust].push_back(clippedRead);
                        }
                        // cout << "Read " << allreads[idxRead].name << " is in cluster " << clust << endl;
                    }
                }
                // cout << endl;

                vector<int> futureHangingLinks;

                cout << "oueopua recomputing depths" << endl;
                unordered_map <int, double> newdepths = recompute_depths(interval.first, interval.second.first, readLimits[backbone], allreads[backbone].depth);
                cout << "recomputed poaiepae" << endl;

                for (auto group : readsPerPart){
                    
                    string toPolish = interval.second.second[group.first];

                    string newcontig = "";
                    if (readsPerPart.size() > 1){

                        // newcontig = local_assembly(group.second);

                        if (newcontig == ""){//if the assembly was not successful for one reason or another
                            //toPolish2 should be polished with a little margin on both sides to get cleanly first and last base
                            string toPolish2 = allreads[backbone].sequence_.str().substr(max(0, interval.first.first - 100), min(interval.first.first, 100)) + toPolish + allreads[backbone].sequence_.str().substr(interval.first.second+1, min(100, int(allreads[backbone].sequence_.str().size())-interval.first.second-1));

                            newcontig = consensus_reads(toPolish2, group.second, thread_id);
                        }
                        EdlibAlignResult result = edlibAlign(toPolish.c_str(), toPolish.size(),
                                    newcontig.c_str(), newcontig.size(),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

                        newcontig = newcontig.substr(max(0,result.startLocations[0]), min(result.endLocations[0]-result.startLocations[0]+2, int(newcontig.size())-result.startLocations[0]));

                        edlibFreeAlignResult(result);
                        
                        
                    }
                    else {
                        string extract = allreads[backbone].sequence_.str().substr(interval.first.first, interval.first.second-interval.first.first+1);
                        EdlibAlignResult result = edlibAlign(extract.c_str(), extract.size(),
                                    singlepolish.c_str(), singlepolish.size(),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                        newcontig = singlepolish.substr(result.startLocations[0], result.endLocations[0]-result.startLocations[0]);

                    }

                    Read r(newcontig);
                    r.name = allreads[backbone].name + "_"+ to_string(interval.first.first)+ "_" + to_string(group.first);
                    r.depth = newdepths[group.first];

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
                    if (omp_get_thread_num() == 0 && DEBUG){
                        cout << "created the contig " << r.name << endl;
                    }
                    local_log_text += "   - " + r.name + "\n";
                }
                hangingLinks = futureHangingLinks;
                n += 1;
            }
            //now wrap up the right of the contig
            int left = partitions[backbone][partitions[backbone].size()-1].first.second+1; //rightmost interval
            string right = allreads[backbone].sequence_.str().substr(left, allreads[backbone].sequence_.size()-left);
            std::pair<int,int> limits= std::make_pair(left, allreads[backbone].sequence_.size()-1);
            newdepths = recompute_depths(limits , partition1, readLimits[backbone], allreads[backbone].depth);

            string contig;
            if (right.size() > 0){
                EdlibAlignResult result = edlibAlign(right.c_str(), right.size(),
                    singlepolish.c_str(), singlepolish.size(),edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                contig = singlepolish.substr(result.startLocations[0]-1, result.endLocations[0]-result.startLocations[0]+2);
                edlibFreeAlignResult(result);
                }
            else{
                contig = "";
            }
            Read r (contig);
            r.name = allreads[backbone].name + "_"+ to_string(left)+ "_" + to_string(0);
            r.depth = newdepths[1];

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

            //now re-build the links right of backbone
            for (int linkIdx : allreads[backbone].get_links_right()){
                if (allLinks[linkIdx].neighbor1 == backbone && allLinks[linkIdx].end1 == 1){
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
            if (omp_get_thread_num() == 0 && DEBUG){
                cout << "now creating the different contigs : " << r.name << endl;
            }
            local_log_text += " - Between positions " + to_string(left) + " and " + to_string(allreads[backbone].sequence_.size()) + " of the contig, I've created these contigs:\n";
            local_log_text +=  "   - " + r.name + "\n\n";

            allreads[backbone].name = "delete_me"; //output_gfa will understand that and delete the contig

        }
        else{
            local_log_text += "Nothing to do\n\n";
        }
        #pragma omp critical
        {
            log_text += local_log_text;
        }

        //free up memory by deleting the sequence of the reads used there
        for (auto n : allreads[backbones_reads[b]].neighbors_){
            if (allOverlaps[n].sequence1 != backbones_reads[b]){
                allreads[allOverlaps[n].sequence1].delete_sequence();
            }
            else{
                allreads[allOverlaps[n].sequence2].delete_sequence();
            }
        }
    }
    std::ofstream o("output.txt");
    o << log_text << endl;
}

//input : list of links from first par to second par and reciproqually
//output : for each part of par, the set of parts of neighbor to which it should be linked
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

//input : an interval, the list of the limits of the reads on the backbone, the depth of the contig of origin
//output : the recomputed read coverage for each of the new contigs, (scaled so that the total is the original depth)
std::unordered_map<int, double> recompute_depths(std::pair<int,int> &limits, std::vector<int> &partition, std::vector<std::pair<int,int>>& readBorders, double originalDepth){

    unordered_map <int, double> newCoverage;
    int lengthOfInterval = limits.second-limits.first+1; //+1 to make sure we do not divide by 0

    for (auto c = 0 ; c < partition.size() ; c++){

        if (newCoverage.find(partition[c]) == newCoverage.end()){
            newCoverage[partition[c]] = 0;
        }

        newCoverage[partition[c]] += max(0.0, double(min(limits.second, readBorders[c].second)-max(limits.first, readBorders[c].first))/lengthOfInterval );

    }

    //now scale all the coverages to obtain exactly the original coverage
    if (originalDepth != -1){ //that would mean we do not know anything about the original depth

        double totalCoverage = 0;
        for (auto cov : newCoverage){
            if (cov.first != -1){
                totalCoverage += cov.second;
            }
        }
        
        if (totalCoverage != 0){
            for (auto cov : newCoverage){
                // newCoverage[cov.first] = cov.second * originalDepth/totalCoverage; DEBUG
            }
        }

    }

    return newCoverage;

}

