#include "modify_gfa.h"
#include <algorithm> //for "sort"
#include <fstream>
#include <omp.h>
#include "input_output.h"
#include "tools.h"

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
        #pragma omp critical
        {
            parse_reads_on_contig(readsFile, backbones_reads[b], allOverlaps, allreads);
        }

        //then separate the contigs

        string local_log_text = ""; //text that will be printed out in the output.txt generated in this thread
        if (DEBUG){
            #pragma omp critical
            {
                cout << "Thread " << omp_get_thread_num() << " looking at " << allreads[backbones_reads[b]].name  << endl;
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
            vector<unordered_map <int,set<int>>> stitches(partitions[backbone].size()); //aggregates the information from stitchesLeft and stitchesRight to know what link to keep

            for (int n = 0 ; n < partitions[backbone].size() ; n++){
                //for each interval, go through the different parts and see with what part before and after they fit best
                if (n > 0){
                    std::unordered_map<int, std::set<int>> stitchLeft = stitch(partitions[backbone][n].second.first, 
                                                                            partitions[backbone][n-1].second.first, 
                                                                            partitions[backbone][n].first.first, 
                                                                            readLimits[backbone]);
                    // std::unordered_map<int, std::set<int>> stitchRight = stitch(
                    //                                                         partitions[backbone][n-1].second.first,
                    //                                                         partitions[backbone][n].second.first,
                    //                                                         partitions[backbone][n].first.first, 
                    //                                                         readLimits[backbone]);


                    for (auto s : stitchLeft){
                        stitches[n][s.first] = s.second;
                    }

                    // for (auto s : stitchRight){
                    //     for (int neighbor : s.second){
                    //         stitches[n][neighbor].emplace(s.first);
                    //     }
                    // }

                    //make sure all contigs on the left and right are stitched
                    set<int> all_contigs_left;
                    for (int a : partitions[backbone][n-1].second.first){
                        all_contigs_left.emplace(a);
                    }
                    all_contigs_left.erase(-1);
                    all_contigs_left.erase(-2);

                    set<int> all_contigs_right;
                    for (int a : partitions[backbone][n].second.first){
                        all_contigs_right.emplace(a);
                    }
                    all_contigs_right.erase(-1);
                    all_contigs_right.erase(-2);

                    for (auto s : stitchLeft){
                        if (s.second.size() == 0){
                            stitchLeft[s.first] = all_contigs_left;
                        }
                    }

                    //check if all contigs on the left of the junction are stitched
                    set <int> stitchedContigs;
                    for (auto s : stitches[n]){
                        for (int neighbor : s.second){
                            stitchedContigs.emplace(neighbor);
                        }
                    }
                    for (auto contig : all_contigs_left){
                        if (stitchedContigs.find(contig) == stitchedContigs.end()){
                            for (auto s : stitchLeft){
                                stitches[n][s.first].emplace(contig);
                            }
                        }
                    }
                    // cout << "stitched contigs : "; for(auto i : stitchedContigs) {cout <<i << ";";} cout << endl;
                    // for (auto s : stitchRight){
                    //     if (stitchedContigs.find(s.first) == stitchedContigs.end()){
                    //         // cout << "gluing back" << endl;
                    //         for (auto s2 : stitchLeft){
                    //             stitches[n][s2.first].emplace(s.first);
                    //         }
                    //     }
                    // }

                    // if (partitions[backbone][n].first.first > 41100 && partitions[backbone][n].first.first < 42900 && omp_get_thread_num() == 0){
                    //     cout << "stitching : ";
                    //     for (auto s: stitches[n]){
                    //         cout << s.first << "<->" ;
                    //         for (auto s2 : s.second) {cout << s2 << ",";}
                    //         cout << "  ;  ";
                    //     } 
                    //     cout << endl;

                    //     cout << "here are partitions[backbone][n-1].second.first and partitions[backbone][n].second.first : " << endl;
                    //     for (auto a : partitions[backbone][n-1].second.first){
                    //         cout << a+1 << " ";
                    //     }
                    //     cout << endl;
                    //     for (auto a : partitions[backbone][n].second.first){
                    //         cout << a+1 << " ";
                    //     }
                    //     cout << endl;
                    //     while (true)
                    //     {
                    //         /* code */
                    //     }
                        
                    // }
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

            //compute the depth from the number of aligning reads
            std::pair<int,int> limitsAll= std::make_pair(0, allreads[backbone].sequence_.size()-1);
            vector<int> partition1 (allreads[backbone].neighbors_.size(), 1);
            unordered_map <int, double> newdepths = recompute_depths(limitsAll , partition1, readLimits[backbone], allreads[backbone].depth);


            int n = 0;
            for (auto interval : partitions[backbone]){

                unordered_map<int, vector<string>> readsPerPart; //list of all reads of each part
                if (omp_get_thread_num() == 0 && DEBUG){
                    cout << "in interval " << interval.first.first << " <-> " << interval.first.second << endl;
                }
                local_log_text += " - Between positions " + to_string(interval.first.first) + " and " + to_string(interval.first.second) + " of the contig, I've created these contigs:\n";

                //taking exactly the right portion of read we need
                for (int r = 0 ; r < interval.second.first.size(); r++){
                    if (interval.second.first[r] > -1){
                        auto idxRead = allOverlaps[allreads[backbone].neighbors_[r]].sequence1;

                        int clust = interval.second.first[r];
                        int limitLeft = max(0,allOverlaps[allreads[backbone].neighbors_[r]].position_1_1-20); //the limit of the read that we should use with a little margin for a clean polish
                        int limitRight = min(allOverlaps[allreads[backbone].neighbors_[r]].position_1_2+20, int(allreads[idxRead].sequence_.size()));

                        string clippedRead; //the read we're aligning with good orientation and only the part we're interested in
                        
                        //cout << "cliipplling read: " << allreads[idxRead].name << " " << limitLeft << " " << limitRight << " " << allreads[idxRead].sequence_.size() << endl;

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
                    }
                }
                if (readsPerPart.size() == 0 && interval.second.first.size() > 0 && interval.second.first[0] <= -1){
                    readsPerPart[-1] = {}; //so that it defaults back to the consensus
                }
                // cout << endl;

                vector<int> futureHangingLinks;

                unordered_map <int, double> newdepths = recompute_depths(interval.first, interval.second.first, readLimits[backbone], allreads[backbone].depth);

                for (auto group : readsPerPart){
                    int overhang = 50; //margin we're taking at the ends of the contig to be sure to take exactly the right portion 
                    
                    int overhangLeft = min(interval.first.first, overhang);
                    int overhangRight = min(int(allreads[backbone].sequence_.size())-interval.first.second-1, overhang);
                    //toPolish should be polished with a little margin on both sides to get cleanly first and last base
                    string toPolish = allreads[backbone].sequence_.str().substr(max(0, interval.first.first - overhang), overhang) 
                        + interval.second.second[group.first].substr(max(0, overhang-overhangLeft), interval.second.second[group.first].size()+overhangLeft+overhangRight-2*overhang)
                        + allreads[backbone].sequence_.str().substr(interval.first.second+1+overhangRight-overhang, overhang);

                    // cout << "toPolisssdcvh: " << toPolish << endl;

                    string newcontig = "";
                    if (readsPerPart.size() > 1){

                        //because racon is not very good at polishing the first and last bases of contig, take the consensus as computed before
                        // if (interval.first.first <= 50){
                        //     //remove the first 50 bases of toPolish
                        //     toPolish = toPolish.substr(50, toPolish.size()-50);
                        // }
                        // if (interval.first.second >= allreads[backbone].sequence_.size()-50){
                        //     //remove the last 50 bases of toPolish
                        //     toPolish = toPolish.substr(0, toPolish.size()-50);
                        // }

                        if (newcontig == ""){
                            //if the contig is close to one end, tell racon to not polish the first or last bases
                            newcontig = consensus_reads(toPolish, group.second, overhang-overhangLeft, overhang-overhangRight, thread_id);
                        }

                        EdlibAlignResult result = edlibAlign(toPolish.c_str(), toPolish.size(),
                                    newcontig.c_str(), newcontig.size(),
                                    edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

                        string cig = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
                        string cigar = convert_cigar(cig);
                        // extract the part of newcontig that does not align to the first and last overhang bases of toPolish2
                        int posOnToPolish = 0;
                        int posOnNewContig = result.startLocations[0];
                        int posStartOnNewContig = 0;
                        int posEndOnNewContig = 0;
                        for (char c : cigar){
                            if (c == 'M'){
                                posOnToPolish++;
                                posOnNewContig++;
                            }
                            else if (c == 'D'){
                                posOnNewContig++;
                            }
                            else if (c == 'I'){
                                posOnToPolish++;
                            }
                            if (posOnToPolish == overhangLeft){
                                posStartOnNewContig = posOnNewContig;
                            }
                            if (posOnToPolish == toPolish.size()-overhangRight-1){
                                posEndOnNewContig = posOnNewContig;
                            }
                            // cout << "indices: " << posOnToPolish << " " << posOnNewContig << endl;
                        }
                        
                        newcontig = newcontig.substr(posStartOnNewContig, min(posEndOnNewContig-posStartOnNewContig+1, int(newcontig.size())-posStartOnNewContig));

                        //because racon is not very good at polishing the first and last bases of contig, take the consensus as computed before
                        // if (interval.first.first <= 50){
                        //     //add back the first 50 bases of newcontig
                        //     newcontig = interval.second.second[group.first].substr(0, 50) + newcontig;
                        // }
                        // if (interval.first.second >= allreads[backbone].sequence_.size()-50){
                        //     //add back the last 50 bases of newcontig
                        //     newcontig = newcontig + interval.second.second[group.first].substr(interval.second.second[group.first].size()-50, 50);
                        // }

                        edlibFreeAlignResult(result);                        
                    }
                    else {
                        newcontig = interval.second.second[group.first];
                    }

                    Read r(newcontig);
                    r.name = allreads[backbone].name + "_"+ to_string(interval.first.first)+ "_" + to_string(group.first);
                    r.depth = newdepths[group.first];

                    //now create all the links IF they are compatible with "stitches"  
                    set<int> linksToKeep;

                    if (interval.first.first == 6000 && group.first == 4){
                        cout << "stiretches: " << endl;
                        for (auto s : stitches[n][group.first]){
                            cout << s << " ";
                        }
                        cout << endl;
                    }
        
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
                    #pragma omp critical
                    {
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
                }
                hangingLinks = futureHangingLinks;
                n += 1;
            }
            //now wrap up the right of the contig
            int left = partitions[backbone][partitions[backbone].size()-1].first.second+1; //rightmost interval
            string right = allreads[backbone].sequence_.str().substr(left, allreads[backbone].sequence_.size()-left);
            std::pair<int,int> limits= std::make_pair(left, allreads[backbone].sequence_.size()-1);
            string contig = right;
            
            Read r (contig);
            r.name = allreads[backbone].name + "_"+ to_string(left)+ "_" + to_string(0);
            r.depth = newdepths[1];
            
            #pragma omp critical
            {
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

            // if (allreads[allreads.size()-1].name.substr(0,9) == "edge_13@0"){
            //     cout << "qfdklmdjlccjj " << endl;
            //     exit(1);
            // }

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
                allreads[allOverlaps[n].sequence1].free_sequence();
            }
            else{
                allreads[allOverlaps[n].sequence2].free_sequence();
            }
        }
    }
    std::ofstream o("output.txt");
    o << log_text << endl;
}

/**
 * @brief Tells to which parts of neighbor each part of par should be linked
 * 
 * @param par partitions on the partition
 * @param neighbor partitions on the neighbor
 * @param position position of the stitch in the backbone
 * @param readLimits limits of the reads in the backbone (so that reads that are not on the position of the stitch are not considered)
 * @return unordered_map<int, set<int>> map associating a set of partitions of neighbor matching each partition of par
 */
unordered_map<int, set<int>> stitch(vector<int> &par, vector<int> &neighbor, int position, vector<pair<int,int>> &readLimits){

    unordered_map<int, unordered_map<int,int>> fit_left; //each parts maps to what left part ?
    unordered_map<int, unordered_map<int,int>> fit_right; //each parts maps to what right part ?
    unordered_map<int, int> cluster_size; 
    unordered_map<int,set<int>> stitch;

    for (auto r = 0 ; r < par.size() ; r++){
        if (par[r] > -1 && neighbor[r] > -1 && readLimits[r].first <= position && readLimits[r].second >= position){
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

            if (fit_right.find(neighbor[r]) != fit_right.end()){
                if (fit_right[neighbor[r]].find(par[r]) != fit_right[neighbor[r]].end()){
                    fit_right[neighbor[r]][par[r]] += 1;
                }
                else{
                    fit_right[neighbor[r]][par[r]] = 1;
                }
            }
            else{
                fit_right[neighbor[r]][par[r]] = 1;
            }
        }
    }

    //now give all associations
    for (auto fit : fit_left){
        // find the best fit
        int best_fit = 0;
        for (auto candidate : fit.second){
            if (candidate.second > best_fit){
                best_fit = candidate.second;
            }
        }
        for (auto candidate : fit.second){
            if (candidate.second == best_fit){ //good compatibility
                stitch[fit.first].emplace(candidate.first);
            }
        }
    }

    for (auto fit : fit_right){
        //find the best fit
        int best_fit = 0;
        for (auto candidate : fit.second){
            if (candidate.second > best_fit){
                best_fit = candidate.second;
            }
        }
        for (auto candidate : fit.second){
            if (candidate.second == best_fit){
                stitch[candidate.first].emplace(fit.first);
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

