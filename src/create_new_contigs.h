#ifndef CREATE_NEW_CONTIGS_H
#define CREATE_NEW_CONTIGS_H


#include <vector>
#include <string>
#include <set>
#include <map>
#include <unordered_map>
#include "Partition.h"
#include "read.h"

void parse_split_file(
    std::string& file, 
    std::vector <Read> &allreads,
    std::vector <Overlap> &allOverlaps,
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::vector<int> > > > &partitions);

void modify_GFA(
    std::string readsFile, 
    std::vector <Read> &allreads, 
    std::vector<unsigned long int> &backbones_reads,
    std::vector <Overlap> &allOverlaps, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::vector<int> > > > &partitions,
    std::vector<Link> &allLinks,
    int num_threads,
    std::string &outFolder,
    float errorRate,
    bool polish,
    std::string &techno,
    std::string &MINIMAP, 
    std::string &RACON,
    bool DEBUG);

std::unordered_map<int, std::set<int>> stitch(std::vector<int> &par, std::vector<int> &neighbor, int position);
std::unordered_map<int, double> recompute_depths(std::pair<int,int> &limits, std::vector<int> &partition, double originalDepth);


void output_GAF(
    std::vector <Read> &allreads, 
    std::vector<unsigned long int> &backbone_reads, 
    std::vector<Link> &allLinks, 
    std::vector <Overlap> &allOverlaps, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::vector<int>  > >> &partitions,
    std::string outputGAF);

void merge_intervals(std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::vector<int>  > >> &partitions);

#endif



