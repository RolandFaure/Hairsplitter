#ifndef MODGFA
#define MODGFA

#include "check_overlaps.h"
#include <unordered_map>
#include <set>

void modify_GFA(std::string refFile, std::string readsFile, std::vector <Read> &allreads, std::vector<unsigned long int> &backbones_reads,  std::vector <Overlap> &allOverlaps, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::vector<int>> >> &partitions, std::vector<Link> &allLinks,
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits, int num_threads);

void modify_FASTA(std::string refFile, std::string readsFile, std::vector <Read> &allreads, std::vector<unsigned long int> &backbones_reads,  std::vector <Overlap> &allOverlaps, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::vector<int>> >> &partitions,
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits, int num_threads);

std::unordered_map<int, std::set<int>> stitch(std::vector<int> &par, std::vector<int> &neighbor);
std::unordered_map<int, double> recompute_depths(std::pair<std::pair<int,int>, std::vector<int>> &interval, std::vector<std::pair<int,int>>& readBorders, double originalDepth);

#endif