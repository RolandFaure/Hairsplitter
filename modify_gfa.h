#ifndef MODGFA
#define MODGFA

#include "check_overlaps.h"
#include <unordered_map>

void modify_GFA(std::string refFile, std::vector <Read> &allreads, std::vector<unsigned long int> &backbones_reads,  std::vector <Overlap> &allOverlaps, 
    std::unordered_map<unsigned long int ,std::vector<int>> &partitions, std::string outputFile, std::vector<Link> &allLinks);


#endif