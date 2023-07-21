#ifndef CALL_VARIANTS_H
#define CALL_VARIANTS_H

#include <string>
#include <vector>

#include "read.h"
#include "robin_hood.h"
#include "tools.h"

float generate_msa(
    long int bbcontig, 
    std::vector <Overlap> &allOverlaps, 
    std::vector <Read> &allreads, 
    std::vector<Column> &snps, 
    robin_hood::unordered_map<int, int> &insertionPos,
    int backboneReadIndex, 
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits, 
    std::string &newref,
    std::string &tmpFolder,
    bool DEBUG);

void call_variants(
    std::vector<Column> &snps, 
    std::vector <Read> &allreads,
    std::vector <Overlap> &allOverlaps, 
    long int contig,
    std::string &ref,
    std::vector<size_t> &suspectPostitions,
    float &meanError, 
    std::string &tmpFolder,
    std::string &outputFile,
    bool DEBUG);

#endif