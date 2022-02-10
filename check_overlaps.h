#ifndef GENERATE_MSA
#define GENERATE_MSA

#include "read.h"
#include "Partition.h"
#include <vector>

//small struct to return a slightly complicated result
struct distancePartition{
    int nmatch;
    int nmismatch;
    int nonComparable;
    float chisquare; //chi square test on the similarity of the two partition
};



void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, std::vector<Partition> &partitions);

float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, int backboneReadIndex, std::vector<Partition> &partitions);
Partition separate_reads(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, float minDistance);

distancePartition distance(Partition &par1, std::vector<char> &par2, float errorRate);
bool distance(Partition &par1, Partition &par2, float thresholdChi, int threshold_p);

#endif