#ifndef GENERATE_MSA
#define GENERATE_MSA

#include "read.h"
#include "Partition.h"
#include <vector>

//small struct to return a slightly complicated result
struct distancePartition{
    int nmatch;
    int nmismatch;
    int n00;
    int n01;
    int n10;
    int n11;
    int nonComparable;
    float chisquare; //chi square test on the similarity of the two partition
    short phased; // worth -1 or 1
};



void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, std::vector<std::vector<short>> &partitions);

float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, int backboneReadIndex, Partition &truePar);

std::vector<short>  separate_reads(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, float minDistance);

distancePartition distance(Partition &par1, std::vector<char> &par2, float errorRate);
bool distance(Partition &par1, Partition &par2, float thresholdChi, int threshold_p);

std::vector<short> threadHaplotypes(std::vector<Partition> &listOfFinalPartitions);

#endif