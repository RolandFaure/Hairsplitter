#ifndef GENERATE_MSA
#define GENERATE_MSA

#include "read.h"
#include "Partition.h"
#include <vector>

//small struct to return a slightly complicated result
struct distancePartition{
    int n00;
    int n01;
    int n10;
    int n11;
    float score;
    short phased; // worth -1 or 1
    bool augmented; //to know if the partition was augmented or not
    Column partition_to_augment;
};


void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, std::vector<std::vector<int>> &partitions, bool assemble_on_assembly);

float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<Column> &snps, int backboneReadIndex, Partition &truePar, bool assemble_on_assembly);
std::string consensus_reads(std::string &backbone, std::vector <std::string> &polishingReads);

std::vector<int>  separate_reads(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<Column> &snps, float minDistance, int numberOfReads);

distancePartition distance(Partition &par1, Column &par2);
distancePartition distance(Partition &par1, Partition &par2, int threshold_p);
float computeChiSquare(distancePartition dis);

std::vector<int> threadHaplotypes(std::vector<Partition> &listOfFinalPartitions, int numberOfReads);
std::vector<int> rescue_reads(std::vector<int> &threadedClusters, std::vector<Column> &snps, std::vector<size_t> &suspectPostitions);

#endif