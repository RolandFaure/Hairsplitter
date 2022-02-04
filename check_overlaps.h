#ifndef GENERATE_MSA
#define GENERATE_MSA

#include "read.h"
#include <vector>

//small struct to return a slightly complicated result
struct distancePartition{
    int nmatch;
    int nmismatch;
    int nonComparable;
    float easyMatches; //number of matches you would expect by chance given the nucleotide composition of the two partitions
};

//class partition : a partition of reads, as a consensus of many partitions
class Partition {

public :
    Partition(int size); //size is the total number of reads in the pileup
    Partition(std::vector <char>& snp);
    void print(); //little function to print the consensus
    std::vector<short> getPartition();
    void augmentPartition(std::vector<short> &newPar); 

private :
    std::vector<short> mostFrequentBases; // at each position, 3 possibilities : 1 for allele1, -1 for allele2 and 0 for non-attributed-yet. Int is the number of partition supporting this
    //then the counts of each allele at each position
    std::vector<int> moreFrequence;
    std::vector<int> lessFrequence;
};

void checkOverlaps(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads);

float generate_msa(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps);
void separate_reads(long int read, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<std::vector<char>> &snps, float minDistance);

bool distance(Partition &par1, std::vector<char> &par2, float errorRate, std::vector<std::pair<Partition, int>> &allPartitions);

#endif