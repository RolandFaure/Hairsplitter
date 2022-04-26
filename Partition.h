#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <list>

//structure of a column of an MSA
struct Column{
    std::list<int> readIdxs;
    std::vector<char> content;
};

//class partition : a partition of reads, as a consensus of many partitions
class Partition {

public :
    Partition(int size); //size is the total number of reads in the pileup
    Partition(Column &snp); //to initialize a binary partition
    void print(); //little function to print the consensus
    std::vector<int> getReads();
    std::vector<short> getPartition();
    std::vector<float> getConfidence();
    std::vector<int> getMore();
    std::vector<int> getLess();
    void augmentPartition(Column &newPar);
    void mergePartition(Partition &p, short phased); 
    void mergePartition(Partition &p); 
    bool isInformative(float errorRate, bool lastReadBiased);
    int number();
    float proportionOf1();

private :
    std::vector <int> readIdx; //the list of reads is sparse : here are the filled indices (ordered)
    std::vector<short> mostFrequentBases; // at each position, 3 possibilities : 1 for allele1, -1 for allele2 and 0 for non-attributed-yet
    //then the counts of each allele at each position
    std::vector<int> moreFrequence;
    std::vector<int> lessFrequence;
    int numberOfOccurences; //  Int is the number of partition supporting this
};

#endif