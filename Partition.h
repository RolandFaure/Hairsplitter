#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
//class partition : a partition of reads, as a consensus of many partitions
class Partition {

public :
    Partition(int size); //size is the total number of reads in the pileup
    Partition(std::vector <char>& snp);
    void print(); //little function to print the consensus
    std::vector<short> getPartition();
    std::vector<float> getConfidence();
    std::vector<int> getMore();
    std::vector<int> getLess();
    void augmentPartition(std::vector<short> &newPar);
    void mergePartition(Partition p, short phased); 
    bool isInformative(float errorRate);
    int number();
    int size();

private :
    std::vector<short> mostFrequentBases; // at each position, 3 possibilities : 1 for allele1, -1 for allele2 and 0 for non-attributed-yet.
    //then the counts of each allele at each position
    std::vector<int> moreFrequence;
    std::vector<int> lessFrequence;
    int numberOfOccurences; //  Int is the number of partition supporting this
};

#endif