#ifndef PARTITION_H
#define PARTITION_H

#include <vector>
#include <list>

//structure of a column of an MSA
struct Column{
    std::vector<int> readIdxs;
    std::vector<char> content;
};

//structure for the partition in an arbitrary number of clusters
struct Separation{
    std::vector<int> partition;
    int left;
    int right;
};

//class partition : a partition of reads, as a consensus of many partitions
/**
 * @brief Class discribing an aggregation of several similar-looking columns of the MSA
 * 
 */
class Partition {

public :
    Partition(); //size is the total number of reads in the pileup
    Partition(Column &snp, int pos); //to initialize a binary partition
    void print(); //little function to print the consensus
    std::vector<int> getReads();
    std::vector<short> getPartition();
    std::vector<float> getConfidence();
    std::vector<int> getMore();
    std::vector<int> getLess();
    void augmentPartition(Column &newPar, int pos);
    void mergePartition(Partition &p, short phased); 
    void mergePartition(Partition &p); 
    bool isInformative(float errorRate, bool lastReadBiased);
    int number();
    float proportionOf1();

    float compute_conf();
    float get_conf();

    int get_left();
    int get_right();

    void flipPartition(); //transforms all 1s in 0s and vice-versa

private :
    std::vector <int> readIdx; //the list of reads is sparse : here are the filled indices (ordered)
    std::vector<short> mostFrequentBases; // at each position, 3 possibilities : 1 for allele1, -1 for allele2 and 0 for non-attributed-yet
    //then the counts of each allele at each position
    std::vector<int> moreFrequence;
    std::vector<int> lessFrequence;

    int numberOfOccurences; //  Int is the number of partition supporting this
    float conf_score;

    int pos_left; //position of leftmost position of the partition on the consensus
    int pos_right; //position of the rightmost position of the partition on the consensus
};

#endif