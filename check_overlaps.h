#ifndef GENERATE_MSA
#define GENERATE_MSA

#include "read.h"
#include "Partition.h"
#include "edlib.h"
#include "robin_hood.h"

#include <bindings/cpp/WFAligner.hpp>


#include <vector>
#include <unordered_map>
#include <thread>
#include <string>

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


void checkOverlaps(std::string fileReads, std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, 
            std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions, bool assemble_on_assembly,
            std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits,
            bool polish, int num_threads);

void compute_partition_on_this_contig(std::string fileReads, long int contig,std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::vector<unsigned long int> &backbones_reads, 
            std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions, bool assemble_on_assembly,
            std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits,
            bool polish);

float generate_msa(long int bbcontig, std::vector <Overlap> &allOverlaps, std::vector <Read> &allreads, std::vector<Column> &snps, robin_hood::unordered_map<int, int> &insertionPos,
    int backboneReadIndex, std::string &truePar, bool assemble_on_assembly, 
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits, std::vector<bool>& misalignedReads, bool polish);

std::string consensus_reads(std::string &backbone, std::vector <std::string> &polishingReads, std::string &id);
void compute_consensus_in_partitions(long int contig, std::vector<std::pair<std::pair<int,int>, std::vector<int>> > &partition, std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, 
    std::vector<Column> &snps, robin_hood::unordered_map<int, int> &insertionPositions,
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions);
// std::string local_assembly(std::vector <std::string> &reads);

std::vector< std::pair<std::pair<int,int>, std::vector<int>> > separate_reads(std::string& ref, std::vector<Column> &snps, float minDistance, int numberOfReads);

distancePartition distance(Partition &par1, Column &par2);
distancePartition distance(Partition &par1, Partition &par2, int threshold_p);
float computeChiSquare(distancePartition dis);

void clean_partition(long int backbone, Partition &originalPartition, std::vector <Read> &allreads,std::vector <Overlap> &allOverlaps);

std::vector<Partition> select_compatible_partitions(std::vector<Partition> &partitions, int numberOfReads, float errorRate);
std::vector<Partition> select_confident_partitions(std::vector<Partition> &partitions, std::vector<bool> trimmedListOfFinalPartitionBool, int numberOfReads, float errorRate);

std::vector< std::pair<std::pair<int,int>, std::vector<int>> > threadHaplotypes(std::vector<Partition> &compatiblePartitions, int numberOfReads,
    std::unordered_map <int, std::pair<int,int>> &clusterLimits);
int compatible_partitions(Partition &p1 , Partition &p2);
std::vector<int> threadHaplotypes_in_interval(std::vector<Partition> &listOfFinalPartitions, int numberOfReads);
bool extend_with_partition_if_compatible(std::vector<int> &alreadyThreadedHaplotypes, Partition &extension, int partitionIndex,
        std::unordered_map <int, std::pair<int,int>> &clusterLimits);//auxilliary function of threadHaplotypes

std::vector<int> rescue_reads(std::vector<int> &threadedClusters, std::vector<Column> &snps, std::vector<size_t> &suspectPostitions);

#endif