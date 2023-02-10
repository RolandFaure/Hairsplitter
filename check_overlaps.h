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
    int solid11;
    int solid10;
    int solid01;
    int solid00;
    float score;
    short phased; // worth -1 or 1
    bool augmented; //to know if the partition was augmented or not
    char secondBase; //the second base in the partition apart from the reference
    Column partition_to_augment;
};


void checkOverlaps(
    std::string fileReads,
    std::vector <Read> &allreads, 
    std::vector <Overlap> &allOverlaps, 
    std::vector<unsigned long int> &backbones_reads, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions,
    bool assemble_on_assembly,
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits,
    bool polish, 
    int num_threads);

void compute_partition_on_this_contig(
    std::string fileReads, 
    long int contig,std::vector <Read> &allreads, 
    std::vector <Overlap> &allOverlaps, 
    std::vector<unsigned long int> &backbones_reads, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions, 
    bool assemble_on_assembly,
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits,
    bool polish);

float generate_msa(
    long int bbcontig, 
    std::vector <Overlap> &allOverlaps, 
    std::vector <Read> &allreads, 
    std::vector<Column> &snps, 
    robin_hood::unordered_map<int, int> &insertionPos,
    int backboneReadIndex, 
    std::string &truePar, 
    bool assemble_on_assembly, 
    std::unordered_map <int, std::vector<std::pair<int,int>>> &readLimits, 
    std::vector<bool>& misalignedReads, 
    bool polish,
    std::string &newref);

std::string consensus_reads(
    std::string &backbone, 
    std::vector <std::string> &polishingReads, 
    int overhangLeft,
    int overhangRight,
    std::string &id);

void compute_consensus_in_partitions(
    long int contig, 
    std::vector<std::pair<std::pair<int,int>, std::vector<int>> > &partition, 
    std::vector <Read> &allreads, 
    std::vector <Overlap> &allOverlaps, 
    std::vector<Column> &snps, 
    robin_hood::unordered_map<int, int> &insertionPositions,
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions);
// std::string local_assembly(std::vector <std::string> &reads);

std::vector< std::pair<std::pair<int,int>, std::vector<int>> > separate_reads(
    std::string& ref, 
    std::vector<Column> &snps,
    int numberOfReads);

std::vector<Partition> get_solid_partitions(
    std::string& ref, 
    std::vector<Column> &snps,
    std::vector<bool> &mask,
    std::vector<size_t> &suspectPostitions,
    float &meanError,
    int numberOfReads);

std::vector<Partition> second_pass(
    std::string& ref,
    std::vector<Partition> &partitions, 
    std::vector<Column> &snps,
    std::vector<size_t> &suspectPostitions,
    float &meanError,
    int numberOfReads);

std::vector <std::vector <bool>> create_masks(
    std::vector<Partition> &partitions, 
    int numberOfReads,
    std::vector<std::vector<bool>> &all_mask_ever_tested);

distancePartition distance(Partition &par1, Column &par2, char ref_base);
distancePartition distance(Partition &par1, Partition &par2, int threshold_p);
float computeChiSquare(distancePartition dis);

std::vector<Partition> select_compatible_partitions(
    std::vector<Partition> &partitions, 
    int numberOfReads,
    float errorRate);

std::vector<Partition> select_confident_partitions(
    std::vector<Partition> &partitions, 
    std::vector<bool> trimmedListOfFinalPartitionBool, 
    int numberOfReads,
    float errorRate,
    int numberOfSuspectPositions);

    
int compatible_partitions(Partition &p1 , Partition &p2);

std::vector<int> threadHaplotypes_in_interval(
    std::vector<Partition> &listOfFinalPartitions,
    int numberOfReads);

std::vector<int> rescue_reads(
    std::vector<int> &clusteredReads, 
    std::vector<Column> &snps, 
    int chunk, 
    std::vector<size_t> &suspectPostitions,
    int sizeOfWindow);

bool extend_with_partition_if_compatible(std::vector<int> &alreadyThreadedHaplotypes, Partition &extension, int partitionIndex,
        std::unordered_map <int, std::pair<int,int>> &clusterLimits);//auxilliary function of threadHaplotypes

#endif