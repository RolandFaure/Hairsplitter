#ifndef FILTER_VARIANTS_H
#define FILTER_VARIANTS_H

#include <vector>
#include <string>
#include <unordered_map>

#include "read.h"
#include "Partition.h"
#include "robin_hood.h"

//struct to return a slightly complicated result
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

void parse_column_file(
    std::string file, 
    std::vector<std::vector<Column>> &snps, 
    std::unordered_map<std::string, int>& name_of_contigs, 
    std::unordered_map<int, std::string>& name_of_contigs2,
    std::vector<int> &numberOfReads);

void output_new_column_file(
    std::string initial_column_file, 
    std::vector<std::vector<Column>>  &new_snps, 
    std::unordered_map<std::string, int>& name_of_contigs, 
    std::string &output_file);

void output_new_vcf_file(
    std::string &initial_vcf_file, 
    std::vector<std::vector<Column>>  &new_snps, 
    std::unordered_map<std::string, int>& name_of_contigs, 
    std::string &output_file);


void keep_only_robust_variants(
    std::vector<std::vector<Column>> &snps_in, 
    std::vector<std::vector<Column>>  &snps_out, 
    std::unordered_map<std::string, int>& name_of_contigs,
    float mean_error,
    int num_threads);

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

distancePartition distance(Partition &par1, Column &par2, char ref_base);
distancePartition distance(Partition &par1, Partition &par2, int threshold_p);
float computeChiSquare(distancePartition dis);

#endif