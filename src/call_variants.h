#ifndef CALL_VARIANTS_H
#define CALL_VARIANTS_H

#include <string>
#include <vector>

#include "read.h"
#include "robin_hood.h"
#include "tools.h"
#include "Partition.h"

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

std::vector<Column> call_variants(
    std::vector<Column> &snps, 
    std::vector <Read> &allreads,
    std::vector <Overlap> &allOverlaps, 
    long int contig,
    std::string &ref,
    std::vector<size_t> &suspectPostitions,
    float &meanError, 
    float automatic_snp_threshold,
    std::vector<Column>  &automatic_snps, 
    std::string &tmpFolder,
    bool DEBUG);


void keep_only_robust_variants(
    std::vector<Column> &msa,
    std::vector<Column>  &snps_in, 
    std::vector<Column>  &snps_out, 
    float mean_error,
    std::vector <Partition> &parts);

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

std::vector<Column> rescue_snps(
    std::vector<Column> &msa, 
    float meanDistance, 
    std::vector<Partition> &partitions,
    std::vector<size_t> &suspectPostitions);

distancePartition distance(Partition &par1, Column &par2, char ref_base);
distancePartition distance(Partition &par1, Partition &par2, int threshold_p);
float computeChiSquare(distancePartition dis);

void output_files(std::unordered_map<int, std::vector<Column>> &variants, 
    std::vector<Overlap> &allOverlaps, 
    std::vector<Read> &allreads, 
    std::string col_file, 
    std::string vcf_file);

#endif