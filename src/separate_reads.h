#ifndef SEP_READS_H
#define SEP_READS_H

#include <string>
#include <vector>

#include "read.h"
#include "Partition.h"
#include "robin_hood.h"

void parse_column_file(
    std::string file, 
    std::vector<std::vector<Column>>  &snps, 
    std::unordered_map<std::string, int> &index_of_names,
    std::unordered_map<int, std::string> &name_of_contigs,
    std::vector<std::vector<std::string>> &names_of_reads,
    std::vector<long int> &length_of_contigs,
    std::vector<std::vector<std::pair<int,int>>>& readLimits,
    std::vector<int>& numberOfReads);

void list_similarities_and_differences_between_reads(
    std::vector<Column> &snps, 
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs);

void create_read_graph(
    std::vector <bool> &mask,
    std::vector<Column> &snps, 
    int chunk, 
    int sizeOfWindow,
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs,
    std::vector< std::vector<int>> &adjacency_matrix,
    float &errorRate);

std::vector<int> merge_clusterings(std::vector<std::vector<int>> &localClusters,
    std::vector< std::vector<int>> &adjacency_matrix, std::vector <bool> &mask);


#endif


