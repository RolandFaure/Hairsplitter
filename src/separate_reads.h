#ifndef SEP_READS_H
#define SEP_READS_H

#include <string>
#include <vector>
#include <Eigen/Sparse>

#include "read.h"
#include "Partition.h"
#include "robin_hood.h"
#include "tools.h"

void parse_column_file(
    std::string file, 
    std::vector<std::vector<Column>>  &snps, 
    std::unordered_map<std::string, int> &index_of_names,
    std::unordered_map<int, std::string> &name_of_contigs,
    std::vector<std::vector<std::string>> &names_of_reads,
    std::vector<long int> &length_of_contigs,
    std::vector<double> &coverage_of_contigs,
    std::vector<std::vector<std::pair<int,int>>>& readLimits,
    std::vector<int>& numberOfReads,
    int max_coverage,
    float rarest_strain_abundance);

void list_similarities_and_differences_between_reads(
    std::vector<Column> &snps, 
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs);

void list_similarities_and_differences_between_reads2(
    std::vector<Column> &snps, 
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs);

void list_similarities_and_differences_between_reads3(
    std::vector<Column> &snps, 
    Eigen::SparseMatrix<int>& similarity,
    Eigen::SparseMatrix<int>& difference);

void create_read_graph(
    std::vector <bool> &mask,
    int chunk, 
    int sizeOfWindow,
    std::vector<std::vector<std::pair<int,int>>> &sims_and_diffs,
    std::vector<std::vector<int>> &adjacency_matrix, 
    float &errorRate);

void create_read_graph_low_memory(
    std::vector<Column> &snps,
    std::vector <bool> &mask,
    int chunk,
    int sizeOfWindow,
    std::vector<std::vector<std::pair<int,int>>> &adjacency_matrix, //containing only the 1s
    float &errorRate);

void create_read_graph_matrix(
    std::vector <bool> &mask,
    int chunk,
    int sizeOfWindow,
    Eigen::SparseMatrix<int>& similarity,
    Eigen::SparseMatrix<int>& difference,
    Eigen::SparseMatrix<int> &adjacency_matrix, //containing only the 1s
    float &errorRate);

void finalize_clustering( 
    std::vector<Column> &snps,  
    std::vector<std::vector<int>> &localClusters, 
    std::vector<std::vector<int>> &strengthened_neighbor_list, 
    Eigen::SparseMatrix<int> &adjacency_matrix,
    bool low_memory,
    std::vector<bool> &mask_at_this_position,
    std::vector<int> &haplotypes,
    float errorRate,
    int posstart,
    int posend);

std::vector<int> merge_wrongly_split_haplotypes(
    std::vector<int> &clusteredReads, 
    std::vector<Column> &snps, 
    std::vector<std::vector<int>> &neighbor_list,
    Eigen::SparseMatrix<int> &adjacency_matrix,
    bool low_memory,
    int posstart,
    int posend);

std::vector<int> merge_clusterings(
    std::vector<std::vector<int>> &localClusters,
    std::vector<std::vector<int>> &neighbor_list,
    Eigen::SparseMatrix<int> &adjacency_matrix, 
    bool low_memory, 
    std::vector <bool> &mask);


#endif


