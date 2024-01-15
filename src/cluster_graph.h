#include <vector>
#include <set>
#include <unordered_map>
#include <Eigen/Sparse>

#include "tools.h"
#include "robin_hood.h"

//low-memory version
std::vector<std::vector<int>> strengthen_adjacency_matrix(std::vector<std::vector<int>> &neighbor_list, int size);
//high-memory version
std::vector<std::vector<int>> strengthen_adjacency_matrix_high_memory(std::vector<std::vector<int>> &adjacency_matrix);
//matrix version
Eigen::SparseMatrix<int> strengthen_adjacency_matrix_high_memory2(Eigen::SparseMatrix<int> &adjacency_matrix);


std::vector<int> chinese_whispers(std::vector<std::vector<int>> &neighbor_list, std::vector<int> &initialClusters, std::vector<bool> &mask);
std::vector<int> chinese_whispers_high_memory(Eigen::SparseMatrix<int> &adjacency_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask);
std::vector<int> chinese_whispers_matrix(Eigen::SparseMatrix<int> &adj_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask);

void merge_close_clusters(std::vector<std::vector<int>> &neighbor_list, Eigen::SparseMatrix<int> &adjacency_matrix, bool low_memory, std::vector<int> &clusters, std::vector<bool> &mask);
