#include <vector>
#include <set>
#include <unordered_map>

#include "tools.h"
#include "robin_hood.h"

//low-memory version
std::vector<std::vector<int>> strengthen_adjacency_matrix(std::vector<std::vector<int>> &neighbor_list, int size);
//high-memory version
std::vector<std::vector<int>> strengthen_adjacency_matrix_high_memory(std::vector<std::vector<int>> &adjacency_matrix);


std::vector<int> chinese_whispers(std::vector<std::vector<int>> &neighbor_list, std::vector<int> &initialClusters, std::vector<bool> &mask);
std::vector<int> chinese_whispers_high_memory(std::vector<std::vector<int>> &adjacency_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask);
// std::vector<int> chinese_whispers(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &initialClusters);

void merge_close_clusters(std::vector<std::vector<int>> &neighbor_list, std::vector<std::vector<int>> &adjacency_matrix_high_memory, bool low_memory, std::vector<int> &clusters, std::vector<bool> &mask);
