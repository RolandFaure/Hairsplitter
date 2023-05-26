#include <vector>

std::vector<std::vector<int>> strengthen_adjacency_matrix(std::vector<std::vector<int>> const &adjacency_matrix);

std::vector<int> chinese_whispers(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask);
// std::vector<int> chinese_whispers(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &initialClusters);

void merge_close_clusters(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &clusters, std::vector<bool> &mask);
