#include "cluster_graph.h"
#include <unordered_map>
#include <random>
#include <algorithm>
#include <iostream>

using std::vector;
using std::unordered_map;
using std::pair;
using std::max;
using std::min;
using std::cout;
using std::endl;

/**
 * @brief Strengthen the adjacency matrix by adding edges between nodes that are connected by a path of length 2
 * 
 * @param adjacency_matrix 
 * @return std::vector<std::vector<int>> 
 */
std::vector<std::vector<int>> strengthen_adjacency_matrix(std::vector<std::vector<int>> const &adjacency_matrix){

    //create a copy of the adjacency matrix
    std::vector<std::vector<int>> strengthened_adjacency_matrix(adjacency_matrix);

    //iterate over all nodes
    for (int i = 0; i < adjacency_matrix.size(); i++){
        //iterate over all nodes that are connected to i
        vector<int> neighbors;
        for (int j = 0; j < adjacency_matrix[i].size(); j++){
            //if i and j are connected
            if (adjacency_matrix[i][j] == 1){
                neighbors.push_back(j);
            }
        }
        //add 1 to the adjacency matrix for all nodes that are connected through i
        for (int j = 0; j < neighbors.size(); j++){
            for (int k = j+1; k < neighbors.size(); k++){
                strengthened_adjacency_matrix[neighbors[j]][neighbors[k]] += 1;
                strengthened_adjacency_matrix[neighbors[k]][neighbors[j]] += 1;
            }
        }
    }

    return strengthened_adjacency_matrix;

}

/**
 * @brief 
 * 
 * @param adjacency_matrix Graph
 * @param initialClusters Starting point of the clustering
 * @return std::vector<int> 
 */
// std::vector<int> chinese_whispers(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &initialClusters){
//     vector<bool> mask(adjacency_matrix.size(), true);
//     return chinese_whispers(adjacency_matrix, initialClusters, mask);
// }


/**
 * @brief Implement the chinese whispers clustering algorithm on adjacency matrix
 * 
 * @param adjacency_list Graph
 * @param initialClusters Starting point of the clustering
 * @param mask Masked reads
 * @return std::vector<int> 
 */
std::vector<int> chinese_whispers(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask){

    //at the beginning, all nodes are in their own cluster
    // std::vector<int> clusters(adjacency_matrix.size());
    // for (int i = 0; i < adjacency_matrix.size(); i++){
    //     clusters[i] = i;
    // }
    auto clusters = initialClusters;

    //keep track of the number of changes in the clustering
    int changes = 3;
    //iterate until less than 2 changes are made  
    while (changes > 2){

        changes = 0;
        //iterate over all nodes in a random order
        // for (int i = 0; i < adjacency_matrix.size(); i++){
        // generate a random order
        vector<int> order(adjacency_matrix.size());
        std::iota(order.begin(), order.end(), 0);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(order.begin(), order.end(), g);
    
        for (int i : order){
            //find the most frequent neighbor
            if (!mask[i]){
                continue;
            }          
            vector<int> neighbors (adjacency_matrix.size(), 0);
            for (int j = 0; j < adjacency_matrix[i].size(); j++){
                if (adjacency_matrix[i][j] >= 1){
                    neighbors[clusters[j]]+= adjacency_matrix[i][j];
                }
            }
            neighbors[-1] = 0; //the cluster -1 is not allowed, but nodes with cluster -1 can change to another cluster

            int max_index = 0;
            int max_value = 0;
            for (int j = 0; j < neighbors.size(); j++){
                if (neighbors[j] > max_value){
                    max_value = neighbors[j];
                    max_index = j;
                }
            }

            //if the most frequent neighbor is not the current node, change the cluster
            if (clusters[i] != max_index && max_value > 0){
                clusters[i] = max_index;
                changes++;
            }
        }

    }

    //put -2 for masked reads
    for (int i = 0; i < mask.size(); i++){
        if (!mask[i]){
            clusters[i] = -2;
        }
    }

    //return the clusters
    return clusters;
}





