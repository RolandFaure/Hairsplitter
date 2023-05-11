#include "cluster_graph.h"
#include <unordered_map>
#include <random>
#include <algorithm>
#include <iostream>
#include <set>
#include <assert.h>

using std::vector;
using std::unordered_map;
using std::pair;
using std::max;
using std::min;
using std::cout;
using std::endl;
using std::set;

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
    int number_of_iterations = 0;

    //iterate until less than 2 changes are made  
    while (changes > 2 && number_of_iterations < 15){

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
                if (adjacency_matrix[i][j] >= 1 && clusters[j] >= 0){
                    neighbors[clusters[j]]+= adjacency_matrix[i][j];
                }
            }

            int max_index = 0;
            int max_value = 0;
            for (int j = 0; j < neighbors.size(); j++){
                if (neighbors[j] > max_value){
                    max_value = neighbors[j];
                    max_index = j;
                }
            }

            //choose the new cluster among the clusters that have more than one third of the links 
            if (max_value > 0){
                if (clusters[i] != max_index){
                    changes++;
                }
                clusters[i] = max_index;
            }
        }

        // cout << "dopudqsop number oc changes " << changes << " ";
        // for (auto j = 0 ; j < clusters.size() ; j++){
        //     if (mask[j]){
        //         cout << clusters[j] << " ";
        //     }
        // }
        // cout << endl;
        number_of_iterations += 1;
    }

    // cout << "done clustddeering" << endl;

    //put -2 for masked reads
    for (int i = 0; i < mask.size(); i++){
        if (!mask[i]){
            clusters[i] = -2;
        }
    }

    //return the clusters
    return clusters;
}

/**
 * @brief Merge close clusters by disrupting clusters and applying chinese whispers
 * 
 * @param adjacency_list Graph
 * @param clusters initial clusters
 * @param mask Masked reads
 * @return std::vector<int> merged clusters
 */
void merge_close_clusters(std::vector<std::vector<int>> const &adjacency_matrix, std::vector<int> &clusters, std::vector<bool> &mask){

    std::set<int> setOfTestedClusters;
    vector<int> initialCountOfClusters (clusters.size(), 0);
    for (auto i : clusters){
        if (i >= 0 && i < initialCountOfClusters.size()){
            initialCountOfClusters[i] += 1;
        }
    }

    for (auto node = 0 ; node < clusters.size() ; node++){
        if (setOfTestedClusters.find(clusters[node])==setOfTestedClusters.end() && clusters[node] >= 0){
            int clusterToTest = clusters[node];

            //go through the graph and convert nodes that have less than 2/3 links to the cluster
            auto newclusters = clusters;
            int changes = 3;
            auto countOfClusters = initialCountOfClusters;
            auto number_of_iterations = 0;
            //iterate until less than 2 changes are made  
            while (changes > 0 && number_of_iterations < 10){

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
                    if (!mask[i] || newclusters[i] != clusterToTest){
                        continue;
                    }          
                    vector<int> neighbors (adjacency_matrix.size(), 0);
                    for (int j = 0; j < adjacency_matrix[i].size(); j++){
                        if (adjacency_matrix[i][j] >= 1 && newclusters[j] >= 0){
                            neighbors[newclusters[j]]+= adjacency_matrix[i][j];
                        }
                    }

                    int max_index = 0;
                    int max_value = 0;
                    int second_index = 0;
                    int second_value = 0;
                    for (int j = 0; j < neighbors.size(); j++){
                        if (neighbors[j] > max_value){
                            second_value = max_value;
                            second_index = max_index;
                            max_value = neighbors[j];
                            max_index = j;
                        }
                        else if (neighbors[j] > second_value){
                            second_value = neighbors[j];
                            second_index = j;
                        }
                    }

                    //choose the new cluster among the clusters that have more than one third of the links 
                    if (max_value > 0 && max_index != clusterToTest){
                        countOfClusters[newclusters[i]]--;
                        countOfClusters[max_index]++;
                        changes++;
                        newclusters[i] = max_index;
                    }
                    else if (max_value > 0 && max_value <= 2*second_value) {
                        //then this node is weak, change the value
                        // cout << "max and second values : " << max_value << " " << second_value << " " << max_index << " " << second_index << endl;
                        // assert(newclusters[i] < countOfClusters.size() && "ERROR in cluster_graph.cpp, error code 569987");
                        // assert(newclusters[i] >= 0 && "ERROR in cluster_graph.cpp, error code 569987");
                        countOfClusters[newclusters[i]]--;
                        countOfClusters[second_index]++;
                        newclusters[i] = second_index;
                        changes++;
                    }

                }
                number_of_iterations++;
            }

            setOfTestedClusters.emplace(clusterToTest);
            if (countOfClusters[clusterToTest] == 0){ //the cluster was weak, it must be suppressed
                clusters = newclusters;
                initialCountOfClusters = countOfClusters;
            }
        }
    }
}




