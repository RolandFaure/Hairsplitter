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
using std::make_pair;

/**
 * @brief Strengthen the adjacency matrix by adding edges between nodes that are connected by a path of length 2
 * 
 * @param adjacency_matrix unordered map
 * @param size size of the adjacency matrix
 * @return unordered map
 */
std::vector<std::vector<int>> strengthen_adjacency_matrix(std::vector<std::vector<int>> &neighbor_list, int size){

    //create a copy of the adjacency matrix
    auto strengthened_neighbor_list = neighbor_list;

    //iterate over all nodes
    for (int i = 0; i < size; i++){
        //iterate over all nodes that are connected to i
        vector<int> neighbors;
        for (auto j : neighbor_list[i]){
            //if i and j are connected
            neighbors.push_back(j);
        }
        //add 1 to the adjacency matrix for all nodes that are connected through i
        for (int j = 0; j < neighbors.size(); j++){
            for (int k = j+1; k < neighbors.size(); k++){
                strengthened_neighbor_list[neighbors[j]].push_back(neighbors[k]);
                strengthened_neighbor_list[neighbors[k]].push_back(neighbors[j]);
            }
        }
    }

    //sort each entry of neighbor_list_low_memory and remove duplicates
    for (auto r = 0 ; r < strengthened_neighbor_list.size() ; r++){
        std::sort(strengthened_neighbor_list[r].begin(), strengthened_neighbor_list[r].end());
        //remove duplicate using unique
        strengthened_neighbor_list[r] = vector<int>(strengthened_neighbor_list[r].begin(), std::unique(strengthened_neighbor_list[r].begin(), strengthened_neighbor_list[r].end()));
    }

    return strengthened_neighbor_list;

}

/**
 * @brief Strengthen the adjacency matrix by adding edges between nodes that are connected by a path of length 2
 * 
 * @param adjacency_matrix 
 * @return std::vector<std::vector<int>> 
 */
std::vector<std::vector<int>> strengthen_adjacency_matrix_high_memory(vector<vector<int>> &adjacency_matrix){

    //create a copy of the adjacency matrix
    auto strengthened_adjacency_matrix = adjacency_matrix;

    //iterate over all nodes
    for (int i = 0; i < adjacency_matrix.size(); i++){
        //iterate over all nodes that are connected to i
        vector<int> neighbors;
        for (int j = 0; j < adjacency_matrix.size(); j++){
            //if i and j are connected
            if (adjacency_matrix[i][j] >= 1){
                neighbors.push_back(j);
            }
        }
        //add 1 to the adjacency matrix for all nodes that are connected through i
        for (int j = 0; j < neighbors.size(); j++){
            for (int k = j+1; k < neighbors.size(); k++){
                // if ((j==0 && k==1) || (j==1 && k==0)){
                //     cout << "nodes 0 and 1 are connected through " << i << endl;
                // }
                strengthened_adjacency_matrix[neighbors[j]][neighbors[k]] += 1;
                strengthened_adjacency_matrix[neighbors[k]][neighbors[j]] += 1;
            }
        }
    }

    return strengthened_adjacency_matrix;
}

/**
 * @brief Strengthen the adjacency matrix by adding edges between nodes that are connected by a path of length 2
 * 
 * @param adjacency_matrix 
 * @return 
 */
Eigen::SparseMatrix<int> strengthen_adjacency_matrix_high_memory2(Eigen::SparseMatrix<int> &adjacency_matrix){
    
    //create a triplet list
    std::vector<Eigen::Triplet<int>> tripletList;

    //iterate over all nodes
    for (int r = 0; r < adjacency_matrix.rows(); r++){
        //iterate over all nodes that are connected to i
        vector<int> neighbors;
        for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency_matrix, r); it; ++it){
            //row and col are connected
            neighbors.push_back(it.row());
        }
        //add 1 to the adjacency matrix for all nodes that are connected through i
        for (int j = 0; j < neighbors.size(); j++){
            for (int k = j+1; k < neighbors.size(); k++){
                tripletList.push_back(Eigen::Triplet<int>(neighbors[j], neighbors[k], 1));
                tripletList.push_back(Eigen::Triplet<int>(neighbors[k], neighbors[j], 1));
            }
        }
    }

    //create the new adjacency matrix
    Eigen::SparseMatrix<int> strengthened_adjacency_matrix(adjacency_matrix.rows(), adjacency_matrix.cols());
    strengthened_adjacency_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    return strengthened_adjacency_matrix + adjacency_matrix;
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
std::vector<int> chinese_whispers(std::vector<std::vector<int>> &neighbor_list, std::vector<int> &initialClusters, std::vector<bool> &mask){


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
        vector<int> order(neighbor_list.size());
        std::iota(order.begin(), order.end(), 0);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(order.begin(), order.end(), g);
    
        for (int i : order){
            //find the most frequent neighbor
            if (!mask[i]){
                continue;
            }          
            vector<int> neighbors (neighbor_list.size(), 0);
            for (auto n : neighbor_list[i]){
                if (clusters[n] >= 0){
                    neighbors[clusters[n]]+= 1;
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
 * @brief clusters a graph using the chinese whispers algorithm
 * 
 * @param adjacency_matrix 
 * @param initialClusters 
 * @param mask 
 * @return std::vector<int> 
 */
std::vector<int> chinese_whispers_high_memory(Eigen::SparseMatrix<int> &adjacency_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask){
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
        vector<int> order(initialClusters.size());
        std::iota(order.begin(), order.end(), 0);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(order.begin(), order.end(), g);
    
        for (int i : order){
            //find the most frequent neighbor
            if (!mask[i]){
                continue;
            }          
            vector<int> neighbors (mask.size(), 0);
            for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency_matrix, i); it; ++it){
                if (clusters[it.row()] >= 0){
                    neighbors[clusters[it.row()]]+= 1;
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
 * @brief clusters a graph using the chinese whispers algorithm
 * 
 * @param adjacency_matrix 
 * @param initialClusters 
 * @param mask 
 * @return std::vector<int> 
 */
std::vector<int> chinese_whispers_matrix(Eigen::SparseMatrix<int> &adj_matrix, std::vector<int> &initialClusters, std::vector<bool> &mask){


    std::random_device rd;
    std::mt19937 g(rd());

    vector<int> new_clusters = initialClusters;

    //convert the initial clusters to a sparse matrix, with a 1 on the row of the cluster
    Eigen::SparseMatrix<int> clusters(initialClusters.size(), initialClusters.size());
    for (int i = 0; i < initialClusters.size(); i++){
        if (initialClusters[i] >= 0 && mask[i]){
            clusters.insert(initialClusters[i],i) = 1;
        }
    }

    //keep track of the number of changes in the clustering
    int changes = 3;
    int number_of_iterations = 0;
    while (changes > 2 && number_of_iterations < 15){
        
        changes = 0;
        Eigen::SparseMatrix<int> clusters_of_neighbors = clusters * adj_matrix;
        Eigen::SparseMatrix<int> clusters_two_iterations_ago; //just to check if we're not oscillating
        
        // keep only the maximum value of each col and set it to 1
        for (int c = 0; c < clusters_of_neighbors.outerSize(); ++c){
            int max_value = 0;
            vector<int> max_indices;
            for (Eigen::SparseMatrix<int>::InnerIterator it(clusters_of_neighbors,c); it; ++it){
                if (it.value() > max_value){
                    max_value = it.value();
                    max_indices.clear();
                    max_indices.push_back(it.row());
                }
                else if (it.value() == max_value){
                    max_indices.push_back(it.row());
                }
            }
            //choose max_index randomly among the max_indices
            if (max_indices.size() == 0){
                continue;
            }
            int max_index = max_indices[std::uniform_int_distribution<int>(0, max_indices.size()-1)(g)];

            for (Eigen::SparseMatrix<int>::InnerIterator it(clusters_of_neighbors,c); it; ++it){
                if (it.row() != max_index){
                    it.valueRef() = 0;
                }
                else{
                    it.valueRef() = 1;
                    if (new_clusters[c] != it.row()){
                        changes++;
                        new_clusters[c] = it.row();
                    }
                }
            }
        }

        clusters = clusters_of_neighbors;

        number_of_iterations++;
    }

    //put -2 for masked reads
    for (int i = 0; i < mask.size(); i++){
        if (!mask[i]){
            new_clusters[i] = -2;
        }
    }

    return new_clusters;
}

/**
 * @brief Merge close clusters by disrupting clusters and applying chinese whispers
 * 
 * @param adjacency_list Graph
 * @param clusters initial clusters
 * @param mask Masked reads
 * @return std::vector<int> merged clusters
 */
void merge_close_clusters(std::vector<std::vector<int>> &neighbor_list, Eigen::SparseMatrix<int> &adjacency_matrix, bool low_memory, std::vector<int> &clusters, std::vector<bool> &mask){

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
            //iterate until less than 1 changes are made  
            while (changes > 0 && number_of_iterations < 10){

                changes = 0;
                //iterate over all nodes in a random order
                // for (int i = 0; i < adjacency_matrix.size(); i++){
                // generate a random order
                vector<int> order(clusters.size());
                std::iota(order.begin(), order.end(), 0);
                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(order.begin(), order.end(), g);
            
                for (int i : order){
                    //find the most frequent neighbor
                    if (!mask[i] || newclusters[i] != clusterToTest){
                        continue;
                    }          
                    vector<int> neighbors (mask.size(), 0);
                    if (low_memory){
                        for (int j = 0 ; j < neighbor_list[i].size() ; j++){
                            if (std::binary_search(neighbor_list[i].begin(), neighbor_list[i].end(),j) == true && newclusters[j] >= 0){
                                neighbors[newclusters[j]]+= 1;
                            }
                        }
                    }
                    else{
                        for (Eigen::SparseMatrix<int>::InnerIterator it(adjacency_matrix, i); it; ++it){
                            if (newclusters[it.row()] >= 0){
                                neighbors[newclusters[it.row()]]+= it.value();
                            }
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




