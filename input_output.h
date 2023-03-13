#ifndef IO
#define IO

#include <iostream>
#include <vector>
#include "check_overlaps.h"
#include "robin_hood.h"
#include "read.h"
//#include "Variant.h"

void parse_reads(std::string fileReads, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices);
void parse_reads_on_contig(std::string fileReads, long int backbone, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads);
void parse_assembly(std::string fileAssembly, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices,
     std::vector<unsigned long int> &backbone_reads, std::vector<Link> &allLinks);

void parse_PAF(std::string filePAF, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, 
    std::vector<unsigned long int> &backbones_reads, bool computeBackbones);

void parse_SAM(std::string fileSAM, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices);

// void parse_VCF(std::string fileVCF, robin_hood::unordered_map<std::string, std::vector <Variant>> &allVariants);
// void parseSAM(std::string fileSAM , robin_hood::unordered_map<std::string, std::vector <Variant>> &allVariants);
        
void output_GAF(std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, std::vector<Link> &allLinks, 
    std::vector <Overlap> &allOverlaps, std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions,
    std::string outputGAF);

void output_FASTA(std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, std::string fileOut);
void output_GFA(std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, std::string fileOut, std::vector<Link> &allLinks);
void output_filtered_PAF(std::string fileOut, std::string fileIn, std::vector <Read> &allreads, std::vector<std::vector<int>> &partitions, robin_hood::unordered_map<std::string, unsigned long int> &indices);
void output_readGroups(std::string readGroupsFile, std::vector <Read> &allreads, std::vector<unsigned long int> &backbone_reads, 
    std::unordered_map<unsigned long int ,std::vector< std::pair<std::pair<int,int>, std::pair<std::vector<int>, std::unordered_map<int, std::string>>  > >> &partitions, std::vector <Overlap> &allOverlaps);


//debug functions
void outputGraph(std::vector<std::vector<int>> &adj , std::vector<int> &clusters, std::string fileOut);
void output_simple_graph(std::vector<std::vector<int>> &adj, std::string fileOut);
void outputGraph_several_clusterings(std::vector<std::vector<int>> &adj,std::vector<std::vector<int>> &clusters, std::string fileOut);
void outputMatrix(std::vector<Column> &snps, std::vector<size_t> suspectPostitions, std::string id);
void outputTruePar(Partition truePar, std::string id);
#endif