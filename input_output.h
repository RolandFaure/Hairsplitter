#ifndef IO
#define IO

#include <iostream>
#include <vector>
#include "check_overlaps.h"
#include "robin_hood.h"
#include "read.h"

void parse_reads(std::string fileReads, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices);
void parse_assembly(std::string fileAssembly, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, std::vector<unsigned long int> &backbone_reads);

void parse_PAF(std::string filePAF, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, std::vector<unsigned long int> &backbones_reads, bool computeBackbones);

void output_filtered_PAF(std::string fileOut, std::string fileIn, std::vector <Read> &allreads, std::vector<Partition> &partitions, robin_hood::unordered_map<std::string, unsigned long int> &indices);
#endif