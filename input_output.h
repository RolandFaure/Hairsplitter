#ifndef IO
#define IO

#include <iostream>
#include <vector>
#include "robin_hood.h"
#include "read.h"

void parse_reads(std::string fileReads, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices);
void parse_PAF(std::string filePAF, std::vector <Overlap>& allOverlaps, std::vector <Read> &allreads, robin_hood::unordered_map<std::string, unsigned long int> &indices, std::vector<unsigned long int> &backbones_reads);

#endif