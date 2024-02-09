#ifndef ALIGN_H
#define ALIGN_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "robin_hood.h"

std::string align(std::string &s1, size_t start1, size_t end1, std::string &s2, size_t start2, size_t end2);

std::string convert_CIGAR(std::string &cigar);
std::string polish(std::string &toPolish, std::vector<std::string> &reads, std::string& path_minimap2, std::string& path_racon);

#endif 