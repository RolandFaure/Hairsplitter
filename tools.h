#ifndef TOOLS
#define TOOLS

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdio>

#include "Partition.h"

std::string convert_cigar(std::string &cigar);
std::string reverse_complement(std::string &seq);
void print_alignment(std::string &ref, std::string &read, std::string &cigar, int start, int end);
int compute_edit_distance(std::string &cigar, std::string &ref, std::string &read, int start, int end);

#endif