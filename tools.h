#ifndef TOOLS
#define TOOLS

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdio>

#include "edlib.h"
#include "Partition.h"

std::string convert_cigar(std::string &cigar);
std::string reverse_complement(std::string &seq);
void print_alignment(std::string &ref, std::string &read, std::string &cigar, int start, int end);
int compute_edit_distance(std::string &cigar, std::string &ref, std::string &read, int start, int end);

void convert_GFA_to_FASTA(std::string &gfa_file, std::string &fasta_file);
void convert_FASTA_to_GFA(std::string &fasta_file, std::string &gfa_file);

void rename_reads(std::string &fasta_file, std::string &prefix);

std::string consensus_reads(
    std::string const &backbone, 
    std::vector <std::string> &polishingReads, 
    std::string &id);

#endif