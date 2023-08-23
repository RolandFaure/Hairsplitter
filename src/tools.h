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

void convert_GFA_to_FASTA(std::string &gfa_file, std::string &fasta_file);
void convert_FASTA_to_GFA(std::string &fasta_file, std::string &gfa_file);

void rename_reads(std::string &fasta_file, std::string &prefix);

std::string consensus_reads(
    std::string &backbone, 
    std::string &full_backbone, 
    int start_pos_on_full_backbone,
    int sizeOfWindow,
    std::vector <std::string> &polishingReads, 
    std::string &id,
    std::string &outFolder,
    std::string &techno,
    std::string &MINIMAP, 
    std::string &RACON);

bool check_alignment(std::string &paf_file);

std::string consensus_reads_wtdbg2(
    std::string const &backbone, 
    std::vector <std::string> &polishingReads, 
    std::string &id,
    std::string &outFolder,
    std::string &techno,
    std::string &MINIMAP, 
    std::string &RACON,
    std::string &SAMTOOLS,
    std::string &WTDBG2
);

void assemble_with_wtdbg2(std::string &fileReads, std::string outputFolder, std::string &ref, std::string id, std::string &WTDBG2, std::string &MINIMAP, std::string &SAMTOOLS);

#endif
