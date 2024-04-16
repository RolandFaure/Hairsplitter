#ifndef TOOLS
#define TOOLS

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdio>

#include "Partition.h"

// A hash function used to hash a pair of any kind
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);

        // return size_t(hash1) << 32 | hash2;
 
        if (hash1 != hash2) {
            return hash1 ^ hash2;              
        }
         
        // If hash1 == hash2, their XOR is zero.
        return hash1;
    }
};

struct equal_pair {
    template <class T1, class T2>
    bool operator () (const std::pair<T1, T2>& p1, const std::pair<T1, T2>& p2) const {
        return p1.first == p2.first && p1.second == p2.second;
    }
};


std::string convert_cigar(std::string &cigar);
std::string convert_cigar2(std::string &cigar);

std::string reverse_complement(std::string &seq);
void print_alignment(std::string &ref, std::string &read, std::string &cigar, int start, int end);
int compute_edit_distance(std::string &cigar, std::string &ref, std::string &read, int start, int end);

void convert_GFA_to_FASTA(std::string &gfa_file, std::string &fasta_file);
void convert_FASTA_to_GFA(std::string &fasta_file, std::string &gfa_file);
void cut_GFA(const std::string& input_assembly, const std::string& output_assembly, int max_length);

void rename_reads(std::string &fasta_file, std::string &prefix);

std::string consensus_reads(
    std::string &backbone, 
    std::string &full_backbone, 
    int start_pos_on_full_backbone,
    int sizeOfWindow,
    std::vector <std::string> &polishingReads, 
    std::vector <std::string> &fullReads,
    std::vector <std::pair<std::string,int>> &CIGARs,
    std::string &id,
    std::string &outFolder,
    std::string &techno,
    std::string &MINIMAP, 
    std::string &RACON,
    std::string &path_to_python,
    std::string &path_src);

int check_alignment(std::string &paf_file);
std::string alternative_backbone(std::string &sam_file, std::string &backbone);

std::string consensus_reads_medaka(
    std::string const &backbone, 
    std::vector <std::string> &polishingReads, 
    std::string &id,
    std::string outFolder,
    std::string &MEDAKA,
    std::string &SAMTOOLS,
    std::string &path_to_python,
    std::string &path_src);

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

std::string basic_assembly(std::string read_file, std::string &MINIMAP, std::string &tmp_folder, std::string &id);

void assemble_with_wtdbg2(std::string &fileReads, std::string outputFolder, std::string &ref, std::string id, std::string &WTDBG2, std::string &MINIMAP, std::string &SAMTOOLS);

#endif
