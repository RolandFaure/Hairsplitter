#ifndef REASSEMBLE_UNALIGNED_READS_H
#define REASSEMBLE_UNALIGNED_READS_H

#include <string>
#include <vector>

#include "read.h"

void reassemble_unaligned_reads(std::vector <Read> &allreads, std::vector <Overlap> &allOverlaps, std::string &fileReads,
                    std::vector<unsigned long int> &backbones_reads, std::string outputFolder, int num_threads); 

void assemble_with_wtdbg2(std::string &fileReads, std::string outputFolder, int num_threads);


#endif