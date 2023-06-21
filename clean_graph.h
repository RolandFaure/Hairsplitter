#ifndef CLEAN_GRAPH_H
#define CLEAN_GRAPH_H

#include "input_output.h"
#include "tools.h"
#include <string>

void clean_graph(
    std::string& assemblyFile,
    std::string& outputFile,
    std::string& logFile,
    int num_threads,
    std::string &outFolder,
    int &nb_of_deleted_contigs,
    int &length_of_deleted_contigs);

#endif
