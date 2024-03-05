#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "robin_hood.h"


struct Mapping{
    int length_of_read;
    std::string read;

    //contig on the left of the mapping
    std::string contig1;
    int position_on_contig1;
    bool orientation_on_contig1;
    int pos_on_read1;
    bool orientation_on_read1;
    bool breakpoint1;

    //contig on the right of the mapping
    std::string contig2;
    int position_on_contig2;
    bool orientation_on_contig2;
    int pos_on_read2;
    bool orientation_on_read2;
    bool breakpoint2;
};

struct Bridge{ //a bridge is a bridge between two contigs defined by a read
    std::string contig1;
    std::string contig2;
    int position1;
    int position2;
    bool strand1; //true if bridge toward the right of contig1, false if bridge toward the left
    bool strand2; //true if bridge toward the right of contig2, false if bridge toward the left

    std::string read_name;
    int pos_read_on_contig1;
    int pos_read_on_contig2;
};

struct SolidBridge{ //a bridge is a bridge between two contigs defined by several read, e.g. by agregating Bridge
    std::string contig1;
    std::string contig2;
    int position1;
    int position2;
    bool strand1; //true if bridge toward the right of contig1, false if bridge toward the left
    bool strand2; //true if bridge toward the right of contig2, false if bridge toward the left

    std::vector<std::string> read_names;
    std::vector<int> pos_read_on_contig1;
    std::vector<int> pos_read_on_contig2;
    std::vector<bool> strand;
};

struct Link{ //a link is a link between two contigs, containing a CIGAR and optionnaly an extra sequence to insert between the two contigs
    std::string contig1;
    std::string contig2;
    int position1;
    int position2;
    bool strand1; //true if bridge toward the right of contig1, false if bridge toward the left
    bool strand2; //true if bridge toward the right of contig2, false if bridge toward the left

    std::string cigar;
    std::string extra_sequence;

    int coverage;
};

#endif 