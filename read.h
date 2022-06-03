#ifndef READ_H
#define READ_H

#include <vector>
#include <list>
//#include <thread>
//#include <mutex>

#include "sequence.h"
#include "Partition.h"

struct Overlap
{
    unsigned long int sequence1;
    unsigned long int sequence2;
    int position_1_1;
    int position_1_2;
    int position_2_1;
    int position_2_2;
    bool strand; //false if the two reads are on different strands
};

struct Link{
    std::string CIGAR;
    long int neighbor1;
    short end1; //either 0 or 1
    long int neighbor2;
    short end2; //either 0 or 1

    int group; //optional, useful in modify_gfa
};

class Read
{

public:
    Read();
    Read(std::string s);
    Sequence sequence_;

    std::string name;
    float depth; //coverage of the contig if indicated in the gfa
    std::string comments; //comments in the GFA. We don't want to lose them !

    void add_overlap(long int o);
    size_t size();
    void new_backbone(std::pair<int, int> pair, size_t size);
    std::vector <long int> neighbors_; //list of long int referring to indices in allOverlaps

    void add_link(size_t l, short end); //that's when reads are actually contigs
    std::vector<size_t> get_links_left();
    std::vector<size_t> get_links_right();
//private :
    std::vector <std::pair<int, short>> backbone_seq; //first element is the id of backbone read, second is the indice of neighbor

private : //when reads are actually contigs
    std::vector<size_t> links_left; //indices of links in allLinks vector
    std::vector<size_t> links_right;
};



#endif // READ_H
