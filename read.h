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

class Read
{

public:
    Read(std::string s);
    Sequence sequence_;

    void add_overlap(long int o);
    size_t size();
    void new_backbone(std::pair<int, short> pair, size_t size);
    std::vector <long int> neighbors_; //list of long int referring to indices in allOverlaps

//private :
    std::vector <std::pair<int, short>> backbone_seq; //first element is the id of backbone read, second is the indice of neighbor

};


#endif // READ_H
