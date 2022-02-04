#ifndef READ_H
#define READ_H

#include <vector>
//#include <thread>
//#include <mutex>

#include "sequence.h"

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
//private :

    std::vector <long int> neighbors_; //list of long int referring to indices in allOverlaps

};


#endif // READ_H
