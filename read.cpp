#include "read.h"

using std::vector;
using std::cout;
using std::endl;

Read::Read(std::string s)
{
    sequence_ = Sequence(s);
}

void Read::add_overlap(long int o){

    neighbors_.push_back(o);

}

size_t Read::size(){
    return sequence_.size();
}