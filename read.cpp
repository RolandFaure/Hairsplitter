#include "read.h"

using std::vector;
using std::cout;
using std::endl;

Read::Read(std::string s)
{
    sequence_ = Sequence(s);
    name = "";
}

void Read::add_overlap(long int o){

    neighbors_.push_back(o);

}

size_t Read::size(){
    return sequence_.size();
}

void Read::new_backbone(std::pair<int, int> pair, size_t size){

    if (pair.second < size){
        backbone_seq.push_back(pair);
    }
    else{
        cout << pair.second << " " << size << endl;
        throw std::logic_error("Problem in backbone, too high neighbor index: ");
    }

}