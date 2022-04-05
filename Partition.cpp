#include "Partition.h"

#include "robin_hood.h"
#include <iostream>
#include <cmath>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::min;
using std::max;
using std::begin;
using std::end;
using std::pair;
using std::make_pair;

Partition::Partition(int size){
    for (auto i = 0 ; i < size ; i++){
        mostFrequentBases.push_back(0);
        lessFrequence.push_back(0);
        moreFrequence.push_back(0);
    }
    numberOfOccurences = 0;
}

Partition::Partition(vector<char> & snp){

    robin_hood::unordered_flat_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['-'] = 4;
    bases2content['?'] = 5;

    int content [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for *, 5 for '-'
    for (int c = 0 ; c < snp.size() ; c++){
        if (snp[c] != '?'){
            content[bases2content[snp[c]]] += 1;
        }
    }

    char mostFrequent2 = 'A';
    int maxFrequence2 = content[0];
    char secondFrequent2 = 'C';
    int secondFrequence2 = content[1];
    if (content[0] < content[1]){
        mostFrequent2 = 'C';
        maxFrequence2 = content[1];
        secondFrequent2 = 'A';
        secondFrequence2 = content[0];
    }
    for (auto i = 2 ; i < 5 ; i++){
        if (content[i] > maxFrequence2){
            secondFrequence2 = maxFrequence2;
            secondFrequent2 = mostFrequent2;
            maxFrequence2 = content[i];
            mostFrequent2 = "ACGT-"[i];
        }
        else if (content[i] > secondFrequence2){
            secondFrequent2 = "ACGT-"[i];
            secondFrequence2 = content[i];
        }
    }

    for (auto i = 0 ; i < snp.size() ; i++){
        if (snp[i] == mostFrequent2){
            mostFrequentBases.push_back(1); 
        }
        else if (snp[i] == secondFrequent2){
            mostFrequentBases.push_back(-1);
        }
        else{
            mostFrequentBases.push_back(0);
        }
        lessFrequence.push_back(0);
        moreFrequence.push_back(1);
    }

    numberOfOccurences = 1;

}

//output : whether or not the position is significatively different from "all reads in the same haplotype"
//one exception possible : "all reads except last one in the same haplotype" -> that's because reads may be aligned on the last one, creating a bias
bool Partition::isInformative(float errorRate, bool lastReadBiased){

    int suspiciousReads [2] = {0,0}; //an array containing how many reads seriously deviate from the "all 1" or "all -1" theories

    auto adjust = 0;
    if (lastReadBiased){
        adjust = 1;
    }
    for (auto read = 0 ; read < mostFrequentBases.size()-adjust ; read++){

        int readNumber = moreFrequence[read] + lessFrequence[read];
        float threshold = 0.5*readNumber + 3*sqrt(readNumber*0.5*(1-0.5)); //to check if we deviate significantly from the "random read", that is half of the time in each partition
        threshold = min(threshold, float(readNumber)-1);
        
        if (moreFrequence[read] > threshold){
            if (mostFrequentBases[read] == -1){ //this deviates from the "all 1 theory"
                suspiciousReads[0]++;
            }
            else if (mostFrequentBases[read] == 1){ //this deviates from the "all -1 theory"
                suspiciousReads[1]++;
            }
        }

    }

    //if I have less than two suspicious reads, then this partition is not informative
    if (suspiciousReads[0] < 2 || suspiciousReads[1] < 2){
        return false;
    }
    else{
        return true;
    }

}

//input : a new partition to add to the consensus (the partition must be well-phased)
//output : updated consensus partition
//WARNING : the input partition must be of the same size as the consensus partition
void Partition::augmentPartition(vector<short>& supplementaryPartition){

    for (auto read = 0 ; read < mostFrequentBases.size() ; read++){
        if (supplementaryPartition[read] != 0){
            if (mostFrequentBases[read] == 0){ //the new partition gives some new position
                mostFrequentBases[read] = supplementaryPartition[read];
                moreFrequence[read] += 1;
            }
            else if (supplementaryPartition[read] == mostFrequentBases[read]){ //the new partition agrees
                moreFrequence[read] += 1;
            }
            else{ //the new partition disagrees
                lessFrequence[read] += 1;
                if (lessFrequence[read] > moreFrequence[read]){ //then the consensus has changed !
                    lessFrequence[read] -= 1;
                    mostFrequentBases[read] *= -1;
                    moreFrequence[read] += 1;
                }
            }
        }
    }
    numberOfOccurences += 1;
}

//input : another partition to be merged into this one and short phased (worth 1 or -1) (are the 0 in front of the 0 or the 1 ?)
//output : updated consensus partition
//WARNING like above, the two partitions must be the same size
void Partition::mergePartition(Partition p, short phased){

    auto moreOther = p.getMore();
    auto lessOther = p.getLess();

    auto other = p.getPartition();

    for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
        if (mostFrequentBases[c] == 0){
            mostFrequentBases[c] = other[c];
            moreFrequence[c] = moreOther[c];
            lessFrequence[c] = lessOther[c];
        }
        else if (phased*other[c] == mostFrequentBases[c]){ //the two partitions agree
            moreFrequence[c] += moreOther[c];
            lessFrequence[c] += lessOther[c];
        }
        else if (phased*other[c] != mostFrequentBases[c]){ //the two partitions disagree
            moreFrequence[c] += lessOther[c];
            lessFrequence[c] += moreOther[c];

            //then the most popular base may have changed
            if (lessFrequence[c] > moreFrequence[c]){
                mostFrequentBases[c] *= -1;
                auto stock = moreFrequence[c];
                moreFrequence[c] = lessFrequence[c];
                lessFrequence[c] = stock;
            }
        }
    }

    numberOfOccurences += p.number();
}

//input : another partition to be merged into this one
//output : updated consensus partition
//WARNING like above, the two partitions must be the same size
void Partition::mergePartition(Partition p){

    auto moreOther = p.getMore();
    auto lessOther = p.getLess();

    auto other = p.getPartition();

    //first determine the phase between the two partitions
    float phase = 0.01;
    for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
        phase += mostFrequentBases[c]*other[c]; //going toward + if phased, toward - if unphased
    }
    auto phased = phase / std::abs(phase);

    for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
        if (mostFrequentBases[c] == 0){
            mostFrequentBases[c] = other[c]*phased;
            moreFrequence[c] = moreOther[c];
            lessFrequence[c] = lessOther[c];
        }
        else if (phased*other[c] == mostFrequentBases[c]){ //the two partitions agree
            moreFrequence[c] += moreOther[c];
            lessFrequence[c] += lessOther[c];
        }
        else if (phased*other[c] == -mostFrequentBases[c]){ //the two partitions disagree
            moreFrequence[c] += lessOther[c];
            lessFrequence[c] += moreOther[c];

            //then the most popular base may have changed
            if (lessFrequence[c] > moreFrequence[c]){
                mostFrequentBases[c] *= -1;
                auto stock = moreFrequence[c];
                moreFrequence[c] = lessFrequence[c];
                lessFrequence[c] = stock;
            }
        }
    }

    numberOfOccurences += p.number();
}

//returns the majoritary partition
vector<short> Partition::getPartition(){
    return mostFrequentBases;
}

vector<float> Partition::getConfidence(){
    vector<float> conf;
    for (int i = 0 ; i < moreFrequence.size() ; i++){
        //conf.push_back(moreFrequence[i]+lessFrequence[i]);
        if (moreFrequence[i]+lessFrequence[i] > 0) {
            conf.push_back(float(moreFrequence[i])/(moreFrequence[i]+lessFrequence[i]));
        }
        else {
            conf.push_back(1);
        }
    }

    return conf;
}

vector<int> Partition::getMore(){
    return moreFrequence;
}

vector<int> Partition::getLess(){
    return lessFrequence;
}

int Partition::number(){
    return numberOfOccurences;
}

void Partition::print(){
    for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
        auto ch = mostFrequentBases[c];
        if (moreFrequence[c] == 0){
            cout << "o";
        }
        else if (float(moreFrequence[c])/(moreFrequence[c]+lessFrequence[c]) < 0.7){
            cout << "!";
        }
        else if (ch == 1){
            cout << 1;
        }
        else if (ch == -1){
            cout << 0;
        }
        else {
            cout << '?';
        }
        //cout << c << ",";
    }
    cout << " " << numberOfOccurences << endl;
    // for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
    //     cout << moreFrequence[c] << "/" << lessFrequence[c] << ";";
    // }
    // cout << endl;
}

int Partition::size(){
    return mostFrequentBases.size();
}

float Partition::proportionOf1(){
    float n1 = 0;
    float n0 = 0;
    for (auto r : mostFrequentBases){
        if (r == 1){
            n1++;
        }
        else if (r==-1){
            n0++;
        }
    }
    if (n1+n0 > 0){
        return n1/(n1+n0);
    }
    else{
        return 0.5;
    }
}
