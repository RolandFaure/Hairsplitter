#include "Partition.h"

#include "robin_hood.h"
#include <iostream>

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
    numberOfOccurences = 1;
}

Partition::Partition(vector<char> & snp){

    robin_hood::unordered_flat_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['*'] = 4;
    bases2content['-'] = 5;

    int content [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for *, 5 for '-'
    for (int c = 0 ; c < snp.size() ; c++){
        if (snp[c] != '-'){
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
            mostFrequent2 = "ACGT*"[i];
        }
        else if (content[i] > secondFrequence2){
            secondFrequent2 = "ACGT*"[i];
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

//input : a new partition to add to the consensus
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

//returns the majoritary partition
vector<short> Partition::getPartition(){
    return mostFrequentBases;
}

vector<float> Partition::getConfidence(){
    vector<float> conf;
    for (int i = 0 ; i < moreFrequence.size() ; i++){
        conf.push_back(moreFrequence[i]+lessFrequence[i]);
        //conf.push_back(float(moreFrequence[i]-lessFrequence[i])/(moreFrequence[i]+lessFrequence[i]));
    }

    return conf;
}

int Partition::number(){
    return numberOfOccurences;
}

void Partition::print(){
    for (auto c : mostFrequentBases){
        if (c == 1){
            cout << 1;
        }
        else if (c == -1){
            cout << 0;
        }
        else {
            cout << '*';
        }
        //cout << c << ",";
    }
    cout << endl;
}

int Partition::size(){
    return mostFrequentBases.size();
}

