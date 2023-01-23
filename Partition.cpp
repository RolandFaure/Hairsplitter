#include "Partition.h"

#include "robin_hood.h"
#include <iostream>
#include <cmath>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::list;
using std::vector;
using std::min;
using std::max;
using std::begin;
using std::end;
using std::pair;
using std::make_pair;

Partition::Partition(){
    readIdx = {};
    mostFrequentBases = {}; 
    moreFrequence = {};
    lessFrequence = {};
    numberOfOccurences = 0;
    conf_score = 0;
    pos_left = -1;
    pos_right = -1;
}

Partition::Partition(Column& snp, int pos, char ref_base){

    pos_left = pos;
    pos_right = pos;
    conf_score = 0;

    readIdx = vector<int>(snp.readIdxs.begin(), snp.readIdxs.end());

    robin_hood::unordered_flat_map<char, short> bases2content;
    bases2content['A'] = 0;
    bases2content['C'] = 1; 
    bases2content['G'] = 2;
    bases2content['T'] = 3;
    bases2content['-'] = 4;
    bases2content['?'] = 5;

    int content [5] = {0,0,0,0,0}; //item 0 for A, 1 for C, 2 for G, 3 for T, 4 for *, 5 for '-'
    for (int c = 0 ; c < snp.content.size() ; c++){
        if (snp.content[c] != '?'){
            content[bases2content[snp.content[c]]] += 1;
        }
    }

    /*
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
    */
    char mostFrequent = ref_base;
    auto maxFrequence = content[bases2content[ref_base]];
    //now find the second most frequent base
    char secondFrequent;
    int maxFrequence2 = -1;
    for (auto i = 0 ; i < 5 ; i++){
        if (ref_base != "ACGT-"[i]) {
            if (content[i] > maxFrequence2){
                secondFrequent = "ACGT-"[i];
                maxFrequence2 = content[i];
            }
        }
    }

    for (auto i = 0 ; i < snp.content.size() ; i++){
        if (snp.content[i] == mostFrequent){
            mostFrequentBases.push_back(1); 
        }
        else if (snp.content[i] == secondFrequent){
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
bool Partition::isInformative(bool lastReadBiased){

    int suspiciousReads [2] = {0,0}; //an array containing how many reads seriously deviate from the "all 1" or "all -1" theories

    auto adjust = 0;
    if (lastReadBiased){
        adjust = 1;
    }

    int numberOfReads = 0;
    for (auto read = 0 ; read < mostFrequentBases.size()-adjust ; read++){

        int readNumber = moreFrequence[read] + lessFrequence[read];
        float threshold = 0.5*readNumber + 3*sqrt(readNumber*0.5*(1-0.5)); //to check if we deviate significantly from the "random read", that is half of the time in each partition
        threshold = min(threshold, float(readNumber)-1);
        
        if (moreFrequence[read] > threshold){
            if (mostFrequentBases[read] == -1){ //this deviates from the "all 1 theory"
                suspiciousReads[0]++;
                numberOfReads++;
            }
            else if (mostFrequentBases[read] == 1){ //this deviates from the "all -1 theory"
                suspiciousReads[1]++;
                numberOfReads++;
            }
        }

    }

    //if I have less than two/10% suspicious reads, then this partition is not informative
    float minNumberOfReads = max(2.0, 0.1*numberOfReads);
    if (suspiciousReads[0] < minNumberOfReads || suspiciousReads[1] < minNumberOfReads){
        return false;
    }
    else{
        return true;
    }

}

//input : a new partition to add to the consensus (the partition must be well-phased in 'A' and 'a')
//output : updated consensus partition
void Partition::augmentPartition(Column& supplementaryPartition, int pos){

    //first adjust pos right and pos left if pos changes the limits of the partition
    
    if (pos != -1){
        if (pos < pos_left || pos_left == -1) {pos_left = pos;}
        if (pos > pos_right) {pos_right = pos;}
    }

    //determine the two most frequent bases in supplementaryPartition.content. There can be any character in supplementaryPartition.content
    vector<int> content (256, 0);//one item for each possible character
    for (auto c : supplementaryPartition.content){
        content[c]++;
    }
    //go through the content vector to find the two most frequent bases
    char mostFrequent = 'b';
    char secondFrequent = 'b';
    int maxFrequence = -1;
    //find the index of the two biggest elements of content
    for (auto i = 0 ; i < 256 ; i++){
        if (content[i] > maxFrequence){
            mostFrequent = i;
            maxFrequence = content[i];
        }
    }
    int maxFrequence2 = -1;
    for (auto i = 0 ; i < 256 ; i++){
        if (i != mostFrequent) {
            if (content[i] > maxFrequence2){
                secondFrequent = i;
                maxFrequence2 = content[i];
            }
        }
    }

    //print the two most frequent bases
    //check if the two most frequent bases should be swapped
    int n1 = 0;
    int n2 = 0;
    int swapped = 0;
    for (auto read : readIdx){
        while (supplementaryPartition.readIdxs[n2] < read){
            n2++;
        }
        if (supplementaryPartition.readIdxs[n2] == read){
            if (supplementaryPartition.content[n2] == mostFrequent && mostFrequentBases[n1] == 1){
                swapped += 1;
            }
            else if (supplementaryPartition.content[n2] == mostFrequent && mostFrequentBases[n1] == -1){
                swapped -= 1;
            }
            else if (supplementaryPartition.content[n2] == secondFrequent && mostFrequentBases[n1] == -1){
                swapped += 1;
            }
            else if (supplementaryPartition.content[n2] == secondFrequent && mostFrequentBases[n1] == 1){
                swapped -= 1;
            }
        }
        n1++;
    }
    if (swapped < 0){
        char tmp = mostFrequent;
        mostFrequent = secondFrequent;
        secondFrequent = tmp;
    }


    auto it1 = readIdx.begin();
    vector<int> idxs1_2;
    vector<short> mostFrequentBases_2;
    vector<int> moreFrequence_2;
    vector<int> lessFrequence_2;

    n1 = 0;
    n2 = 0;
    for (auto read : supplementaryPartition.readIdxs){
        while(it1 != readIdx.end() && *it1 < read){ //positions that existed in the old partitions that are not found here
            mostFrequentBases_2.push_back(mostFrequentBases[n1]);
            moreFrequence_2.push_back(moreFrequence[n1]);
            lessFrequence_2.push_back(lessFrequence[n1]);
            idxs1_2.push_back(*it1);
            it1++;
            n1++;
        }
        short s = 0;
        if (supplementaryPartition.content[n2] == secondFrequent){s=-1;}
        if (supplementaryPartition.content[n2] == mostFrequent){s=1;}
        if (it1 == readIdx.end() || *it1 != read){ //then this is a new read
            n1--; //because n1++ further down
            if (s == 1){
                mostFrequentBases_2.push_back(1);
                moreFrequence_2.push_back(1);
                lessFrequence_2.push_back(0);
                idxs1_2.push_back(read);
            }
            else if (s == -1){
                mostFrequentBases_2.push_back(-1);
                moreFrequence_2.push_back(1);
                lessFrequence_2.push_back(0);
                idxs1_2.push_back(read);
            }
        }
        else{ //we're looking at a read that is both old and new
            if (s == 0){ //the new partition does not brign anything
                mostFrequentBases_2.push_back(mostFrequentBases[n1]);
                moreFrequence_2.push_back(moreFrequence[n1]);
                lessFrequence_2.push_back(lessFrequence[n1]);
            }
            else if (mostFrequentBases[n1] == 0){ //no information on previous partition
                mostFrequentBases_2.push_back(s);
                moreFrequence_2.push_back(1);
                lessFrequence_2.push_back(0);
            }
            else if (s == mostFrequentBases[n1]){ //both partitions agree
                mostFrequentBases_2.push_back(mostFrequentBases[n1]);
                moreFrequence_2.push_back(moreFrequence[n1]+1);
                lessFrequence_2.push_back(lessFrequence[n1]);
            }
            else if (s == -mostFrequentBases[n1]) { //the new partitions disagrees
                if (lessFrequence[n1]+1 > moreFrequence[n1]){ //then the consensus has changed !
                    mostFrequentBases_2.push_back(-mostFrequentBases[n1]);
                    moreFrequence_2.push_back(moreFrequence[n1]+1);
                    lessFrequence_2.push_back(lessFrequence[n1]);
                }
                else{
                    mostFrequentBases_2.push_back(mostFrequentBases[n1]);
                    moreFrequence_2.push_back(moreFrequence[n1]);
                    lessFrequence_2.push_back(lessFrequence[n1]+1);
                }
            }
            idxs1_2.push_back(read);
            it1++;
        }
        n1++;
        n2++;
    }
    while(it1 != readIdx.end()){ //positions that existed in the old partitions that are not found in the new one
        mostFrequentBases_2.push_back(mostFrequentBases[n1]);
        moreFrequence_2.push_back(moreFrequence[n1]);
        lessFrequence_2.push_back(lessFrequence[n1]);
        idxs1_2.push_back(*it1);
        it1++;
        n1++;
    }

    mostFrequentBases = mostFrequentBases_2;
    moreFrequence = moreFrequence_2;
    lessFrequence = lessFrequence_2;
    readIdx = idxs1_2;
    numberOfOccurences += 1;
}

//input : another partition to be merged into this one and short phased (worth 1 or -1) (are the 0 in front of the 0 or the 1 ?)
//output : updated consensus partition
void Partition::mergePartition(Partition &p, short phased){

    //first adjust the limits of the partition
    pos_left = std::min(pos_left, p.get_left());
    pos_right = std::max(pos_right, p.get_right());

    //then merge the content of the two partitions
    auto moreOther = p.getMore();
    auto lessOther = p.getLess();

    auto other = p.getPartition();
    auto idx2 = p.getReads();

    vector<int> newIdx;
    vector<short> newMostFrequent;
    vector<int> newMoreFrequence;
    vector<int> newLessFrequence;

    int n1 = 0;
    int n2 = 0;
    while (n1 < readIdx.size() && n2 < idx2.size() ){
        while (readIdx[n1] < idx2[n2] && n1 < readIdx.size()){
            newIdx.push_back(readIdx[n1]);
            newMostFrequent.push_back(mostFrequentBases[n1]);
            newMoreFrequence.push_back(moreFrequence[n1]);
            newLessFrequence.push_back(lessFrequence[n1]);
            n1++;
        }

        while (readIdx[n1] > idx2[n2] && n2 < idx2.size()){
            newIdx.push_back(idx2[n2]);
            newMostFrequent.push_back(other[n2]*phased);
            newMoreFrequence.push_back(moreOther[n2]);
            newLessFrequence.push_back(lessOther[n2]);
            n2++;
        }

        if (n1 < readIdx.size() && n2 < idx2.size() && readIdx[n1] == idx2[n2]){ //the read is present on both partitions
            newIdx.push_back(readIdx[n1]);
            if (mostFrequentBases[n1] == 0){
                newMostFrequent.push_back(other[n2]*phased);
                newMoreFrequence.push_back(moreOther[n2]);
                newLessFrequence.push_back(lessOther[n2]);
            }
            else if (other[n2] == 0){
                newMostFrequent.push_back(mostFrequentBases[n1]);
                newMoreFrequence.push_back(moreFrequence[n1]);
                newLessFrequence.push_back(lessFrequence[n1]);
            }
            else if (phased*other[n2] == mostFrequentBases[n1]){ //the two partitions agree
                newMostFrequent.push_back(mostFrequentBases[n1]);
                newMoreFrequence.push_back(moreFrequence[n1]+moreOther[n2]);
                newLessFrequence.push_back(lessFrequence[n1]+lessOther[n2]);
            }
            else if (phased*other[n2] == -mostFrequentBases[n1]){ //the two partitions disagree
                newMostFrequent.push_back(mostFrequentBases[n1]);
                newMoreFrequence.push_back(moreFrequence[n1]+lessOther[n2]);
                newLessFrequence.push_back(lessFrequence[n1]+moreOther[n2]);

                //then the most popular base may have changed
                if (newLessFrequence[newLessFrequence.size()-1] > newMoreFrequence[newMoreFrequence.size()-1] ){
                    newMostFrequent[newMostFrequent.size()-1] *= -1;
                    auto stock = newMoreFrequence[newMoreFrequence.size()-1];
                    newMoreFrequence[newMoreFrequence.size()-1] = newLessFrequence[newLessFrequence.size()-1];
                    newLessFrequence[newLessFrequence.size()-1] = stock;
                }
            }
            n1++;
            n2++;
        }
    }
    while (n2 < idx2.size()){
        newIdx.push_back(idx2[n2]);
        newMostFrequent.push_back(other[n2]*phased);
        newMoreFrequence.push_back(moreOther[n2]);
        newLessFrequence.push_back(lessOther[n2]);
        n2++;
    }
    while (n1 < readIdx.size()){
        newIdx.push_back(readIdx[n1]);
        newMostFrequent.push_back(mostFrequentBases[n1]);
        newMoreFrequence.push_back(moreFrequence[n1]);
        newLessFrequence.push_back(lessFrequence[n1]);
        n1++;
    }

    readIdx = newIdx;
    mostFrequentBases = newMostFrequent;
    lessFrequence = newLessFrequence;
    moreFrequence = newMoreFrequence;

    numberOfOccurences += p.number();
}

//input : another partition to be merged into this one
//output : updated consensus partition
//WARNING like above, the two partitions must be the same size
void Partition::mergePartition(Partition &p){

    auto moreOther = p.getMore();
    auto lessOther = p.getLess();

    auto other = p.getPartition();
    auto idx2 = p.getReads();

    //first determine the phase between the two partitions
    int n1 = 0;
    int n2 = 0;
    float phase = 0.01; //not 0 to be sure phase is not 0 at the end (by default we choose +1 for phased)
    while (n1<readIdx.size() && n2 < idx2.size()){
        if (readIdx[n1] < idx2[n2]){
            n1++;
        }
        else if (readIdx[n1] > idx2[n2]){
            n2++;
        }
        else{
            phase += mostFrequentBases[n1]*other[n2];
            n1++;
            n2++;
        }
    }
    auto phased = phase / std::abs(phase);

    this->mergePartition(p, phased);
}

/**
 * @brief Extend the partition with the partition composed of all positions correlated with this partition aggregated
 * 
 * @param p 
 */
void Partition::extend_with_partition(Partition &p){
    numberOfOccurences = p.number();

    //iterate over the reads of the partition
    pos_left = std::min(pos_left, p.get_left());
    pos_right = std::max(pos_right, p.get_right());

    //then merge the content of the two partitions priorizing the current partition
    auto moreOther = p.getMore();
    auto lessOther = p.getLess();

    auto other = p.getPartition();
    auto idx2 = p.getReads();

    vector<int> newIdx;
    vector<short> newMostFrequent;
    vector<int> newMoreFrequence;
    vector<int> newLessFrequence;

    int n1 = 0;
    int n2 = 0;
    while (n1 < readIdx.size() && n2 < idx2.size() ){
        while (readIdx[n1] < idx2[n2] && n1 < readIdx.size()){
            newIdx.push_back(readIdx[n1]);
            newMostFrequent.push_back(mostFrequentBases[n1]);
            newMoreFrequence.push_back(moreFrequence[n1]);
            newLessFrequence.push_back(lessFrequence[n1]);
            n1++;
        }

        while (readIdx[n1] > idx2[n2] && n2 < idx2.size()){ //if read is present on the other partition but not on this one
            newIdx.push_back(idx2[n2]);
            newMostFrequent.push_back(other[n2]);
            newMoreFrequence.push_back(moreOther[n2]);
            newLessFrequence.push_back(lessOther[n2]);
            n2++;
        }

        if (n1 < readIdx.size() && n2 < idx2.size() && readIdx[n1] == idx2[n2]){ //the read is present on both partitions, priorize the current partition
            newIdx.push_back(readIdx[n1]);
            if (mostFrequentBases[n1] == 0){
                newMostFrequent.push_back(other[n2]);
                newMoreFrequence.push_back(moreOther[n2]);
                newLessFrequence.push_back(lessOther[n2]);
            }
            else{
                newMostFrequent.push_back(mostFrequentBases[n1]);
                newMoreFrequence.push_back(moreFrequence[n1]);
                newLessFrequence.push_back(lessFrequence[n1]);
            }
            n1++;
            n2++;
        }
    }
    while (n2 < idx2.size()){
        newIdx.push_back(idx2[n2]);
        newMostFrequent.push_back(other[n2]);
        newMoreFrequence.push_back(moreOther[n2]);
        newLessFrequence.push_back(lessOther[n2]);
        n2++;
    }
    while (n1 < readIdx.size()){
        newIdx.push_back(readIdx[n1]);
        newMostFrequent.push_back(mostFrequentBases[n1]);
        newMoreFrequence.push_back(moreFrequence[n1]);
        newLessFrequence.push_back(lessFrequence[n1]);
        n1++;
    }

    readIdx = newIdx;
    mostFrequentBases = newMostFrequent;
    lessFrequence = newLessFrequence;
    moreFrequence = newMoreFrequence;

}

/**
 * @brief Strengthen the partition with the partition p
 * 
 * @param p 
 */
void Partition::strengthen_partition(Partition &p){

    //check if the partitiin has the same phase
    int phase = 0;
    int n = 0;
    int n1 = 0;
    auto preads = p.getReads();
    auto ppartition = p.getPartition();
    for (auto n1 = 0 ; n1 < readIdx.size() ; n1++){
        while (preads[n] < readIdx[n1] && n < preads.size()){
            n++;
        }
        if (preads[n] == readIdx[n1]){
            if (ppartition[n] != mostFrequentBases[n1]){
                phase -= 1;
            }
            else {
                phase += 1;
            }
            n++;
        }
    }
    if (phase < 0){
        this->flipPartition();
    }

    n = 0;
    preads = p.getReads();
    auto pmore = p.getMore();
    auto pless = p.getLess();
    ppartition = p.getPartition();
    auto pconf = p.getConfidence();

    auto allConfidences = this->getConfidence();

    for (n1 = 0 ; n1 < readIdx.size() ; n1++){
        while (preads[n] < readIdx[n1] && n < preads.size()){
            n++;
        }
        if (preads[n] == readIdx[n1]){
            if ((allConfidences[n1] < 0.8 || numberOfOccurences < 5) && (pconf[n]>0.7 && p.number()>5)){
                mostFrequentBases[n1] = ppartition[n];
                moreFrequence[n1] = pmore[n];
                lessFrequence[n1] = pless[n];
            }
            n++;
        }
    }

}

//input : nothing except the partition itself
//output : a confidence score of the partition stored
float Partition::compute_conf(){
    double conf = 1;
    int numberReads = 0;
    auto confidences = this->getConfidence();
    for (auto c = 0 ; c < confidences.size() ; c++) {
        if (moreFrequence[c] > 1){
            conf*=confidences[c];
            numberReads++;
        }   
    }
    if (conf == 1){ //do not divide by 0 when calculating the score
        conf = 0.99;
    }
    conf_score = pow(1/(1-exp(log(conf)/numberReads)), 2)*this->number(); //exp(log(conf)/numberReads) is the geometrical average confidence

    return conf_score;
}

//returns the majoritary partition
vector<short> Partition::getPartition(){
    return mostFrequentBases;
}

/**
 * @brief Returns the partition, taking into account the fact that -1 and 1 are not equiprobable, so switching unsure reads towards the less confident base
 * 
 * @return vector<short> 
 */
vector<short> Partition::get_tweaked_partition(){

    vector<short> tweaked;

    //you could do a very smart probabilistic computation here if we knew how many partitions participated to each read, but we don't

    //compute the average confidence for the -1 and 1
    double conf1 = 0;
    double conf2 = 0;
    int nb1 = 0;
    int nb2 = 0;
    for (auto i = 0 ; i < mostFrequentBases.size() ; i++){
        if (mostFrequentBases[i] == 1){
            conf1+=float(moreFrequence[i])/(moreFrequence[i]+lessFrequence[i]);
            nb1++;
        }
        else if (mostFrequentBases[i] == -1){
            conf2+=float(moreFrequence[i])/(moreFrequence[i]+lessFrequence[i]);
            nb2++;
        }
    }
    if (nb1 > 0){
        conf1/=nb1;
    }
    if (nb2 > 0){
        conf2/=nb2;
    }

    //compute p1 such that (1-conf1)/(1-conf2) = p1/(1-p1)
    double ratio = (1-conf1)/(1-conf2);
    double p1 = ratio/(1+ratio);

    for (auto i = 0 ; i < mostFrequentBases.size() ; i++){
        if (mostFrequentBases[i] == 0){
            tweaked.push_back(0);
        }
        else if (mostFrequentBases[i] == 1){
            if (float(moreFrequence[i])/(moreFrequence[i]+lessFrequence[i]) > 1-p1){
                tweaked.push_back(1);
            }
            else {
                tweaked.push_back(-1);
            }
        }
        else if (mostFrequentBases[i] == -1){
            if (float(moreFrequence[i])/(moreFrequence[i]+lessFrequence[i]) > p1){
                tweaked.push_back(-1);
            }
            else {
                tweaked.push_back(1);
            }
        }
    }

    cout << "ufqoidp: " << p1 << " " << conf1 << " " << conf2 << endl;
    cout << "ufqoidp here is the tweaked partition : " << endl;
    for (auto i = 0 ; i < tweaked.size() ; i++){
        cout << tweaked[i] << " ";
    }

    return tweaked;
}

vector<int> Partition::getReads(){
    return readIdx;
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

    int c = 0;
    int n = 0;
    auto it = readIdx.begin();


    while(it != readIdx.end()){
        if (*it > c){
            cout << "?";
        }
        else{
            auto ch = mostFrequentBases[n];
            if (moreFrequence[n] == 0){
                cout << "o";
            }
            else if (float(moreFrequence[n])/(moreFrequence[n]+lessFrequence[n]) < 0.6){
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
            it++;
            n++;
        }
        c++;
    }

    cout << " " << numberOfOccurences << " " << pos_left << " <-> " << pos_right << endl;
    //compute the average confidence of positions containing -1 and of positions containing 1
    // float conf1 = 0;
    // float conf0 = 0;
    // int n1 = 0;
    // int n0 = 0;
    // for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
    //     if (mostFrequentBases[c] == 1){
    //         conf1 += float(moreFrequence[c])/(moreFrequence[c]+lessFrequence[c]);
    //         n1++;
    //     }
    //     else if (mostFrequentBases[c] == -1){
    //         conf0 += float(moreFrequence[c])/(moreFrequence[c]+lessFrequence[c]);
    //         n0++;
    //     }
    // }
    // if (n1 > 0){
    //     conf1 /= n1;
    // }
    // if (n0 > 0){
    //     conf0 /= n0;
    // }
    // cout << "conf 1 : " << conf1 << " conf 0 : " << conf0 << endl;

    // for (auto c = 0 ; c < mostFrequentBases.size() ; c++){
    //     cout << moreFrequence[c] << "/" << lessFrequence[c] << ";";
    // }
    // cout << endl;
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

float Partition::get_conf(){
    if (conf_score == 0){
        this->compute_conf();
    }
    return conf_score;
}

int Partition::get_left(){
    return pos_left;
}

int Partition::get_right(){
    return pos_right;
}

//transforms all -1 in 1 and vice-versa
void Partition::flipPartition(){
    for (auto r = 0 ; r < mostFrequentBases.size() ; r++){
        mostFrequentBases[r] *= -1;
    }
}

/**
 * @brief Mark all the positions that were excluded by the mask as -2s in the partition
 * 
 * @param mask 
 */
void Partition::apply_mask(std::vector<bool> &mask){
    //set the partition to 1 for all positions where the mask is false

    vector<int> newIdxReads;
    vector<short> newMostFrequent;
    vector<int> newMore;
    vector<int> newLess;

    int indexIdx = 0;
    for (auto r = 0 ; r < mask.size() ; r++){

        if (indexIdx < readIdx.size() && readIdx[indexIdx] == r){
            newIdxReads.push_back(r);
            newMostFrequent.push_back(mostFrequentBases[indexIdx]);
            newMore.push_back(moreFrequence[indexIdx]);
            newLess.push_back(lessFrequence[indexIdx]);
            indexIdx++;
        }
        else if (!mask[r]){
            newIdxReads.push_back(r);
            newMostFrequent.push_back(-2);
            newMore.push_back(1);
            newLess.push_back(0);
        }
    }

    readIdx = newIdxReads;
    mostFrequentBases = newMostFrequent;
    moreFrequence = newMore;
    lessFrequence = newLess;
}

/**
 * @brief Changes the partition of reads: useful after making a correction to the partition
 * 
 * @param newPartition 
 */
void Partition::new_corrected_partition(std::vector<short> newPartition, std::vector<int> newIdxs, std::vector<int> more, std::vector<int> less){
    readIdx = newIdxs;
    mostFrequentBases = newPartition;
    moreFrequence = more;
    lessFrequence = less;
}
