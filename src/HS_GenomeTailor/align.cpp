#include "align.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::to_string;
using std::stoi;
using std::stof;
using std::stod;
using std::pair;
using std::make_pair;
using std::min;
using std::max;

string align(string &s1, size_t start1, size_t end1, string &s2, size_t start2, size_t end2){

    int index1 = start1;
    int index2 = start2;
    int number_of_consecutive_matches = 0;

    int num_match = 0;
    int num_indel = 0;

    int lost_path_1 = 0;
    int lost_path_2 = 0;

    int length_of_boaring = 1;
    int num_indel_in_boaring = 0;
    bool now_insert = false; //if true first insert, if false first delete

    string cigar = "";

    while (index1 < end1 && index2 < end2){
        if (s1[index1] == s2[index2]){
            index1++;
            index2++;
            number_of_consecutive_matches++;
            if (number_of_consecutive_matches > 3){
                num_match++;
                cigar += "=";
            }
            else if (number_of_consecutive_matches == 3){ //entering a new stretch
                // cout << "entering a new stretch at index " << index1 << " " << index2 << endl;
                // cout << "lost path is " << lost_path_1 << " " << lost_path_2 << endl;
                //link to lost index
                int local_index1, local_index2;
                if (lost_path_1 - index1 > lost_path_2 -index2){
                    num_indel += (lost_path_1 - index1) - (lost_path_2 -index2);
                    for (auto i = 0; i < (lost_path_1 - index1) - (lost_path_2 -index2); i++){
                        cigar += "I";
                    }
                    local_index1 = lost_path_1;
                    local_index2 = index2 + (lost_path_1 - index1);
                }
                else{
                    num_indel += (lost_path_2 - index2) - (lost_path_1 -index1);
                    for (auto i = 0; i < (lost_path_2 - index2) - (lost_path_1 -index1); i++){
                        cigar += "D";
                    }
                    local_index1 = index1 + (lost_path_2 - index2);
                    local_index2 = lost_path_2;
                }
                while(local_index1 < index1 && local_index2 < index2){
                    if (s1[local_index1] == s2[local_index2]){
                        num_match++;
                        cigar += "=";
                    }
                    else{
                        num_indel++;
                        cigar += "X";
                    }
                    local_index1++;
                    local_index2++;
                }

                int length_of_boaring = 1;
                int num_indel_in_boaring = 0;
                bool now_insert = false;
            }
        }
        else{
            if (number_of_consecutive_matches >= 3){
                lost_path_1 = index1;
                lost_path_2 = index2;
            }

            //now boar through the matrix
            if (now_insert){
                index2++;
                num_indel_in_boaring++;
            }
            else{
                index1++;
                num_indel_in_boaring++;
            }
            if (num_indel_in_boaring >= length_of_boaring){
                now_insert = !now_insert;
                num_indel_in_boaring = 0;
                length_of_boaring*=2;
            }            
            number_of_consecutive_matches = 0;
        }
    }

    num_indel += (end1 - index1) + (end2 - index2);
    for (auto i = index1; i < end1; i++){
        cigar += "D";
    }
    for (auto i = index2; i < end2; i++){
        cigar += "I";
    }

    return cigar;
}

//3M5I2D -> MMMIIIIIDD
std::string convert_CIGAR(std::string &cigar){
    std::string newCIGAR = "";
    int pos = 0;
    int num = 0;
    while (pos < cigar.size()){
        if (cigar[pos] >= '0' && cigar[pos] <= '9'){
            while (cigar[pos] >= '0' && cigar[pos] <= '9'){
                num = num * 10 + cigar[pos] - '0';
                pos++;
            }
        } else {
            for (int i = 0; i < num; i++){
                newCIGAR += cigar[pos];
            }
            num = 0;
            pos++;
        }
    }
    return newCIGAR;
}

/**
 * @brief Polish a sequence with reads 
 * @details Will not polish the edges of the sequence
 * 
 * @param toPolish 
 * @param reads 
 * @return std::string 
 */
std::string polish(std::string &toPolish, std::vector<std::string> &reads, std::string& path_minimap2, std::string& path_racon, std::string tmp_folder){

    //start by running minimap2 by aligning the reads to the sequence

    //create a temporary file to store the reads
    string readsFile = tmp_folder +"reads.tmp_polish.fa";
    ofstream readsStream(readsFile);
    for (int i = 0; i < reads.size(); i++){
        readsStream << ">read" << i << endl;
        readsStream << reads[i] << endl;
    }
    readsStream.close();

    //create a temporary file to store the seq to polish
    string seqFile = tmp_folder +"seq.tmp_polish.fa";
    ofstream seqStream(seqFile);
    seqStream << ">seq" << endl;
    seqStream << toPolish << endl;
    seqStream.close();

    //create a temporary file to store the alignment
    string alignmentFile = tmp_folder +"alignment.tmp_polish.paf";

    //run minimap2
    string command = path_minimap2 + " -x map-ont -t 1 " + seqFile + " " + readsFile + " > " + alignmentFile + " 2> "+ tmp_folder +"logminimap.tmp.txt";
    system(command.c_str());

    //run racon
    string polished = tmp_folder +"polished.tmp_polish.fa";
    string command2 = path_racon + " " + readsFile + " " + alignmentFile + " " + seqFile + " > " + polished + " 2> "+ tmp_folder +"logracon.tmp.txt";
    int racon_success = system(command2.c_str());

    string polishedSeq = toPolish;
    if (racon_success == 0){
        //read the polished sequence
        ifstream polishedStream(polished);
        string header;
        getline (polishedStream, header);
        getline(polishedStream, polishedSeq);
        polishedStream.close();
    }

    //delete the temporary files
    system(("rm "+ tmp_folder +"*.tmp_polish.*").c_str());

    return polishedSeq;
}