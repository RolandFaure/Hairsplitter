#include "tools.h"

#include <iostream>

using std::cout;
using std::endl;
using std::stoi;
using std::string;
using std::vector;

//input : a CIGAR
//output : an alignment string where 1 letter = 1 base. i.e. 5M1D1M -> MMMMMIM
std::string convert_cigar(std::string &cigar){

    string res;
    string num = "";
    for (auto c : cigar){
        if ((int)c - '0' >= 0 && (int)c - '0' <= 9){
            num += c;
        }
        else{
            int n = std::stoi(num);
            for (auto i = 0 ; i < n ; i++){
                res += c;
            }
            num = "";
        }
    }

    return res;
}

string reverse_complement(string &seq){
    string res = "";
    for (auto c : seq){
        if (c == 'A'){
            res = 'T' + res;
        }
        else if (c == 'T'){
            res = 'A' + res;
        }
        else if (c == 'C'){
            res = 'G' + res;
        }
        else if (c == 'G'){
            res = 'C' + res;
        }
        else{
            res = c + res;
        }
    }
    return res;
}

void print_alignment(std::string &ref, std::string &read, std::string &cigar, int start, int end){

    string ref_aligned;
    string read_aligned;

    int ref_pos = start;
    int read_pos = 0;

    for (auto i = 0 ; i < start ; i++){
        ref_aligned += '-';
    }

    for (auto i = 0 ; i < cigar.size() ; i++){
        if (cigar[i] == 'M' || cigar[i] == 'X' || cigar[i] == '='){
            ref_aligned += ref[ref_pos];
            read_aligned += read[read_pos];
            ref_pos++;
            read_pos++;
        }
        else if (cigar[i] == 'D'){
            ref_aligned += '-';
            read_aligned += read[read_pos];
            read_pos++;
        }
        else if (cigar[i] == 'I'){
            ref_aligned += ref[ref_pos];
            read_aligned += '-';
            ref_pos++;
        }
    }

    // for (auto i = ref_pos ; i < end-1 ; i++){
    //     ref_aligned += '-';
    // }

    cout << ref_aligned << endl;
    cout << read_aligned << endl;

}

int compute_edit_distance(std::string &cigar, std::string &ref, std::string &read, int start, int end){
    
        int res = 0;
    
        int ref_pos = start;
        int read_pos = 0;
    
        for (auto i = 0 ; i < cigar.size() ; i++){
            if (cigar[i] == 'M' || cigar[i] == 'X' || cigar[i] == '='){
                if (ref[ref_pos] != read[read_pos]){
                    res++;
                }
                ref_pos++;
                read_pos++;
            }
            else if (cigar[i] == 'D'){
                res++;
                read_pos++;
            }
            else if (cigar[i] == 'I'){
                res++;
                ref_pos++;
            }
        }
    
        return res;
}

