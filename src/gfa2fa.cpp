#include <iostream>
#include <fstream>
#include <string>

using std::cout;
using std::string;

int main(int argc, char *argv[])
{  
    //in a streaming fashion, convert a gfa file to a fasta file
    //the gfa file is given as an argument
    //the fasta file is written to stdout

    std::string gfafile = argv[1];
    std::ifstream infile(gfafile);
    std::string line;
    std::string name;
    std::string sequence;
    int i = 0;
    while (std::getline(infile, line))
    {
        if (line[0] == 'S'){

            //the name is the second field
            string name = line.substr(2, line.find('\t', 2) - 2);
            name = name.substr(0, name.find(' '));

            //the sequence is the third field
            string seq = line.substr(line.find('\t', 2) + 1, line.find('\t', line.find('\t', 2) + 1) - line.find('\t', 2) - 1);
            if (seq.size() == 0){
                continue;
            }
            cout << ">" << name << "\n";
            cout << seq << "\n";
        }
    }

}