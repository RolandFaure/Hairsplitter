#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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
            string name, seq, tags;
            string field;
            int i = 0;
            std::istringstream line2(line);
            while (std::getline(line2, field, '\t')){
                if (i == 1){
                    name = field;
                }
                if (i == 2){
                    seq = field;
                }
                if (i == 3){
                    tags = field;
                }
                if (i > 3){
                    tags += "\t" + field;
                }
                i++;
            }

            if (seq.size() == 0){
                continue;
            }

            cout << ">" << name << " " << tags << "\n";
            cout << seq << "\n";
        }
    }

}