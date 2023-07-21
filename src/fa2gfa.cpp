#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{   
    //in a streaming fashion, convert a fasta file to a gfa file
    //the fasta file is given as an argument
    //the gfa file is written to stdout

    std::string fastafile = argv[1];
    std::ifstream infile(fastafile);
    std::string line;
    std::string name;
    std::string sequence;
    int i = 0;
    while (std::getline(infile, line))
    {
        if (line[0] == '>')
        {
            if (i > 0)
            {
                std::cout << "S\t" << name << "\t" << sequence << "\n";
            }
            name = line.substr(1);
            sequence = "";
            i++;
        }
        else
        {
            sequence += line;
        }
    }
    if (i > 0){
        std::cout << "S\t" << name << "\t" << sequence << "\n";
    }
    return 0;
}
