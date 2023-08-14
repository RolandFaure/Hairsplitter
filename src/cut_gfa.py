'''
Function to cut the long contigs of a gfa in smaller contigs
'''

import argparse
import numpy as np

def parse_args():
    parser=argparse.ArgumentParser(description="cuts the long contigs of a gfa in smaller contigs")
    parser.add_argument("--assembly", "-a", required=True, default=None, help="GFA assembly file")
    parser.add_argument("--length", "-l", required=True, help="Maximal length of the outputted contigs")
    parser.add_argument("--output", "-o", required=True, help="Output file")
    return parser.parse_args()

if __name__ == "__main__":
    
    args = parse_args()
    assembly = args.assembly
    length = int(args.length)
    output = args.output

    #dict associating each contig with its length
    length_of_contigs = {}

    #list all the L lines
    L_lines = []

    #open the output file
    with open(output, "w") as out:
        with open(assembly, "r") as f:
            #first go through all the lines starting with a S
            for line in f:
                if line.startswith("S"):
                    
                    ls = line.strip().split("\t")
                    length_of_contigs[ls[1]] = len(ls[2])
                    for chunk in range (int(np.floor((len(ls[2])-1)/length)+1)): #-1 to avoid creating a contig of length 0
                        if chunk*length < len(ls[2]):
                            # print("segment of length ", chunk*length, " to ", len(ls[2]), " of contig ", ls[2][chunk*length:len(ls[2])+1])
                            out.write("S\t" + ls[1] + "@" + str(chunk) + "\t" + ls[2][chunk*length:min((chunk+1)*length, len(ls[2])+1)] + "\t" + "\t".join(ls[3:]).strip("\n")+"\n")

                            
                            #link the slice with the previous one
                            if chunk > 0 :
                                out.write("L\t" + ls[1] + "@" + str(chunk-1) + "\t+\t" + ls[1] + "@" + str(chunk) + "\t+\t0M\n")
                
                elif line.startswith("L"):
                    L_lines.append(line)
                
        
            for line in L_lines:
                ls = line.strip("\n").split("\t")
                if ls[2] == "+":
                    out.write("L\t" + ls[1] + "@" + str(int(np.floor((length_of_contigs[ls[1]]-1)/length))) + "\t+\t")
                else:
                    out.write("L\t" + ls[1] + "@0\t-\t")

                if ls[4] == "-":
                    out.write(ls[3] + "@" + str(int(np.floor((length_of_contigs[ls[3]]-1)/length))) + "\t-\t"+ls[5])
                else:
                    out.write(ls[3] + "@0\t+\t"+ls[5])
                out.write("\n")

    out.close()
