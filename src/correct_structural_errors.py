#/usr/bin python3
#-*- coding:utf-8 -*-
#author = Roland Faure

"""
This module is used to correct structural errors in the gfa file based on the read following a simple rule: every read must align end-to-end on the assembly graph.
"""

import argparse
import numpy as np
import os
import re
import sys

resolution = 100

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return "".join(complement[base] for base in reversed(seq))

def parse_args():
    #input arguments: assembly graph, reads file, output file, and optionnaly a path to minigraph
    parser=argparse.ArgumentParser(description="corrects structural errors in the gfa file based on the read following a simple rule: every read must align end-to-end on the assembly graph")
    parser.add_argument("--assembly", "-a", required=True, default=None, help="GFA assembly file")
    parser.add_argument("--reads", "-r", required=True, default=None, help="Reads file")
    parser.add_argument("--output", "-o", required=True, help="Output file")
    parser.add_argument("--threads", "-t", required=False, default=1, help="Number of threads to use")
    parser.add_argument("--folder", required=False, default="tmp_correct_structural_errors", help="Folder where to store the temporary files")
    parser.add_argument("--minigraph", "-m", required=False, default="minigraph", help="Path to minigraph")
    parser.add_argument("--racon", "-c", required=False, default="racon", help="Path to racon")
    parser.add_argument("--minimap2", "-p", required=False, default="minimap2", help="Path to minimap2")

    return parser.parse_args()

#creating a class contig, storing the contigs and their links
class Contig:
    def __init__(self, name, sequence):
        self.name = name
        self.length = len(sequence)
        self.sequence = sequence

        self.neighbors = set() #set of tuple (position on this contig, orientation on this contig (direction of the unknown), other contig, position on the other contig, orientation on the other contig, CIGAR)
        self.potential_new_links = set() #just like neighbors, but not definitvie. (position on this contig, orientation on this contig (direction of the unknown), other contig, position on the other contig, orientation on the other contig, nameOfRead, orientationOnRead, posStartOnRead, posEndOnRead)
        self.future_links = {} #will be a collection of nearly identical potential_new_links

    #go through the potential new links and select those that are present more than twice
    def find_new_links(self):
        
        # if self.name !=  "edge_49" :
        #     return
        #go through the links and see if the link has not been seen before
        links = {} #associating a link with its number of occurences
        for link in self.potential_new_links:
            
            # print("looking at link ", link)
            #go through the already compyted link and compute a distance to see if it is the same
            added = False
            for already_link in links.keys() :
                # if link[2] == "edge_30":
                #     print("Comparing ", already_link, " and ", link)
                if (link[2] == "*" and already_link[2] == "*") and abs(link[0]-already_link[0]) < resolution*2 and link[1] == already_link[1] :
                    links[already_link].append(link[-4:])
                    added = True
                    break
                elif link[2] != "*" and abs(link[0]-already_link[0]) < resolution*2 and link[1] == already_link[1] and link[2]==already_link[2] and abs(link[3]-already_link[3]) < resolution*2 and link[4] == already_link[4] :  
                    links[already_link].append(link[-4:])
                    added = True
                    break

            if not added :
                links[link] = [link[-4:]]
               
        #now go through the links and see if they are present more than twice
        # print("All the links for contig ", self.name, " are ")
        for k in links.keys() :
            if len(links[k]) > 2 :
                # print(k)
                self.future_links[k] = links[k]
                # print(k, " ", links[k])
        # print(self.future_links)
        # if (len(links) > 0) :
        #     sys.exit(1)
       
#function that takes as input the new link and the contigs, and generates the new link
def generate_new_link(contigName, new_link, reads, contigs, read_file, read_index):

    # print("Generating new link ", contigName, " ", new_link)
    #inventoriate the reads and output them in a file
    reads_file = "tmp_reads.fa"
    toPolish_file = "tmp_toPolish.fa"
    fi = open(read_file, "r")
    fo = open(reads_file, "w")
    foundToPolish = False

    # if contigName == "edge_4" and new_link[2] == "edge_4" :
    #     print("oozzozozozoz ", reads)
    #     sys.exit()

    beginning_toPolish = -1
    end_toPolish = -1
    nreads = len(reads)
    nread = 0
    for read in reads:
        read_name = read[0]
        orientation_on_read = read[1]
        start = int(read[2])
        end = int(read[3])

        seq = ""
        fi.seek(read_index[read_name])
        seq = fi.readline().strip()
        fi.seek(0)

        if orientation_on_read == "-" :
            seq = reverse_complement(seq)
            start, end = len(seq)-end, len(seq)-start

        seq = seq[max(0,start-300):min(end+300, len(seq))] #take a little margin on both sides
        if (foundToPolish == False and start > 300 and end < len(seq)-300) or (foundToPolish == False and nread == nreads-1):
            foundToPolish = True
            beginning_toPolish = min(300, start)
            end_toPolish = max(len(seq)-300, len(seq)-end)
            foPolish = open(toPolish_file, "w")
            foPolish.write(">" + read_name + "\n" + seq + "\n")
            foPolish.close()
        fo.write(">" + read_name + "\n" + seq + "\n")
        nread += 1

    fo.close()

    #polish the joint using racon
    command = "minimap2 -x map-pb tmp_toPolish.fa tmp_reads.fa > tmp.paf 2> trash.txt"
    minimap = os.system(command)
    if minimap != 0 :
        print("Error while running minimap2: " + command + "\n")
        sys.exit(1)

    command = "racon tmp_reads.fa tmp.paf tmp_toPolish.fa > tmp_repolished.fa 2>trash.txt"
    racon = os.system(command)
    if racon != 0 :
        print("Error while running racon: " + command + "\n")
        sys.exit(1)

    #see where it aligns left and right of the junction
    #create a file with the sequence of the contigs left and right of the junction
    start_left = int(max(0, new_link[0]-300))
    end_left = int(min(new_link[0]+300, contigs[contigName].length))
    # print("start_left ", start_left, " end_left ", end_left, " length ", contigs[contigName].length)
    seq_left = contigs[contigName].sequence[start_left:end_left] #take a little margin on both sides
    if new_link[1] == "-" :
        seq_left = reverse_complement(seq_left) #reverse complement the sequence -> the interesting coordinate will always be the rightmost

    start_right = int(max(0, new_link[3]-300))
    end_right = int(min(new_link[3]+300, contigs[new_link[2]].length))
    # print("start_right ", start_right, " end_right ", end_right, " length ", contigs[new_link[2]].length)
    seq_right = contigs[new_link[2]].sequence[start_right:end_right] #take a little margin on both sides
    if new_link[4] == "+" :
        seq_right = reverse_complement(seq_right) #reverse complement the sequence -> the interesting coordinate will always be the leftmost

    #output left and right to a file and align them on the polished sequence
    left_and_right = "tmp_left_and_right.fa"
    fo = open(left_and_right, "w")
    fo.write(">left\n" + seq_left + "\n")
    fo.write(">right\n" + seq_right + "\n")
    fo.close()

    command = "minimap2 -cx asm20 tmp_repolished.fa tmp_left_and_right.fa > tmp_left_and_right.paf 2> trash.txt"
    minimap = os.system(command)
    if minimap != 0 :
        print("Error while running minimap2: " + command + "\n")
        sys.exit(1)

    #parse the paf file and get the CIGARs
    orientation_repolished = ">"
    #the four followning values should be adjusted, but in case there are problems during alignment, give them a default value
    pos_on_read_left = beginning_toPolish
    pos_on_read_right = end_toPolish
    pos_on_contig_left = new_link[0]
    pos_on_contig_right = new_link[3]
    with open("tmp_left_and_right.paf", "r") as f:
        for line in f:
            if line[0] != "@" :
                ls = line.strip().split("\t")
                if ls[0] == "left" :
                    if ls[4] == "-":
                        orientation_repolished = "<"

                    if new_link[1] == "+":
                        pos_on_contig_left = int(ls[3]) + start_left
                    else :
                        pos_on_contig_left = int(ls[1]) - int(ls[3]) + start_left

                    if ls[4] == "+" :
                        pos_on_read_left = int(ls[8])
                    else :
                        pos_on_read_left = int(ls[7])
                    
                elif ls[0] == "right" :
                    if new_link[4] == "+":
                        pos_on_contig_right = int(ls[1]) - int(ls[2]) + start_right
                    else :
                        pos_on_contig_right = int(ls[2]) + start_right
                    if ls[4] == "+" :
                        pos_on_read_right = int(ls[7])
                    else :
                        pos_on_read_right = int(ls[8])

    # print("The new link will go from pos ", pos_on_contig_left, " (",contigs[contigName].length,") of ", contigName, " to pos ", pos_on_contig_right, " of ", new_link[2], " (",contigs[new_link[2]].length,") and between that stands the polished sequence coord ", pos_on_read_left, " to ", pos_on_read_right)

    #to avoid extremely small hanging contigs
    if new_link[1] == "-" and pos_on_contig_left < 20 :
        pos_on_contig_left = 0 
    elif new_link[1] == "+" and pos_on_contig_left > contigs[contigName].length - 20 :
        pos_on_contig_left = contigs[contigName].length
    if new_link[4] == "-" and pos_on_contig_right < 20 :
        pos_on_contig_right = 0
    elif new_link[4] == "+" and pos_on_contig_right > contigs[new_link[2]].length - 20 :
        pos_on_contig_right = contigs[new_link[2]].length

    if pos_on_read_right > pos_on_read_left : #then retrieve the sequence
        seq_repolished = ""
        with open("tmp_repolished.fa", "r") as f:
            for line in f:
                if line[0] != ">" :
                    seq_repolished += line.strip()
        seq_repolished = seq[pos_on_read_left:pos_on_read_right]
        nameOfNewContig = "join_" + contigName +"_"+ str(pos_on_contig_left) + "_" + new_link[2]+"_"+str(pos_on_contig_right)
        new_contig = Contig(nameOfNewContig, seq_repolished)
        if orientation_repolished == ">" :
            new_contig.neighbors.add((0, "-", contigName, pos_on_contig_left, new_link[1], "0M"))
            contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], nameOfNewContig, 0, "-", "0M"))

            new_contig.neighbors.add((len(seq_repolished), "+", new_link[2], pos_on_contig_right, new_link[4], "0M"))
            contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4], nameOfNewContig, len(seq_repolished), "+", "0M"))

        else :
            new_contig.neighbors.add((len(seq_repolished), "-", contigName, pos_on_contig_left, new_link[1], "0M"))
            contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], nameOfNewContig, len(seq_repolished), "-", "0M"))

            new_contig.neighbors.add((0, "+", new_link[2], pos_on_contig_right, new_link[4], "0M"))
            contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4],nameOfNewContig, 0, "+", "0M"))

        contigs[nameOfNewContig] = new_contig

        # print("New contig ", new_contig.name, " has length ", new_contig.length, " and neighbors ", new_contig.neighbors)

    else : #means that the two contigs actually overlap, just add a link between the two
        length_of_CIGAR = pos_on_read_left - pos_on_read_right
        contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], new_link[2], pos_on_contig_right, new_link[4], str(length_of_CIGAR) + "M"))
        contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4], contigName, pos_on_contig_left, new_link[1], str(length_of_CIGAR) + "M"))
    

    #this remains to be tested
    # if contigName == "edge_33" and new_link[2] == "edge_36" and new_link[0] > 342000 :    
    #     print("ooooop")
    #     sys.exit(1)

#function that takes as input the new link towards the unknown and generated a new contig
def create_new_contig(contigName, new_link, reads, contigs, read_file, read_index):

    #begin by generating the new contig
    # print("Generating new contig ", contigName, " ", new_link)
    #inventoriate the reads and output them in a file
    reads_file = "tmp_reads.fa"
    toPolish_file = "tmp_toPolish.fa"
    fi = open(read_file, "r")
    fo = open(reads_file, "w")

    beginning_toPolish = -1
    end_toPolish = -1
    nreads = len(reads)
    nread = 0
    best_read = reads[0]
    best_read_sequence = ""
    best_read_length= 0
    for read in reads:
        read_name = read[0]
        orientation_on_read = read[1]
        start = int(read[2])
        end = int(read[3])

        seq = ""
        fi.seek(read_index[read_name])
        seq = fi.readline().strip()
        fi.seek(0)

        if orientation_on_read == "-" :
            seq = reverse_complement(seq)

        if (end-start)*(len(seq)-start)>best_read_length :
            best_read = read
            best_read_sequence = seq[max(0,start-300):min(end+300, len(seq))]
            best_read_length = (end-start)*(len(seq)-start)

        seq = seq[max(0,start-300):min(end+300, len(seq))] #take a little margin on both sides
        
        fo.write(">" + read_name + "\n" + seq + "\n")
        nread += 1

    fo.close()

    #write the best read to polish
    foPolish = open(toPolish_file, "w")
    foPolish.write(">" + best_read[0] + "\n" + best_read_sequence + "\n")
    foPolish.close()

    #polish the joint using racon
    command = "minimap2 -x map-pb tmp_toPolish.fa tmp_reads.fa > tmp.paf 2> trash.txt"
    minimap = os.system(command)
    if minimap != 0 :
        print("Error while running minimap2: " + command + "\n")
        sys.exit(1)

    command = "racon tmp_reads.fa tmp.paf tmp_toPolish.fa > tmp_repolished.fa 2>trash.txt"
    racon = os.system(command)
    if racon != 0 :
        print("Error while running racon: " + command + "\n")
        sys.exit(1)

    #see where it aligns left and right of the junction
    #create a file with the sequence of the contigs left and right of the junction
    start_left = int(max(0, new_link[0]-300))
    end_left = int(min(new_link[0]+300, contigs[contigName].length))
    # print("start_left ", start_left, " end_left ", end_left, " length ", contigs[contigName].length)
    seq_left = contigs[contigName].sequence[start_left:end_left] #take a little margin on both sides
    if new_link[1] == "-" :
        seq_left = reverse_complement(seq_left) #reverse complement the sequence -> the interesting coordinate will always be the rightmost

    #output left and right to a file and align them on the polished sequence
    left_and_right = "tmp_left_and_right.fa"
    fo = open(left_and_right, "w")
    fo.write(">left\n" + seq_left + "\n")
    fo.close()

    command = "minimap2 -cx asm20 tmp_repolished.fa tmp_left_and_right.fa > tmp_left_and_right.paf 2> trash.txt"
    minimap = os.system(command)
    if minimap != 0 :
        print("Error while running minimap2: " + command + "\n")
        sys.exit(1)

    #parse the paf file and get the CIGARs
    orientation_repolished = ">"
    #the four followning values should be adjusted, but in case there are problems during alignment, give them a default value
    pos_on_read_left = beginning_toPolish
    pos_on_contig_left = new_link[0]
    with open("tmp_left_and_right.paf", "r") as f:
        for line in f:
            if line[0] != "@" :
                ls = line.strip().split("\t")
                if ls[0] == "left" :
                    if ls[4] == "-":
                        orientation_repolished = "<"

                    if new_link[1] == "+":
                        pos_on_contig_left = int(ls[3]) + start_left
                    else :
                        pos_on_contig_left = int(ls[1]) - int(ls[3]) + start_left

                    if ls[4] == "+" :
                        pos_on_read_left = int(ls[8])
                    else :
                        pos_on_read_left = int(ls[7])

    # print("The new link will go from pos ", pos_on_contig_left, " (",contigs[contigName].length,") of ", contigName, " to pos ", pos_on_contig_right, " of ", new_link[2], " (",contigs[new_link[2]].length,") and between that stands the polished sequence coord ", pos_on_read_left, " to ", pos_on_read_right)

    #to avoid extremely small hanging contigs
    if new_link[1] == "-" and pos_on_contig_left < 20 :
        pos_on_contig_left = 0 
    elif new_link[1] == "+" and pos_on_contig_left > contigs[contigName].length - 20 :
        pos_on_contig_left = contigs[contigName].length

    seq_repolished = ""
    with open("tmp_repolished.fa", "r") as f:
        for line in f:
            if line[0] != ">" :
                seq_repolished += line.strip()
    seq_repolished = seq[min(pos_on_read_left, len(seq_repolished)):]
    if len(seq_repolished) < 100: #if the sequence is too short, don't create a new contig
        return

    nameOfNewContig = "new_" + contigName +"_"+ str(pos_on_contig_left) + "_*"
    new_contig = Contig(nameOfNewContig, seq_repolished)
    if orientation_repolished == ">" :
        new_contig.neighbors.add((0, "-", contigName, pos_on_contig_left, new_link[1], "0M"))
        contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], nameOfNewContig, 0, "-", "0M"))

    else :
        new_contig.neighbors.add((len(seq_repolished), "-", contigName, pos_on_contig_left, new_link[1], "0M"))
        contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], nameOfNewContig, len(seq_repolished), "-", "0M"))

    contigs[nameOfNewContig] = new_contig


    # print("eeicjizlxlz ")
    # sys.exit(1)

#function to output the gfa file from the contigs
def output_gfa(contigs, nameOfFile):

    #list the breakpoins
    new_contigs = {}
    name_of_new_contigs = {} #associating a (contigName, position) to the names of the two new contigs
    for contig in contigs.keys():
        breakpoints = []
        for neighbor in contigs[contig].neighbors:
            breakpoints.append(neighbor[0])

        breakpoints = list(set(breakpoints))
        breakpoints.sort()

        # if contig == "edge_8":
        #     print("bp of edge_8 : ", breakpoints)

        if len(breakpoints) == 0 or breakpoints[0] != 0 :
            breakpoints = [0] + breakpoints
        if breakpoints[-1] != contigs[contig].length :
            breakpoints = breakpoints + [contigs[contig].length]

        for i in range(len(breakpoints)-1):

            new_contigs[contig + "_" + str(breakpoints[i]) + "_" + str(breakpoints[i+1])] = Contig(contig + "_" + str(breakpoints[i]) + "_" + str(breakpoints[i+1]), contigs[contig].sequence[breakpoints[i]:breakpoints[i+1]])
            
            name_of_new_contigs[(contig, breakpoints[i])] = [contig + "_" + str(breakpoints[i]) + "_" + str(breakpoints[i+1])]
            if i > 0 :
                name_of_new_contigs[(contig, breakpoints[i])].append(contig + "_" + str(breakpoints[i-1]) + "_" + str(breakpoints[i]))
                #link the two pieces of the contig
                contigLeftName = contig + "_" + str(breakpoints[i-1]) + "_" + str(breakpoints[i])
                contigRightName = contig + "_" + str(breakpoints[i]) + "_" + str(breakpoints[i+1])
                new_contigs[contigLeftName].neighbors.add((new_contigs[contigLeftName].length, "+", contigRightName, 0, "-", "0M"))
                # print("xxixui adding link : ", contigLeftName, " ", (new_contigs[contigLeftName].length, "+", contigRightName, 0, "-", "0M"))

        # print("iddidizf ", contig, " ", breakpoints)
        name_of_new_contigs[(contig, breakpoints[-1])] = ["*", contig + "_" + str(breakpoints[-2]) + "_" + str(breakpoints[-1])]

    #create the links
    for contig in contigs.keys():
        for neighbor in contigs[contig].neighbors:
            #test to make sure not to write duplicate links
            if contig <= neighbor[2] or (contig == neighbor[2] and neighbor[0] < neighbor[3]) :
                #let's find the name of the contigs
                nameOfThisContig, nameOfOtherContig = "", ""
                positionOnThisContig, positionOnOtherContig = 0, 0
                if neighbor[1] == "+":
                    nameOfThisContig = name_of_new_contigs[(contig, neighbor[0])][1]
                    positionOnThisContig = new_contigs[nameOfThisContig].length
                else :
                    nameOfThisContig = name_of_new_contigs[(contig, neighbor[0])][0]
                    positionOnThisContig = 0
                
                if neighbor[4] == "+" :
                    nameOfOtherContig = name_of_new_contigs[(neighbor[2], neighbor[3])][1]
                    positionOnOtherContig = new_contigs[nameOfOtherContig].length
                else :
                    nameOfOtherContig = name_of_new_contigs[(neighbor[2], neighbor[3])][0]
                    positionOnOtherContig = 0

                new_contigs[nameOfThisContig].neighbors.add((positionOnThisContig, neighbor[1], nameOfOtherContig, positionOnOtherContig, neighbor[4], neighbor[5]))
                new_contigs[nameOfOtherContig].neighbors.add((positionOnOtherContig, neighbor[4], nameOfThisContig, positionOnThisContig, neighbor[1], neighbor[5]))
            
    #now write the gfa
    with open(nameOfFile, "w") as f:
    
        #first write all the segments
        for contig in new_contigs.keys():
            f.write("S\t" + contig + "\t" + new_contigs[contig].sequence + "\n")

        #now output all the links
        for contig in new_contigs.keys():
            for neighbor in new_contigs[contig].neighbors:
                if neighbor[0] == new_contigs[contig].length and neighbor[3] == 0 :
                    f.write("L\t" + contig + "\t+\t" + neighbor[2] + "\t+\t" + neighbor[5] + "\n")
                elif neighbor[0] == new_contigs[contig].length and neighbor[3] == new_contigs[neighbor[2]].length :
                    f.write("L\t" + contig + "\t+\t" + neighbor[2] + "\t-\t" + neighbor[5] + "\n")
                elif neighbor[0] == 0 and neighbor[3] == 0 :
                    f.write("L\t" + contig + "\t-\t" + neighbor[2] + "\t+\t" + neighbor[5] + "\n")
                elif neighbor[0] == 0 and neighbor[3] == new_contigs[neighbor[2]].length :
                    f.write("L\t" + contig + "\t-\t" + neighbor[2] + "\t-\t" + neighbor[5] + "\n")
                else :
                    # print("What is this link ? ", neighbor)
                    # print("contig ", contig, " has length ", contigs[contig].length, " and neighbor ", neighbor[2], " has length ", contigs[neighbor[2]].length)
                    sys.exit(1)

class Mapping:
    def __init__(self):
        self.length_of_read = 0
        self.read = ""

        self.contig1 = ""
        self.position_on_contig1 = 0
        self.orientation_on_contig1 = 0
        self.pos_on_read1 = 0
        self.orientation_on_read1 = 0
        self.breakpoint1 = False

        self.contig2 = ""
        self.position_on_contig2 = 0
        self.orientation_on_contig2 = 0
        self.pos_on_read2 = 0
        self.orientation_on_read2 = 0
        self.breakpoint2 = False

    def __str__(self):
        return self.read + " " + self.contig1 + " " + self.contig2 + " " + str(self.breakpoint1) + " " + str(self.breakpoint2) \
            + " " + str(self.pos_on_read1) + " " + str(self.pos_on_read2) + " " + str(self.position_on_contig1) + " " \
            + str(self.position_on_contig2) + " " + str(self.orientation_on_contig1) + " " + str(self.orientation_on_contig2) \
            + " " + str(self.orientation_on_read1) + " " + str(self.orientation_on_read2)

def main():

    #parse the agruments
    args = parse_args()
    assembly = args.assembly
    reads = args.reads
    output = args.output
    threads = args.threads
    minigraph = args.minigraph
    tmp_folder = args.folder.strip("/") + "/"

    #check the dependencies
    minigraph_res = os.system(minigraph + " --version >& /dev/null")
    if minigraph_res != 0:
        print("Error: could not find/run minigraph at specified location: ", minigraph)
        exit(1)

    racon_res = os.system(args.racon + " --version >& /dev/null")
    if racon_res != 0:
        print("Error: could not find/run racon at specified location: ", args.racon)
        exit(1)

    minimap2_res = os.system(args.minimap2 + " --version >& /dev/null")
    if minimap2_res != 0:
        print("Error: could not find/run minimap2 at specified location: ", args.minimap2)
        exit(1)

    #parse the gfa and fill the contigs list
    contigs = {} #this is the name of the ocntigs, associated to the contigs in the assembly graph
    with open(assembly, "r") as f:
        for line in f:
            if line.startswith("S"):
                ls = line.strip().split("\t")
                contigs[ls[1]] = Contig(ls[1], ls[2])

            if line.startswith("L"):
                ls = line.strip().split("\t")
                pos_this_contig = contigs[ls[1]].length
                orientation_on_this_contig = "+"
                if ls[2] == "-":
                    pos_this_contig = 0
                    orientation_on_this_contig = "-"

                orientation_on_other_contig = "-"
                pos_other_contig = 0
                if ls[4] == "-":
                    pos_other_contig = contigs[ls[3]].length
                    orientation_on_other_contig = "+"
                contigs[ls[1]].neighbors.add((pos_this_contig, orientation_on_this_contig, ls[3], pos_other_contig, orientation_on_other_contig, ls[5]))
                contigs[ls[3]].neighbors.add((pos_other_contig, orientation_on_other_contig, ls[1], pos_this_contig, orientation_on_this_contig, ls[5]))

    #now map the reads on the assembly graph using minigraph
    alignmentFile = tmp_folder + "tmp.gaf"
    minigraph_res = os.system(minigraph + " -c -N 0 -t " + str(threads) + " " + assembly + " " + reads + " > " + alignmentFile + " 2>" + tmp_folder + "log_minigraph.txt")
    if minigraph_res != 0:
        print("Error: could not run minigraph, was trying to run ", minigraph + " -c -t " + threads + " " + assembly + " " + reads + " > " + alignmentFile)
        exit(1)

    #parse the alignment file : associate to each read the contig it aligns on (index only reads that do not align end-to-end)
    #reads_breakpoints = {}
    contig_mapping = {}
    with open(alignmentFile, "r") as f:
        for line in f:
            ls = line.strip().split("\t")
            if ls[4] == "-" :
                print("WHOA a .gaf with strand -, I'm panicking now... SOS, error code SOS !")
                sys.exit(1)
            
            min_length_for_breakpoint = min(0.2*int(ls[1]), 100)
            if (int(ls[2]) > min_length_for_breakpoint or int(ls[1])-int(ls[3]) > min_length_for_breakpoint) and ls[11] == "60" : #read is not mapped end to end
                # if ls[0] == "12009_@3372ceef-1e80-79d4-2f1a-87a9f49b30dc" :
                #     print("read ", ls[0], " is not mapped end to end on contig ", line)
                if ls[0] not in contig_mapping :
                    contig_mapping[ls[0]] = []
                all_contigs = re.split('[><]' , ls[5])
                del all_contigs[0]
                orientations = "".join(re.findall("[<>]", ls[5]))

                m = Mapping()
                m.read = ls[0]
                m.length_of_read = int(ls[1])

                m.contig1 = all_contigs[0]
                m.position_on_contig1 = int(ls[7])
                m.pos_on_read1 = int(ls[2])
                m.orientation_on_contig1 = "-"
                m.orientation_on_read1 = "-"
                if orientations[0] == "<":
                    m.position_on_contig1 = contigs[all_contigs[0]].length - int(ls[7])
                    m.orientation_on_contig1 = "+"

                m.contig2 = all_contigs[-1]
                m.position_on_contig2 = int(ls[8]) - int(np.sum([contigs[i].length for i in all_contigs[:-1]]))
                m.pos_on_read2 = int(ls[3])
                m.orientation_on_contig2 = "+"
                m.orientation_on_read2 = "+"
                if orientations[-1] == "<":
                    m.position_on_contig2 = contigs[all_contigs[-1]].length - (int(ls[8]) - int(np.sum([contigs[i].length for i in all_contigs[:-1]])))
                    m.orientation_on_contig2 = "-"
                
                if int(ls[2]) > min_length_for_breakpoint : #then breakpoint towards the beginning of the read
                    m.breakpoint1 = True

                if int(ls[1])-int(ls[3]) > min_length_for_breakpoint: #then breakpoint towards the end of the read
                    m.breakpoint2 = True
                
                contig_mapping[ls[0]].append(m)

                # #let's look if we have a breakpoint towards the beginning of the read
                # if int(ls[2]) > min_length_for_breakpoint : #then breakpoint towards the beginning of the read
                #     contig = all_contigs[0]
                #     orientation_on_contig = "-" #indicates the direction of the unknown
                #     position_on_contig = int(ls[7])
                #     if orientations[0] == "<":
                #         position_on_contig = contigs[contig].length - int(ls[7])
                #         orientation_on_contig = "+"
                    
                #     orientation_on_read = "-" #indicates the direction of the unknown
                #     pos_on_read = int(ls[2])
                #     length_of_read = int(ls[1])

                #     reads_breakpoints[ls[0]] += [(contig, position_on_contig, orientation_on_contig, ls[0], pos_on_read, orientation_on_read, length_of_read)]

                # #let's look if we have a breakpoint towards the end of the read
                # if int(ls[1])-int(ls[3]) > min_length_for_breakpoint: #then breakpoint towards the end of the read
                #     contig = all_contigs[-1]
                #     orientation_on_contig = "+" #indicates the direction of the unknown
                #     position_on_contig = int(ls[8]) - int(np.sum([contigs[i].length for i in all_contigs[:-1]]))
                #     if orientations[-1] == "<":
                #         position_on_contig = contigs[contig].length - (int(ls[8]) - int(np.sum([contigs[i].length for i in all_contigs[:-1]])))
                #         orientation_on_contig = "-"
                    
                #     orientation_on_read = "+"
                #     pos_on_read = int(ls[3])
                #     length_of_read = int(ls[1])

                #     reads_breakpoints[ls[0]] += [(contig, position_on_contig, orientation_on_contig, ls[0], pos_on_read, orientation_on_read, length_of_read)]

    for read in contig_mapping.keys():
        contig_mapping[read].sort(key=lambda x: int(x.pos_on_read1))

        #see if there are links between the breakpoints
        b = 0
        for b in range (len(contig_mapping[read])):

            if contig_mapping[read][b].breakpoint1 == True : #breakpoint into the unknown
                contig = contig_mapping[read][b].contig1
                position_on_contig = contig_mapping[read][b].position_on_contig1
                orientation_on_contig = contig_mapping[read][b].orientation_on_contig1
                position_on_read = contig_mapping[read][b].pos_on_read1
                orientation_on_read = contig_mapping[read][b].orientation_on_read1
                length_of_read = contig_mapping[read][b].length_of_read

                if position_on_read > max(300, 0.2*length_of_read) :
                    #in this case store in the "other contig" (which does not exist), the read
                    if orientation_on_read == "+" :
                        contigs[contig].potential_new_links.add((position_on_contig, orientation_on_contig, "*", "*", "*", \
                                                                 read, "-", length_of_read-position_on_read, length_of_read))
                        
            elif b < len(contig_mapping[read])-1 and contig_mapping[read][b].breakpoint2 == True and contig_mapping[read][b+1].breakpoint1 \
                and abs(contig_mapping[read][b+1].pos_on_read1-contig_mapping[read][b].pos_on_read2) < 300: #then the two breakpoints are linked !!

                contig_mapping[read][b+1].breakpoint1 = False

                contig = contig_mapping[read][b].contig2
                position_on_contig = contig_mapping[read][b].position_on_contig2
                orientation_on_contig = contig_mapping[read][b].orientation_on_contig2
                position_on_read = contig_mapping[read][b].pos_on_read2
                other_contig = contig_mapping[read][b+1].contig1
                position_on_other_contig = contig_mapping[read][b+1].position_on_contig1
                orientation_on_other_contig = contig_mapping[read][b+1].orientation_on_contig1
                other_position_on_read = contig_mapping[read][b+1].pos_on_read1

                contigs[contig].potential_new_links.add((position_on_contig, orientation_on_contig, other_contig, position_on_other_contig, orientation_on_other_contig,\
                                                        read, "+", position_on_read, other_position_on_read))
                contigs[other_contig].potential_new_links.add((position_on_other_contig, orientation_on_other_contig, contig, position_on_contig, orientation_on_contig,\
                                                        read, "-", other_position_on_read, position_on_read))
            
            elif contig_mapping[read][b].breakpoint2 == True : #then breakpoint into the unknown
                contig = contig_mapping[read][b].contig2
                position_on_contig = contig_mapping[read][b].position_on_contig2
                orientation_on_contig = contig_mapping[read][b].orientation_on_contig2
                position_on_read = contig_mapping[read][b].pos_on_read2
                orientation_on_read = contig_mapping[read][b].orientation_on_read2
                length_of_read = contig_mapping[read][b].length_of_read

                if length_of_read - position_on_read > max(300,0.2*length_of_read) :
                    contigs[contig].potential_new_links.add((position_on_contig, orientation_on_contig, "*", "*", "*", \
                                                                read, "+", position_on_read, length_of_read))


        
    # for read in reads_breakpoints:
    #     #sort the alignments by position on the read, from left to right
    #     reads_breakpoints[read].sort(key=lambda x: int(x[4]))
        
        # if read == "12009_@a8ec3e4e-9869-21c9-590a-7690551170bd" :
        #for each breakpoint, see if it is linked with the next breakpoint (missing link in GFA) or if it is a leap in the unknown (missing contig in GFA)
        # b = 0
        # while b < len(reads_breakpoints[read]):

        #     contig = reads_breakpoints[read][b][0]
        #     position_on_contig = reads_breakpoints[read][b][1]
        #     orientation_on_contig = reads_breakpoints[read][b][2]

        #     #check if the next breakpoint is right next on the read
        #     # print(b, " ", reads_breakpoints[read][b], " ", len(reads_breakpoints[read]))
        #     if b < len(reads_breakpoints[read]) -1 and reads_breakpoints[read][b][5] == "+" and reads_breakpoints[read][b+1][5] == "-" and reads_breakpoints[read][b+1][4]-reads_breakpoints[read][b][4] < 300: #then the two breakpoints are linked are linked !!
        #         other_contig = reads_breakpoints[read][b+1][0]
        #         position_on_other_contig = reads_breakpoints[read][b+1][1]
        #         orientation_on_other_contig = reads_breakpoints[read][b+1][2]

        #         #add the potential new link to the contigs
        #         contigs[contig].potential_new_links.add((position_on_contig, orientation_on_contig, other_contig, position_on_other_contig, orientation_on_other_contig,\
        #                                                  read, "+", reads_breakpoints[read][b][4], reads_breakpoints[read][b+1][4]))
        #         contigs[other_contig].potential_new_links.add((position_on_other_contig, orientation_on_other_contig, contig, position_on_contig, orientation_on_contig,\
        #                                                  read, "-", reads_breakpoints[read][b+1][6]-reads_breakpoints[read][b+1][4], reads_breakpoints[read][b+1][6]-reads_breakpoints[read][b][4]))

        #         # if contig == "edge_31" or other_contig == "edge_31" :
        #         #     print("Adding a link because of read ", read, " ", contig, " ", position_on_contig, " ", orientation_on_contig, " ", other_contig, " ", position_on_other_contig, " ", orientation_on_other_contig)

        #         b += 1 #because the next breakpoint is already linked to this one
            
        #     else : #means a breakpoint going in an unknown contig 
        #         position_on_read = reads_breakpoints[read][b][4]
        #         orientation_on_read = reads_breakpoints[read][b][5]
        #         length_of_read = reads_breakpoints[read][b][6]
        #         #check if there is a significant overhang of the read
        #         if (orientation_on_read == "-" and position_on_read > 300) or (orientation_on_read == "+" and length_of_read - position_on_read > 300) :
        #             #in this case store in the "other contig" (which does not exist), the read
        #             # print("qfoidpdfsp")
        #             if orientation_on_read == "+" :
        #                 contigs[contig].potential_new_links.add((position_on_contig, orientation_on_contig, "*", "*", "*", \
        #                                                          read, "+", position_on_read, length_of_read))
        #             else :
        #                 contigs[contig].potential_new_links.add((position_on_contig, orientation_on_contig, "*", "*", "*", \
        #                                                          read, "-", length_of_read-position_on_read, length_of_read))
        
            # b += 1

    #build a fastq_index, i.e. a dict associating each read with its position in the fastq file
    fastq_index = {}
    offset = 0
    with open(reads, "r") as f:
        for i, line in enumerate(f):
            offset += len(line)
            if i % 4 == 0:
                fastq_index[line.strip().split()[0][1:]] = offset

    #go through the contigs and create all the new links
    all_contig_names = list(contigs.keys()) #doing this as not to iterate on contigs, which size will change
    for contig in all_contig_names:
        if contig == "edge_8" or True:
            contigs[contig].find_new_links()
            #now create the new links
            for link in contigs[contig].future_links.keys():
                # print("piypio ", link)
                if link[2] != "*" and (contig < link[2] or (contig == link[2] and link[0] < link[3])) : #create an order on the links so that we do not create the same link twice
                    generate_new_link(contig, link, contigs[contig].future_links[link], contigs, reads, fastq_index)
                # elif link[2] == "*" :
                #     create_new_contig(contig, link, contigs[contig].future_links[link], contigs, reads, fastq_index)

    #now output the gfa
    # print("Now outputting the graph")
    output_gfa(contigs, output)

if __name__ == "__main__":
    main()
