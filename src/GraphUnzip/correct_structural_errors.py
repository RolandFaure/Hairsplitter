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

from input_output import load_gfa
from input_output import export_to_GFA
from finish_untangling import merge_adjacent_contigs
from finish_untangling import trim_overlaps
from simple_unzip import simple_unzip
from repolish import repolish_contigs
import segment as sg


resolution = 50

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return "".join(complement[base] for base in reversed(seq))

#function that detach and destroy short tips
#input : a list of segments
#output : the same list of segments, but with the short tips detached
def detach_tips(segments):
                                              
    #detach a short dead ends
    changes = True
    max_tip_length = 1000
    to_delete = set()
    while changes :
        changes = False
        for s, seg in enumerate(segments):
            
            for end in range(2) :

                if len(seg.links[end]) > 1 : #one of the branch may be a short dead end    
                    extended_lengths = [extended_length(seg.links[end][i], seg.otherEndOfLinks[end][i], max_tip_length*5, 20) for i in range(len(seg.links[end]))]
                    max_length = max(extended_lengths)
                    toDelete = set()
                    for n in range(len(seg.links[end])) :
                        if 5*extended_lengths[n] < max_length and max_length > 10000 and extended_lengths[n] < max_tip_length :
                            toDelete.add((seg, end, seg.links[end][n], seg.otherEndOfLinks[end][n]))
                            changes = True 
                    for seg, end, neighbor, otherEnd in toDelete :
                        sg.delete_link(seg, end, neighbor, otherEnd, warning=False)
                        #also delete the neighbor and all its neighbors since it is a dead end
                        if len(neighbor.links[1-otherEnd]) == 0 and len(neighbor.links[otherEnd]) == 0 :
                            to_delete.add(neighbor)
                            # if neighbor.names == ['edge_339_0_8617']:
                            #     print("Deleting ", neighbor.names, " because it is a dead end from ", seg.names)
                            #     sys.exit()

        #now delete the segments that have been marked for deletion
        for seg in to_delete :
            seg.cut_all_links()
            segments.remove(seg)
        to_delete = set()

    return segments

#a function returning the longest path you can get with neighbors of neighbors of neighbors of... (up to threshold) 
def extended_length(segment, end, thresholdLength, thresholdContigs) :
    
    #print("Extended length called with threshold ", thresholdLength, " on segment , ", sg.names)
    
    if thresholdContigs == 0 or thresholdLength <= 0 :
        return segment.length
    
    #start by looking down the longest contig, it will be fastest
    longestContig = [i for i in range(len(segment.links[1-end]))]
    longestContig.sort(key= lambda x : segment.links[1-end][x].length, reverse = True)
    
    #print("Longest contigs : ", longestContig, [segment.links[1-end][i].length for i in longestContig])
    
    maxLength = 0
    for n in longestContig :
        neighbor = segment.links[1-end][n]
        l = extended_length(neighbor, segment.otherEndOfLinks[1-end][n], thresholdLength-segment.length, thresholdContigs-1)
        if l > maxLength :
            maxLength = l
        if maxLength > thresholdLength :
            break
        
    return maxLength + segment.length

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
                length_on_read1 = abs(already_link[8]-already_link[7])
                length_on_read2 = abs(link[8]-link[7])
                if (link[2] == "*" and already_link[2] == "*") and abs(link[0]-already_link[0]) < 200 and link[1] == already_link[1] :
                    links[already_link].append(link[-4:])
                    added = True
                    break
                elif link[2] != "*" and link[1] == already_link[1] and link[2]==already_link[2] and link[4] == already_link[4] \
                    and abs(link[0]-already_link[0]) < max(length_on_read1, length_on_read2)+200 and abs(link[3]-already_link[3]) < max(length_on_read1, length_on_read2)+200\
                    and length_on_read1 < 1000 and length_on_read2 < 1000:

                    links[already_link].append(link[-4:])
                    added = True
                    if length_on_read2 < length_on_read1: #then take this read as the best one, since it is the one that does the smallest leap
                        links[link] = [link[-4:]]+links[already_link].copy()[:-1] #put this read in the first place, it will be the easiest to polish
                        del links[already_link]
                    break

            if not added :
                links[link] = [link[-4:]]
               
        #now go through the links and see if they are present at least 10 times
        # print("All the links for contig ", self.name, " are ")
        for k in links.keys() :
            # print("linkd : ", k, " ", links[k])
            if len(links[k]) > 4:
                self.future_links[k] = links[k]
                # print("yeeess")
                

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
    for r,read in enumerate(reads):
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
            start, end = len(seq)-start, len(seq)-end

        seq = seq[max(0,start-300):min(end+300, len(seq))] #take a little margin on both sides
        if (foundToPolish == False and len(seq) > 600) or (foundToPolish == False and nread == nreads-1):
            foundToPolish = True
            beginning_toPolish = min(300, start)
            end_toPolish = max(len(seq)-300, len(seq)-end)
            foPolish = open(toPolish_file, "w")
            foPolish.write(">" + read_name+ "_" + str(r) + "\n" + seq + "\n")
            foPolish.close()

        if len(seq) == 0 :
            print("Error, empty sequence for read ", read_name, " ", start, " ", end, " ", read)
            sys.exit(1)
        fo.write(">" + read_name+ "_" + str(r) + "\n" + seq + "\n")
        nread += 1

    fo.close()

    #polish the joint using racon
    command = "minimap2 -x map-ont tmp_toPolish.fa tmp_reads.fa > tmp.paf 2> trash.txt"
    minimap = os.system(command)
    if minimap != 0 :
        print("Error while running minimap2 1: " + command + "\n")
        sys.exit(1)

    command = "racon tmp_reads.fa tmp.paf tmp_toPolish.fa > tmp_repolished.fa 2>trash.txt"
    racon = os.system(command)
    if racon != 0:
        print()
        print("Error while running racon: " + command)
        print("Trying to build the link ", new_link)
        # print(reads)
        # sys.exit(1)

        os.system("cp tmp_toPolish.fa tmp_repolished.fa")
        # pos_on_contig_left = new_link[0]
        # pos_on_contig_right = new_link[3]
        # contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], new_link[2], pos_on_contig_right, new_link[4], "0M"))
        # contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4], contigName, pos_on_contig_left, new_link[1], "0M"))
        # return #give up on this link
        # print("Trying to build the link ", new_link, " between ", contigName, " and ", new_link[2])
        # sys.exit(1)

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
        print("Error while running minimap2 2: " + command + "\n")
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

    if pos_on_read_right > pos_on_read_left : #then retrieve the sequence
        seq_repolished = ""
        with open("tmp_repolished.fa", "r") as f:
            for line in f:
                if line[0] != ">" :
                    seq_repolished += line.strip()
        seq_repolished = seq[pos_on_read_left:pos_on_read_right]
        nameOfNewContig = "join_" + contigName +"_"+ str(pos_on_contig_left) + "_" + new_link[2]+"_"+str(pos_on_contig_right)
        if nameOfNewContig in contigs.keys() : #oops, same link twice
            return
        new_contig = Contig(nameOfNewContig, seq_repolished)

        if orientation_repolished == ">" :
            new_contig.neighbors.add((0, "-", contigName, pos_on_contig_left, new_link[1], "0M"))
            contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], nameOfNewContig, 0, "-", "0M"))

            new_contig.neighbors.add((len(seq_repolished), "+", new_link[2], pos_on_contig_right, new_link[4], "0M"))
            contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4], nameOfNewContig, len(seq_repolished), "+", "0M"))

        else :
            new_contig.neighbors.add((len(seq_repolished), "+", contigName, pos_on_contig_left, new_link[1], "0M"))
            contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], nameOfNewContig, len(seq_repolished), "+", "0M"))

            new_contig.neighbors.add((0, "-", new_link[2], pos_on_contig_right, new_link[4], "0M"))
            contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4],nameOfNewContig, 0, "-", "0M"))

        contigs[nameOfNewContig] = new_contig

        # print("New contig ", new_contig.name, " has length ", new_contig.length, " and neighbors ", new_contig.neighbors)

    else : #means that the two contigs actually overlap, just add a link between the two
        #just create a new link between the two contigs that wont overlap
        length_of_overlap = pos_on_read_left - pos_on_read_right
        if length_of_overlap == 0 or True: #that should be improved in the future, because there is no guarantee they are really M and moreover the overlap can be longer than the contig
            contigs[contigName].neighbors.add((pos_on_contig_left, new_link[1], new_link[2], pos_on_contig_right, new_link[4], str(length_of_overlap) + "M"))
            contigs[new_link[2]].neighbors.add((pos_on_contig_right, new_link[4], contigName, pos_on_contig_left, new_link[1], str(length_of_overlap) + "M"))


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
        print("Error while running minimap2 3: " + command + "\n")
        sys.exit(1)

    command = "racon tmp_reads.fa tmp.paf tmp_toPolish.fa > tmp_repolished.fa 2>trash.txt"
    racon = os.system(command)
    if racon != 0 :
        print("Error while running racon: " + command + "\n")
        #remove the temporary files
        os.system("rm tmp_reads.fa tmp_toPolish.fa tmp.paf tmp_repolished.fa trash.txt")
        return best_read[0] #give up on this contig
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
        print("Error while running minimap2 4: " + command + "\n")
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

                # print("Link between ", neighbor[0], " ", neighbor[1], " and ", neighbor[3], " ", neighbor[4], " ", neighbor)

                #let's find the name of the contigs
                nameOfThisContig, nameOfOtherContig = "", ""
                positionOnThisContig, positionOnOtherContig = 0, 0
                if neighbor[1] == "+":
                    nameOfThisContig = name_of_new_contigs[(contig, neighbor[0])][1]
                    positionOnThisContig = new_contigs[nameOfThisContig].length
                else :
                    nameOfThisContig = name_of_new_contigs[(contig, neighbor[0])][0]
                    positionOnThisContig = 0
                
                # print("Building link between ", contig, " ", neighbor[0], " and ", neighbor[2], " ", neighbor[3], " with lengths ", contigs[contig].length, " ", contigs[neighbor[2]].length)
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
    minigraph_res = os.system(minigraph + " --version > /dev/null")
    if minigraph_res != 0:
        print("Error: could not find/run minigraph at specified location: ", minigraph)
        exit(1)

    racon_res = os.system(args.racon + " --version > /dev/null")
    if racon_res != 0:
        print("Error: could not find/run racon at specified location: ", args.racon)
        exit(1)

    minimap2_res = os.system(args.minimap2 + " --version > /dev/null")
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
    #first map the reads without the CIGAR to see what reads are not mapped end-to-end
    print("First alignment of the reads on the assembly graph...")
    alignmentFile = tmp_folder + "tmp0.gaf"
    # minigraph_res = os.system(minigraph + " -N 0 -t " + str(threads) + " " + assembly + " " + reads + " > " + alignmentFile + " 2>" + tmp_folder + "log_minigraph.txt")
    not_mapped_end_to_end = set()
    with open(alignmentFile, "r") as f:
        for line in f:
            ls = line.strip().split("\t")            
            min_length_for_breakpoint = min(0.2*int(ls[1]), 100)
            if (int(ls[2]) > min_length_for_breakpoint or int(ls[1])-int(ls[3]) > min_length_for_breakpoint) and ls[11] == "60" :
                not_mapped_end_to_end.add(ls[0])

    #create a new file with only the reads that are not mapped end-to-end
    reads_file = tmp_folder + "reads_not_mapped_end_to_end."+reads.split(".")[-1]
    fi = open(reads, "r")
    fo = open(reads_file, "w")
    write = False
    for line in fi:
        if line[0] == ">" or line[0] == "@" :
            read_name = line.strip().split()[0][1:]
            # print("read name ", read_name)
            if read_name in not_mapped_end_to_end :
                fo.write(line)
                fo.write(fi.readline())
                write = True
            else :
                write = False
        elif write == True :
            fo.write(line)
    fo.close()
    fi.close()
    print("Now aligning the reads that are not mapped end-to-end with basepair resolution...")
    alignmentFile = tmp_folder + "tmp.gaf"
    # minigraph_res = os.system(minigraph + " -c --secondary=no -t " + str(threads) + " " + assembly + " " + reads_file + " > " + alignmentFile + " 2>" + tmp_folder + "log_minigraph.txt")
    # if minigraph_res != 0:
    #     print("Error: could not run minigraph, was trying to run ", minigraph + " -c -t " + threads + " " + assembly + " " + reads_file + " > " + alignmentFile)
    #     exit(1)

    print("Moving on to the breakpoint detection and graph adjustments...")

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

    for read in contig_mapping.keys():
        contig_mapping[read].sort(key=lambda x: int(x.pos_on_read1))

        # if read == "0_c31beb2a-ef35-076b-cb90-a1c8ee84bfba":
        #     print("here iq mf mappign ", [str(i) for i in contig_mapping[read]])

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
                        
                        
            if b < len(contig_mapping[read])-1 and contig_mapping[read][b].breakpoint2 == True and contig_mapping[read][b+1].breakpoint1 :
                #and abs(contig_mapping[read][b+1].pos_on_read1-contig_mapping[read][b].pos_on_read2) < 300: #then the two breakpoints are linked !!

                contig_mapping[read][b+1].breakpoint1 = False

                contig = contig_mapping[read][b].contig2
                position_on_contig = contig_mapping[read][b].position_on_contig2
                orientation_on_contig = contig_mapping[read][b].orientation_on_contig2
                position_on_read = contig_mapping[read][b].pos_on_read2
                other_contig = contig_mapping[read][b+1].contig1
                position_on_other_contig = contig_mapping[read][b+1].position_on_contig1
                orientation_on_other_contig = contig_mapping[read][b+1].orientation_on_contig1
                other_position_on_read = contig_mapping[read][b+1].pos_on_read1

                # if read == "0_c31beb2a-ef35-076b-cb90-a1c8ee84bfba":
                #     print("qspoiten tnew longkqd ", contig, " ", position_on_contig, " ", orientation_on_contig, " ", position_on_read, " ", other_contig, " ", position_on_other_contig, " ", orientation_on_other_contig, " ", other_position_on_read)

                if position_on_read < other_position_on_read : #else the read align on two contigs at once, which is weird but maybe not a link
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
        if contig == "edge_181" or True:
            contigs[contig].find_new_links()
            #now create the new links
            for link in contigs[contig].future_links.keys():
                # print("piypio ", link, " ", contigs[contig].future_links[link])
                if link[2] != "*" and (contig < link[2] or (contig == link[2] and link[0] < link[3])) : #create an order on the links so that we do not create the same link twice
                    # print("New link between ", contig, " and ", link)
                    # print(contigs[contig].future_links[link])
                    potential_new_link = generate_new_link(contig, link, contigs[contig].future_links[link], contigs, reads, fastq_index)
                    n = 0
                    while potential_new_link != None  and n < 3:
                        potential_new_link = generate_new_link(contig, potential_new_link, contigs[contig].future_links[link], contigs, reads, fastq_index)
                        n+= 1

                # elif link[2] == "*" :
                #     create_new_contig(contig, link, contigs[contig].future_links[link], contigs, reads, fastq_index)

    #trim the small dead ends resulting from the stitching

    #now output the gfa
    print("Now outputting the graph")
    output_stitched_gfa = tmp_folder + "stitched.gfa"
    output_gfa(contigs, output_stitched_gfa)

    segments, names = load_gfa(output_stitched_gfa)
    segments = detach_tips(segments)
    segments = merge_adjacent_contigs(segments)

    copies = sg.compute_copiesNumber(segments)

    export_to_GFA(segments, copies, gfaFile=output_stitched_gfa, exportFile=output, merge_adjacent_contigs=True, rename_contigs=False)

    print("Done correcting the graph, the output is in ", output)

if __name__ == "__main__":
    main()

