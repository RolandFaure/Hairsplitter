#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

import numpy as np
import sys
#import matplotlib.pyplot as plt

import input_output as io
from bisect import bisect_left #to look through sorted lists

from copy import deepcopy
import os

from transform_gfa import check_segments
import segment as s
from segment import add_link
from segment import Segment

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

#this function detects and breaks up long (>length) chimeric contigs
def break_up_chimeras(segments, names, interactionMatrix, length) :
    
    allsegments = []
    allXs = []
    for s, segment in enumerate(segments) :
        
        if segment.length > length :
            
            interactions = []
            X = []
            
            for axis in range(1, len(segment.names)) :
                
                interaction = 0
                
                for nameLeft in segment.names[axis:] :
                    for nameRight in segment.names[:axis] :

                        
                        interaction += interactionMatrix[names[nameLeft], names[nameRight]]
                        
                interactions += [interaction]
                
                if axis > 1 :
                    X += [ X[-1] + segment.lengths[axis-1]/segment.length ]
                else :
                    X = [segment.lengths[0]/segment.length]
                
            allsegments += [interactions]
            
            #plt.plot (X, interactions)
            
            inSlump = False
            localMinimums = []
            for axis in range(1, len(interactions)-1) :
                
                if interactions[axis] < 0.7*np.max(interactions[:axis]) and interactions[axis] < 0.7*np.max(interactions[axis:]):
                    if not inSlump :
                        inSlump = True
                        
                    localMinimums += [axis]
                        
                
                else :
                    if inSlump :
                        inSlump = False
                        
                        loin = [interactions[i] for i in localMinimums].index( np.min([interactions[i] for i in localMinimums]) )
                        
                        print('Breaking up contig ', segment.names, ' between ', segment.names[localMinimums[loin]-1], ' and ', segment.names[localMinimums[loin]], ' because it looks like a chimeric contig')
                                
                        #Now break the contig where it should
                        newSegment1, newSegment2 = segment.break_contig(localMinimums[loin])
                        segments[s] = newSegment1
                        segments.append(newSegment2)
                        
                        localMinimums = []

        
    #plt.show()    
    return segments


# input : one supercontig to be joined with a neighbor at one end
# output : actualized listOfSegments with the two contigs merged
def merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments):

    if len(segment.links[endOfSegment]) != 1:
        print("ERROR : trying to merge simply two contigs that cannot be merged simply")
        return -1, -1
    
    neighbor = segment.links[endOfSegment][0]
    endOfSegmentNeighbor = segment.otherEndOfLinks[endOfSegment][0]

    if len(neighbor.links[endOfSegmentNeighbor]) != 1 :
        print('ERROR : trying to merge simply two contigs that cannot be merged simply')
        return -1,-1
        
    if neighbor == segment :  # then do not merge a contig with itself
        return -1, -1

    # add the new segment
    s.merge_two_segments(segment, endOfSegment, neighbor, listOfSegments)
    

    # delete links going towards the two ex-segments
    otherEnd = 1 - endOfSegment
    otherEndNeighbor = 1 - endOfSegmentNeighbor

    for i, n in enumerate(segment.links[otherEnd]) :
        n.remove_end_of_link(segment.otherEndOfLinks[otherEnd][i], segment, otherEnd)
        
    for i, n in enumerate(neighbor.links[otherEndNeighbor]) :
        #print('Removing ', neighbor.names, ' from ', n.names, ' and adding the new contig',listOfSegments[-1].names, ' at end ', neighbor.otherEndOfLinks[otherEndNeighbor][i])
        n.remove_end_of_link(neighbor.otherEndOfLinks[otherEndNeighbor][i], neighbor, otherEndNeighbor)

    # delete the ex-segments
    listOfSegments.remove(segment)
    listOfSegments.remove(neighbor)
    
    return listOfSegments

# input : one supercontig to be joined with a list of neighbors at one end (if they can be merged simply, i.e. if they have only one neighbor at the other end)
# output : actualized listOfSegments with the contigs merged
def merge_simply_all_contigs_on_one_side_of_this_one(segment, endOfSegment, listOfSegments): 

    if len(segment.links[endOfSegment]) == 0:
        return listOfSegments

    list_of_segments_to_merge = [segment]
    list_of_names = segment.names[::2*endOfSegment-1]
    if endOfSegment == 1:
        list_of_orientations = segment.orientations
    else:
        list_of_orientations = [1-segment.orientations[i] for i in range(len(segment.orientations)-1,-1,-1)]
    list_of_lengths = segment.lengths[::2*endOfSegment-1]
    list_of_insideCIGARs = segment.insideCIGARs[::2*endOfSegment-1]
    list_of_read_coverages = segment.depths[::2*endOfSegment-1]
    end_of_contig_now = endOfSegment
    contig_now = segment
    already_merged_contigs = set([segment.ID])
    while len(contig_now.links[end_of_contig_now]) == 1 and len(contig_now.links[end_of_contig_now][0].links[contig_now.otherEndOfLinks[end_of_contig_now][0]]) == 1 and contig_now.links[end_of_contig_now][0].ID not in already_merged_contigs:
        
        list_of_insideCIGARs += contig_now.CIGARs[end_of_contig_now]

        new_end_of_contig = 1-contig_now.otherEndOfLinks[end_of_contig_now][0]
        contig_now = contig_now.links[end_of_contig_now][0]
        end_of_contig_now = new_end_of_contig
        already_merged_contigs.add(contig_now.ID)

        list_of_segments_to_merge.append(contig_now)
        list_of_names += contig_now.names[::2*end_of_contig_now-1]
        list_of_lengths += contig_now.lengths[::2*end_of_contig_now-1]
        list_of_read_coverages += contig_now.depths[::2*end_of_contig_now-1]
        list_of_insideCIGARs += contig_now.insideCIGARs[::2*end_of_contig_now-1]

        if end_of_contig_now == 1:
            list_of_orientations += contig_now.orientations
        else:
            list_of_orientations += [1-contig_now.orientations[i] for i in range(len(contig_now.orientations)-1,-1,-1)]

    #create a new segment that is the merge of all the segments in list_of_segments_to_merge
    newSegment = Segment(list_of_names, \
                        list_of_orientations, \
                        list_of_lengths, \
                        list_of_insideCIGARs, \
                        segLinks = [list_of_segments_to_merge[0].links[1-endOfSegment], list_of_segments_to_merge[-1].links[end_of_contig_now]], \
                        segOtherEndOfLinks = [list_of_segments_to_merge[0].otherEndOfLinks[1-endOfSegment], list_of_segments_to_merge[-1].otherEndOfLinks[end_of_contig_now]], \
                        segCIGARs = [list_of_segments_to_merge[0].CIGARs[1-endOfSegment], list_of_segments_to_merge[-1].CIGARs[end_of_contig_now]], \
                        lock = True, \
                        HiCcoverage = sum([i.HiCcoverage for i in list_of_segments_to_merge]), \
                        readCoverage=list_of_read_coverages)    

    #create the links from the neighbors of the new segment
    for i, neighbor in enumerate(newSegment.links[0]):
        neighbor.add_end_of_link(newSegment.otherEndOfLinks[0][i], newSegment, 0, newSegment.CIGARs[0][i])
    for i, neighbor in enumerate(newSegment.links[1]):
        neighbor.add_end_of_link(newSegment.otherEndOfLinks[1][i], newSegment, 1, newSegment.CIGARs[1][i])

    
    #add the new segment to the list of segments
    listOfSegments.append(newSegment)


    #delete all the segments that have been merged
    for s in list_of_segments_to_merge:
        s.cut_all_links()
        listOfSegments.remove(s)

    return listOfSegments

# a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):

    print("merge_adjacent_contigs")
    goOn = True
    while goOn:
        goOn = False
        for segment in listOfSegments:

            alreadyDidThisOne = False # if the segment is deleted when looking at its first end, you don't want it to look at its other end, since it does not exist anymore
            for endOfSegment in range(2):
                
                if not alreadyDidThisOne:
                    if len(segment.links[endOfSegment]) == 1 and len(segment.links[endOfSegment][0].links[segment.otherEndOfLinks[endOfSegment][0]])== 1 :  # then merge
                        alreadyDidThisOne = True
                        if segment.ID != segment.links[endOfSegment][0].ID:
                            goOn = True
                            listOfSegments = merge_simply_all_contigs_on_one_side_of_this_one(segment, endOfSegment, listOfSegments)



    return listOfSegments

# input : a list of segments
# output : a list of segments, with the segments that can be duplicated duplicated
def duplicate_contigs(segments):

    print("Duplicating contigs")

    continueDuplication = True
    while continueDuplication:
        continueDuplication = False
        toDelete = []
        alreadyDuplicated = set()
        originalLength = len(segments)
        for se in range (originalLength) :
            contig = segments[se]
            #check if the segment can be duplicated by one of its end
            for end in range(2):
                if contig.ID not in alreadyDuplicated:

                    if len(contig.links[end]) > 1 and len(contig.links[1-end]) > 1 \
                        and all([len(contig.links[end][i].links[contig.otherEndOfLinks[end][i]]) == 1 for i in range(len(contig.links[end]))]) \
                        and all([len(contig.links[1-end][i].links[contig.otherEndOfLinks[1-end][i]]) == 1 for i in range(len(contig.links[1-end]))]) \
                        and (contig.depth > 0.7*np.sum([contig.links[end][i].depth for i in range(len(contig.links[end]))]) or contig.length < 1000) \
                        and all([neighbor.depth > 0.2*contig.depth for neighbor in contig.links[end]]) \
                        and contig.ID not in [i.ID for i in contig.links[end]] :

                        #then duplicate the segment !
                        alreadyDuplicated.add(contig.ID)
                        numberofcopies = len(contig.links[end])
                        toDelete += [se]
                        totalNeighborCoverage = np.sum([contig.links[end][i].depth for i in range(len(contig.links[end]))])
                        if totalNeighborCoverage == 0:
                            totalNeighborCoverage = 1
                        for n in range(numberofcopies):
                            percentageOfDepth = contig.links[end][n].depth/totalNeighborCoverage
                            newSegment = Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i*percentageOfDepth for i in contig.depths])
                            segments.append(newSegment)
                            add_link(segments[-1], end, contig.links[end][n] , contig.otherEndOfLinks[end][n], "0M")
                            for on, otherneighbor in enumerate(contig.links[1-end]) :
                                add_link(newSegment, 1-end, otherneighbor, contig.otherEndOfLinks[1-end][on], "0M")

        for s in sorted(toDelete, reverse=True):
            continueDuplication = True
            segments[s].cut_all_links()
            del segments[s]

        # segments = merge_adjacent_contigs(segments)

    return segments

#input: segments
#output : segments trimmed to reduce the number of overlaps
def trim_overlaps(segments):

    something_changes = True
    while something_changes :
        something_changes = False

        for s in segments :
            min_overlap_left = 0
            if len(s.CIGARs[0]) != 0 :
                #check that there are nothing else than Ms in the CIGAR before trimming
                only_M = True
                for cig in s.CIGARs[0] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    min_overlap_left = min([int(i.strip('M')) for i in s.CIGARs[0]])
            min_overlap_right = 0
            if len(s.CIGARs[1]) != 0 :
                only_M = True
                for cig in s.CIGARs[1] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    min_overlap_right = min([int(i.strip('M')) for i in s.CIGARs[1]])

            max_overlap_left = 0
            if len(s.CIGARs[0]) != 0 :
                only_M = True
                for cig in s.CIGARs[0] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    max_overlap_left = max([int(i.strip('M')) for i in s.CIGARs[0]]+[0])
            max_overlap_right = 0
            if len(s.CIGARs[1]) != 0 :
                only_M = True
                for cig in s.CIGARs[1] :
                    for l in cig :
                        if (l > '9' or l < '0') and l != 'M':
                            only_M = False
                if only_M :
                    max_overlap_right = max([int(i.strip('M')) for i in s.CIGARs[1]]+[0])

            already_trimmed = s.get_trim()

            trim_left = min(min_overlap_left, s.length-max_overlap_right)
            trim_right = min(min_overlap_right, s.length-max_overlap_left)


            s.set_trim([already_trimmed[0]+trim_left, already_trimmed[1]+trim_right])

            if (trim_left > 0 or trim_right > 0) :
                something_changes = True

            leftCIGARs = [str(int(i.strip('M'))-trim_left)+'M' for i in s.CIGARs[0]]
            for l in range(len(s.links[0])) :
                newcigar = leftCIGARs[l]
                s.links[0][l].set_CIGAR(s.otherEndOfLinks[0][l], s, 0, newcigar)
                s.set_CIGAR(0, s.links[0][l], s.otherEndOfLinks[0][l], newcigar)

            rightCIGARs = [str(int(i.strip('M'))-trim_right)+'M'for i in s.CIGARs[1]]
            for l in range(len(s.links[1])) :
                newcigar = rightCIGARs[l]
                s.links[1][l].set_CIGAR(s.otherEndOfLinks[1][l], s, 1, newcigar)
                s.set_CIGAR(1, s.links[1][l], s.otherEndOfLinks[1][l], newcigar)

            # if (trim_left > 0 or trim_right > 0) and s.names == ['edge_165@0@0', 'edge_167@0@0', 'edge_169@0@0'] :
            #     print("qfdjqsdjlm after changeds : ", s.CIGARs[0], " ", s.CIGARs[1])

            

    return segments

#function that merges all adjacent contigs but that start from a GFA and not a list of segments and output a GFA
#does not use the complex contig/link data structure defined above, as it takes a lot of time
def merge_adjacent_contigs_GFA(gfa_in, gfa_out):

    #go through the GFA and index all contigs, and all links
    segments = [] #just a list of segment names
    depths = {} #associates a segment with its depth
    lengths = {} #associates a segment with its length
    links = {} #associates a segment with two lists of tuples ([(neighbor, end, CIGAR), ...], [(neighbor, end, CIGAR), ...])
    location_of_segments_in_GFA_file = {} #associates a segment with its location in the GFA file (to recover the sequence quickly)

    with open(gfa_in, 'r') as f:

        line = f.readline()
        while line :

            if line[0] == 'S' :
                length_line = len(line)
                line = line.strip().split('\t')
                location_of_segments_in_GFA_file[line[1]] = f.tell()-length_line
                segments.append(line[1])
                if line[1] not in links :
                    links[line[1]] = [[], []]
                lengths[line[1]] = len(line[2])
                depths[line[1]] = 0
                #look for the depth
                for field in line[3:] :
                    if field[:2] == "DP" :
                        depths[line[1]] = float(field[5:])
                        break

            elif line[0] == 'L' :
                line = line.split('\t')
                if line[1] not in links :
                    links[line[1]] = [[], []]
                if line[3] not in links :
                    links[line[3]] = [[], []]

                side1 = 1
                side2 = 0
                if line[2] == '-' :
                    side1 = 0
                if line[4] == '-' :
                    side2 = 1

                links[line[1]][side1].append((line[3], side2, line[5].strip()))
                if line[3] != line[1] or side1 != side2 : #to add only one link if it is a self link
                    links[line[3]][side2].append((line[1], side1, line[5].strip()))

            line = f.readline()

    #merge the contigs
    old_segments_to_new_segments = {} #associates the old segment name with (the new segment name, endOfTheNewSegment) #if it is not at an end put -1 we dont care there will be no link
    new_segments = [] #each new segment is a list [(segment, orientation, CIGAR), ...]
    for segment in segments :

        for end in range(2):
            if segment not in old_segments_to_new_segments :
                # print("Merging ", segment, " at end ", end)
                #check if this is the end of a new contig
                if len(links[segment][end]) != 1 or len(links[links[segment][end][0][0]][links[segment][end][0][1]]) != 1 :

                    #then start a new contig and go on until end of new contig
                    new_segment = [(segment, 1-end, "0M")]
                    endNow = 1-end #end from which we are trying to extend
                    segmentNow = segment
                    already_merged_segments = set([segment])
                    while len(links[segmentNow][endNow]) == 1 and len(links[links[segmentNow][endNow][0][0]][links[segmentNow][endNow][0][1]]) == 1 and links[segmentNow][endNow][0][0] not in already_merged_segments:

                        new_segment.append((links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1], links[segmentNow][endNow][0][2]))

                        segmentNow, endNow = links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1]
                        already_merged_segments.add(segmentNow)

                    #now the new segment is created, fill the old_segments_to_new_segments
                    new_segment_name = "_".join([i[0] for i in new_segment])
                    for i, s in enumerate(new_segment):
                        if i== 0 and i == len(new_segment)-1 :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, 2)
                        elif i == 0 :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, 0)
                        elif i == len(new_segment)-1 :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, 1)
                        else :
                            old_segments_to_new_segments[s[0]] = (new_segment_name, -1)

                    new_segments.append(new_segment)

    #some contigs were in loops with only contigs that have one neighbor left and right and still have not been merged, so do this now
    for segment in segments :

        if segment not in old_segments_to_new_segments : #then it is in a loop
            
            end = 0
            #then start a new contig and go on until end of new contig
            new_segment = [(segment, 1-end, "0M")]
            endNow = 1-end #end from which we are trying to extend
            segmentNow = segment
            already_merged_segments = set([segment])
            while len(links[segmentNow][endNow]) == 1 and len(links[links[segmentNow][endNow][0][0]][links[segmentNow][endNow][0][1]]) == 1 and links[segmentNow][endNow][0][0] not in already_merged_segments:

                new_segment.append((links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1], links[segmentNow][endNow][0][2]))

                segmentNow, endNow = links[segmentNow][endNow][0][0], 1-links[segmentNow][endNow][0][1]
                already_merged_segments.add(segmentNow)

            #now the new segment is created, fill the old_segments_to_new_segments
            new_segment_name = "_".join([i[0] for i in new_segment])
            for i, s in enumerate(new_segment):
                if i== 0 and i == len(new_segment)-1 :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, 2)
                elif i == 0 :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, 0)
                elif i == len(new_segment)-1 :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, 1)
                else :
                    old_segments_to_new_segments[s[0]] = (new_segment_name, -1)

            new_segments.append(new_segment)
        

    #write the new GFA
    # print("outputting gfa")
    gfa = open(gfa_in, 'r')
    with open(gfa_out, 'w') as f:

        L_lines = []
        number_out = 0

        for new_segment in new_segments :

            # print("Writing new segment ", number_out)
            number_out += 1

            #compute the sequence of the new segment
            seq = ""
            depth_total = 0
            length_total = 0
            # print("Looking for sequence of ", new_segment)
            for seg in new_segment :

                # print("Looking for sequence of ", seg[0])
                length_total += lengths[seg[0]]
                depth_total += depths[seg[0]]*lengths[seg[0]]

                line_subsegment = ""
                gfa.seek(location_of_segments_in_GFA_file[seg[0]])
                line_subsegment = gfa.readline()
                line_subsegment = line_subsegment.strip().split('\t')
                if len(line_subsegment) < 3 :
                    seq_subsegment = ""
                else :
                    seq_subsegment = line_subsegment[2]
                gfa.seek(0)

                length_overlap = int(seg[2].strip('M'))

                if seg[1] == 0 :
                    seq_subsegment = reverse_complement(seq_subsegment)

                seq += seq_subsegment[length_overlap:]

                # print("After having added ", seg[0], " ", seq)

            #write the new segment
            new_segment_name = "_".join([i[0] for i in new_segment])
            f.write("S\t"+new_segment_name+"\t"+seq+"\tDP:f:"+ str(depth_total/(length_total+1)) +"\n")

            #now check if the new segment has links at its ends
            #first left end
            for l in links[new_segment[0][0]][1-new_segment[0][1]] :
                old_neighbor = l[0]
                new_neighbor_name = old_segments_to_new_segments[old_neighbor][0]
                new_neighbor_end = old_segments_to_new_segments[old_neighbor][1]
                if new_neighbor_end == 2: #means this is a non-merged segment, fall back to the original orientation
                    new_neighbor_end = l[1]

                if new_neighbor_end == -1 :
                    print("ERROR : new segment has a link at its left end with a segment that is not at an end")
                    sys.exit(1)

                if new_segment_name < new_neighbor_name : #to avoid writing the same link twice
                    L_line = "L\t"+new_segment_name+"\t-\t" + new_neighbor_name + "\t"+ str("+-"[new_neighbor_end]) +"\t" + l[2] + "\n"
                    L_lines.append(L_line)

            #then right end
            for l in links[new_segment[-1][0]][new_segment[-1][1]] :
                old_neighbor = l[0]
                new_neighbor_name = old_segments_to_new_segments[old_neighbor][0]
                new_neighbor_end = old_segments_to_new_segments[old_neighbor][1]
                if new_neighbor_end == 2: #means this is a non-merged segment, fall back to the original orientation
                    new_neighbor_end = l[1]

                if new_neighbor_end == -1 :
                    print("ERROR : new segment has a link at its right end with a segment that is not at an end")
                    sys.exit(1)

                if new_segment_name <= new_neighbor_name : #to avoid writing the same link twice
                    L_lines.append("L\t"+new_segment_name+"\t+\t" + new_neighbor_name + "\t"+ "+-"[new_neighbor_end] +"\t" + l[2] + "\n")

        for l in L_lines :
            f.write(l)

    gfa.close()
    f.close()

                







              

    

