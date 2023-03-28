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

# a loop to merge all adjacent contigs
def merge_adjacent_contigs(listOfSegments):

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
                            listOfSegments = merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments)


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

                    if len(contig.links[end]) > 1 \
                        and all([len(contig.links[end][i].links[contig.otherEndOfLinks[end][i]]) == 1 for i in range(len(contig.links[end]))]) \
                        and (contig.depth > 0.7*np.sum([contig.links[end][i].depth for i in range(len(contig.links[end]))]) or contig.length < 1000) \
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

        segments = merge_adjacent_contigs(segments)

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





              

    

