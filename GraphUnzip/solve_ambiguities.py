#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:30:45 2020

File dedicated to the algorithm af making bigger contigs, including solving bubbles
"""

import numpy as np
#import matplotlib.pyplot as plt

import input_output as io
from bisect import bisect_left #to look through sorted lists

from copy import deepcopy
import os

from transform_gfa import check_segments
import segment as s
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
                        
                        # nameLeft = nameLeft.split('-')[0]
                        # nameRight = nameRight.split('-')[0]
                        
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

# this function measures the intensity of interactions between one supercontig
# and several candidate, including without taking account of the common parts
# of the supercontigs

def intensity_of_interactions(
    segment,
    endOfSegment,
    listOfNeighborIndices,
    listOfSegments,
    interactionMatrix,
    names, 
    copiesnumber,
    debugDir = '',
):
    
    candidatesSegments = [segment.links[endOfSegment][i] for i in listOfNeighborIndices]
    listOfTouchingEnds = [segment.otherEndOfLinks[endOfSegment][i] for i in listOfNeighborIndices]
    
    for candidate in candidatesSegments :
        if candidate == segment : #small loop, don't solve that !
            return [-1]
    
    ##first compute all contigs common to all candidates, to take them out

    #commonContigs, neighborsOfNeighborsUsed = compute_commonContigs(segment, candidatesSegments, listOfTouchingEnds, depthOfCommonContigs)
    commonContigs = compute_commonContigs2(segment, endOfSegment, listOfNeighborIndices, 2000000)
    
    if debugDir != '':
        f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
        f.write('common contigs: ' + str(commonContigs)+'\n')
        f.close()
    

    ##Now compute the score of each candidates    

    # bestsignature decides what contig is most characterisic of the segment
    bestSignature = np.min([copiesnumber[x] for x in segment.names])
    # return for each supercontig its absolute score and its relative score (wihtout the common parts)
    
    relativeScores = []
    returnRelativeScore = True

    for c in candidatesSegments:
                
        absoluteScore, relativeScore, depthHere = c.interaction_with_contigs(segment, interactionMatrix, names, copiesnumber, commonContigs, bestSignature)
            
        if all([i in commonContigs for i in c.names]) :
            returnRelativeScore = False #if a contig is entirely contained in commonContigs, don't say that his score is 0, that would lead to errors
    
        relativeScores.append(relativeScore)

    # if 'edge_229' in segment.names:    
        # print('At contig ', segment.names, ' choosing between ',  [i.names for i in candidatesSegments], ' and the result is ', relativeScores, absoluteScores)
        # #print('Best signature : ', bestSignature, ' and the signatures are : ', [copiesnumber[x] for x in segment.names])
        # print('Common contigs : ', commonContigs, '\n')
    
    if returnRelativeScore :
        return relativeScores
    else :
        return [-1]

#input : a segment, a end of segment a list of neighboring contig, a depth (length in nucleotides) beyond which we do not need to look
#output : a list of contig's names that are present in the prolongation of only one of the neighbor
def compute_commonContigs2(segment, endOfSegment, listOfNeighborIndices, depth) :
    
    listOfTouchedContig = []
    
    for n in listOfNeighborIndices :
        vicinityContigs = set()
        propagate_vicinity(segment.links[endOfSegment][n], segment.otherEndOfLinks[endOfSegment][n], vicinityContigs, 0, depth)
        listOfTouchedContig += [set()]
        for comm in vicinityContigs :
            listOfTouchedContig[-1].add(comm.split('$:')[0])
        
    commonContigs = set.intersection(*listOfTouchedContig)
    
    for name in segment.names :
        commonContigs.add(name) #to ensure no interaction is computed with itself
    
    return commonContigs

#a recursive function to find all contigs within the limitDepth of the end of a contig
def propagate_vicinity(segment, endOfSegment, vicinityContigs, depth, limitDepth) :
    
    if depth > limitDepth :
        return 0
    else :
        for name in segment.names :
            vicinityContigs.add(name+"$:"+str(endOfSegment)+"$:"+str(segment.ID))
        for n, neighbor in enumerate(segment.links[1-endOfSegment]):
            if neighbor.names[0]+"$:"+str(segment.otherEndOfLinks[1-endOfSegment][n])+"$:"+str(neighbor.ID) not in vicinityContigs :
                propagate_vicinity(neighbor, segment.otherEndOfLinks[1-endOfSegment][n], vicinityContigs, depth+segment.length , limitDepth)
        return 0

def compute_commonContigs(segment, candidatesSegments, listOfTouchingEnds, depth) :

    commonContigs = []
    if depth == 0 :
        return commonContigs
        
    if depth == 2 :
    #first compute the list of common contigs counting neighbors and neighbors of neighbors of segment
        potentialCommonContigs = candidatesSegments[0].names.copy()
    
        for i in candidatesSegments[0].links[1-listOfTouchingEnds[0]] :
            potentialCommonContigs += i.names
    
    
        for contig in potentialCommonContigs:
            
    
            for n in range(1, len(candidatesSegments)):
    
                candidate = candidatesSegments[n]
                presentInNeighbor = contig in candidate.names
    
                for i in candidate.links[1-listOfTouchingEnds[n]] :                      
                    presentInNeighbor = presentInNeighbor or (contig in i.names)
                    
    
            if presentInNeighbor:
                commonContigs += [contig]
                
        #check if that list of common contigs is not too big
        neighborOfNeighborActivated = True
        for c in candidatesSegments :
            if all([elem in commonContigs for elem in c.names]) : #this is not good, commoncontigs is too big
                neighborOfNeighborActivated = False
        
        if neighborOfNeighborActivated :
            return commonContigs, True #the True value is to signifie that neighbors of neighbors were used
    
    
    
    #recompute common contigs but only with neighbors
    commonContigs = []
    potentialCommonContigs = candidatesSegments[0].names.copy()
        
    for contig in potentialCommonContigs:
    
        common = True

        for n, candidate in enumerate(candidatesSegments[1:]):

            common = common and (contig in candidate.names)

        if common:
            commonContigs += [contig]
        
    return commonContigs, False #the False value is to signifie that neighbors of neighbors were not used


# small function to look in a sorted list l if x is present in logarithmic time
def isPresent(l, x):
    i = bisect_left(l, x)
    if i != len(l) and l[i] == x:
        return True
    return False


# here we look specifically at one contig and its immediate surroundings (can
# return -1 if fails in short loop)
def duplicate_around_this_end_of_contig(
    segment, endOfSegment, listOfSuperContigs, copiesnumber
):  # endOfSegment should be 0 if it's the left end and1 if it's the right end

    if (
        segment in segment.links[endOfSegment]
    ):  # if a segment loops on itself, another module would be needed
        return 0
    if any(
        segment.links[endOfSegment].count(i) >= 2 for i in segment.links[endOfSegment]
    ):  # if segment.links[endOfSegment] has two copies of the same segment, it means one link going towards each end, that is not solvable
        return 0

    for i in segment.names:
        copiesnumber[i] += len(segment.links[endOfSegment]) - 1

    # add all the new supercontigs
    segment.divide_depths(len(segment.links[endOfSegment]))
    for neighbor in segment.links[endOfSegment]:
        s.merge_two_segments(
            segment, endOfSegment, neighbor, listOfSuperContigs
        )  # the merged segment is appended at the end of listOfSuperContigs

    # now delete the merged supercontigs
    # start by deleting the links that linked the merged supercontigs to the outside

    deletedContigs = [segment.ID]
    otherEnd = 1 - endOfSegment

    for i, neighbor in enumerate(segment.links[otherEnd]):

        neighbor.remove_end_of_link(
            segment.otherEndOfLinks[otherEnd][i], segment, otherEnd
        )

    for m, merged in enumerate(segment.links[endOfSegment]):

        if (
            len(merged.links[segment.otherEndOfLinks[endOfSegment][m]]) == 1
        ):  # then the original copy is fully integrated in the supercontig
            deletedContigs.append(merged.ID)

            otherEnd = 1 - segment.otherEndOfLinks[endOfSegment][m]
            for i, neighbor in enumerate(merged.links[otherEnd]):
                try:
                    neighbor.remove_end_of_link(
                        merged.otherEndOfLinks[otherEnd][i], merged, otherEnd
                    )
                except ValueError:  # that means we're in a small loop which we can't solve
                    print(
                        "There is merging difficulty around the far end of "
                        + str(merged.names)
                        + " from "
                        + str(segment.names)
                        + " . Please check that there is indeed a loop there."
                    )
                    return 0

        else:  # then the original contig still exists by itself, just delete the link going toward segment
            try:
                merged.remove_end_of_link(
                    segment.otherEndOfLinks[endOfSegment][m], segment, endOfSegment
                )
            except ValueError:  # that means we're in a small loop which whe can't solve
                print("There is merging difficulty around the near end of "+ str(merged.names)+ " from "+ str(segment.names)+ " . Please check that there is indeed a loop there.")
                return 0

    # delete all segments that should be
    deletedContigs.sort()
    for i in range(len(listOfSuperContigs) - 1, -1, -1):
        h = listOfSuperContigs[i].ID
        if isPresent(
            deletedContigs, h
        ):  # in other words, if h is in deletedContigs (written like that because it has logarithmic efficiency)
            del listOfSuperContigs[i]

    # lock all the segments that have been duplicated, so that they are not duplicated by both ends
    segment.lockNode(endOfSegment)

# similar to the function above, but simpler: put in one supercontig two smaller supercontig linked by a link unambinguous at both ends
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
                    
                    if len(segment.links[endOfSegment]) == 1and len(segment.links[endOfSegment][0].links[segment.otherEndOfLinks[endOfSegment][0]])== 1:  # then merge
                        alreadyDidThisOne = True
                        if segment != segment.links[endOfSegment][0]:
                            goOn = True
                            listOfSegments = merge_simply_two_adjacent_contig(segment, endOfSegment, listOfSegments)

    return listOfSegments

#merge_contigs looks at choices endofsegment by endofsegment, and duplicates all the necessary contigs
def merge_contigs(listOfSegments, copiesnumber, multiplicities, names, verbose = False):
    # each contig can be multiplied only once in this function (to ensure it is not multiplied by both ends) : once it is duplicated, it gets locked

    #look at both ends of each segment sequentially
    for segment in listOfSegments:
        
            
        for endOfSegment in range(2) :      

            l = segment.links[endOfSegment]
 
            if len(l) > 1 and not segment.freezed[endOfSegment]:
                
                startMerging = not segment.locked
                for i in l:
                    startMerging = startMerging and (not i.locked)
                    
                   
                if startMerging:  # if nothing is locked for now
                    
                    #check that the read coverage is adequate with duplication of the contig
                    #if  segment.depths[-endOfSegment] >= max([segment.links[endOfSegment][n].depths[-segment.otherEndOfLinks[endOfSegment][n]] for n in range(len(segment.links[endOfSegment]))]) :
                    if multiplicities[names[segment.names[-endOfSegment]]] >= len(segment.links[endOfSegment])* copiesnumber[segment.names[-endOfSegment]]:
                        duplicate_around_this_end_of_contig(segment,endOfSegment,listOfSegments,copiesnumber)
                    
                    if verbose :
                        print('Duplicating contig ', segment.names, ' around its end touching ', [i.names for i in listOfSegments])
                    
    #now that the duplicating is done, unlock all the segments
    for segment in listOfSegments :
        segment.locked = False
        
    # now just merge all two contigs that are next to each other
    listOfSegments = merge_adjacent_contigs(listOfSegments)
    
    return listOfSegments, copiesnumber
                 
#input : a graph. This function takes out all links that are not confirmed by long reads
def check_all_links(segments, lr_links) :
    
    
    for segment in segments :
        
        for endOfSegment in range(2) :
            
            toRemove = set()
            
            for n, neighbor in enumerate(segment.links[endOfSegment]) :
                
                
                endOfNeighbor = segment.otherEndOfLinks[endOfSegment][n]
                
                #define the link in the format of the lr_links
                link = (segment.names[-endOfSegment], (segment.orientations[-endOfSegment] == endOfSegment), neighbor.names[-endOfNeighbor], (neighbor.orientations[-endOfNeighbor] == endOfNeighbor))
                linkb = (link[2], link[3], link[0], link[1])
                
                                
                if not link in lr_links and not linkb in lr_links : #then the link is not confirmed by long reads
                    
                    toRemove.add((segment, endOfSegment, neighbor, endOfNeighbor))
            
            for i in toRemove :
                
                i[0].remove_end_of_link(i[1], i[2], i[3])
                i[2].remove_end_of_link(i[3], i[0], i[1])
    
#get_rid_of_bad_links compare links using HiC contact informations when there is a choice and delete links that are not supported by long-range evidence
def get_rid_of_bad_links(listOfSegments, interactionMatrix, tagInteractionMatrix, names, copiesnumber,thresholdRejected,thresholdAccepted, debugDir = '', neighborsOfNeighbors = True, verbose = False):

    HiCmatrix = (interactionMatrix.count_nonzero() > 0) #a boolean value to tell if there is need to use the Hi-C interaction matrix

    #loop through all segments inspecting the robustness of all links.
    c = 0
    
    #then compute the intensity of interactions knowing the common contigs

    for segment in listOfSegments:
        
        c += 1

        for endOfSegment in range(2):
                
            if len(segment.links[endOfSegment]) >= 2 : #then it means that there is a choice to be made at one end of the segment. Let's see how HiC contacts confirm those links
                    
                    
                # comparison pairwise of the links, those that should be deleted are deleted
                toRemove = []
                for n1 in range(len(segment.links[endOfSegment]) - 1):
                    n2 = n1 + 1
                    while n2 < len(segment.links[endOfSegment]):
                                                     
                        #first compute using linked reads    
                        if tagInteractionMatrix.count_nonzero()>0 and (linksStrength == [-1] or (all([i>1 for i in linksStrength]) or all([i<=1 for i in linksStrength]))) :
                            linksStrength = intensity_of_interactions(segment, endOfSegment, [n1,n2],listOfSegments, tagInteractionMatrix, names, copiesnumber, debugDir = debugDir)
                                    
                        #if it is not enough, use Hi-C
                        linksStrength = [-1]
                        if (linksStrength == [-1] or (all([i>1 for i in linksStrength]) or all([i<=1 for i in linksStrength]))) and HiCmatrix:
                            linksStrength = intensity_of_interactions(segment, endOfSegment, [n1,n2], listOfSegments, interactionMatrix, names, copiesnumber, debugDir = debugDir)
                        
                        if debugDir != '':
                            f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
                            f.write('I have to decide, at '+'_'.join(segment.names)+ ' between '+ '_'.join(segment.links[endOfSegment][n1].names)+ ' and '+'_'.join(segment.links[endOfSegment][n2].names) + ' with these values : '+ str(linksStrength)+'\n')
                            f.close()
                            
                        if verbose : #or 'edge_289' in segment.names :
                            print('I have to decide, at '+'_'.join(segment.names)+ ' between '+ '_'.join(segment.links[endOfSegment][n1].names)+ ' and '+'_'.join(segment.links[endOfSegment][n2].names) + ' with these values : '+ str(linksStrength))

                        if linksStrength == [-1]: #means that the configuration does not enable the algorithm to compare the two interactions
                            segment.freezeNode(endOfSegment) 

                        elif any([i>1 for i in linksStrength]): #the condition is to prevent too much duplicating if there is no mapping or almost  
                            
                            # if '262' in segment.names :
                            #     print('I have to decide, at '+'_'.join(segment.names)+ ' between '+ '_'.join(segment.links[endOfSegment][n1].names)+ ' and '+'_'.join(segment.links[endOfSegment][n2].names) + ' with these values : '+ str(linksStrength)+'\n')
                            #     print([i.names for i in segment.links[0]])
                            if linksStrength[0] > linksStrength[1]:
                                if (linksStrength[1] <= linksStrength[0] * thresholdRejected) or (linksStrength[1] == 1 and linksStrength[0] > 2):  # then it means that probably link does not exist
                                
                                    #however, we'd like not to create any dead end, hence a little precaution:
                                    if len(segment.links[endOfSegment][n2].links[segment.otherEndOfLinks[endOfSegment][n2]]) == 1:
                                        segment.freezeNode(endOfSegment) #n2 is weak, but this edge is the only outgoing edge from the contig at the other end
                                    else:
                                        if verbose :
                                            print('\nRemoving link from ', segment.links[endOfSegment][n2].names, ' to ', segment.names, '\n')
                                        if n2 not in toRemove :
                                            toRemove += [n2]
                                    
                                elif (linksStrength[1] < linksStrength[0] * thresholdAccepted):  # then it's not clear, the link is freezed
                                    segment.freezeNode(endOfSegment)

                            else:
                                if linksStrength[0] < linksStrength[1] * thresholdRejected or (linksStrength[0] == 1 and linksStrength[1] > 2):  # then decide that the link does not exist
                                    #however, we'd like not to create any dead end, hence a little precaution:
                                    if len(segment.links[endOfSegment][n1].links[segment.otherEndOfLinks[endOfSegment][n1]]) == 1:
                                        segment.freezeNode(endOfSegment) #n2 is weak, but this edge is the only outgoing edge from the contig at the other end
                                    else:
                                        if verbose :
                                            print('\nRemoving link from ', segment.links[endOfSegment][n1].names, ' to ', segment.names, '\n')
                                        if n1 not in toRemove :
                                            toRemove += [n1]
                                    
                                    
                                elif linksStrength[0] < linksStrength[1] * thresholdAccepted:  # then it's not clear, the link is freezed
                                    segment.freezeNode(endOfSegment)
                        else : #linksStrength <= [1,1]
                            segment.freezeNode(endOfSegment)
                            # print('get_rid_of_bad_links, ...  freeznoding2 : ' + '\t'.join( ['_'.join(segment.links[endOfSegment][n1].names), '_'.join(segment.links[endOfSegment][n2].names)])+'\n')

                        n2+=1
                        
                #Remove all links that have been marked as removable
                toRemove.sort()
                toRemove.reverse()
                
                for n in toRemove :
                    segment._links[endOfSegment][n].remove_end_of_link(segment._otherEndOfLinks[endOfSegment][n], segment, endOfSegment)
                    segment.remove_end_of_link(endOfSegment, segment._links[endOfSegment][n], segment._otherEndOfLinks[endOfSegment][n])
                    
        
    return listOfSegments                        

def stats_on_thresholds(segments, names, interactionMatrix, copiesNumber) :

    ratios = []
    for segment in segments :
        
        for endOfSegment in range(2) :
            
            if len(segment.links[endOfSegment]) >= 2 : #then it means that there is a choice to be made at one end of the segment. Let's see how HiC contacts confirm those links
                    
                # comparison pairwise of the links, those that should be deleted are deleted
                    for n1 in range(len(segment.links[endOfSegment]) - 1):
                        n2 = n1 + 1
                        while n2 < len(segment.links[endOfSegment]):
                            
                            d = 2
                               
                            absoluteLinksStrength, linksStrength, neighborsOfNeighborsUsed = intensity_of_interactions(segment, [segment.links[endOfSegment][n1], segment.links[endOfSegment][n2]],\
                                                                                             [segment.otherEndOfLinks[endOfSegment][n1], segment.otherEndOfLinks[endOfSegment][n2]],\
                                                                                             segments, interactionMatrix, names, copiesNumber, depthOfCommonContigs = d)
                                                                                             
                            n2 += 1
                            
                            if len(linksStrength) == 2 :
                                ratios += [ np.min(linksStrength)/np.max(linksStrength) ]
                                
    # plt.hist(ratios)
    # plt.xlabel('i(X)/i(Y) ratio')
    # plt.ylabel('Number of ambiguities having this value')
    # plt.show()
    
    return ratios

def solve_ambiguities(listOfSegments, interactionMatrix, tagInteractionMatrix, multiplicities, names, stringenceReject, stringenceAccept, steps, copiesNumber = {}, useNeighborOfNeighbor = True, debugDir = '', verbose = False):
    

    if debugDir != '' :
        if not os.path.isdir(debugDir) :
            os.mkdir(debugDir)
        f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'w')
        f.close()
    
    if copiesNumber == {} :
        for segment in listOfSegments :
            for name in segment.names :
                cn[name] = 1
    
    listOfSegments = merge_adjacent_contigs(listOfSegments)
    print('Merged adjacent contigs for the first time')

   # s.check_if_all_links_are_sorted(listOfSegments)

    for i in range(steps):
        get_rid_of_bad_links(listOfSegments, interactionMatrix, tagInteractionMatrix, names, copiesNumber, stringenceReject, stringenceAccept, debugDir = debugDir, neighborsOfNeighbors = useNeighborOfNeighbor, verbose = verbose)
                
        #solve_small_loops(listOfSegments, names, repeats, lr_links, check_links)
        
        #solve_l_loops(listOfSegments, lr_links)
            
        print('Got rid of bad links')

        # for se in listOfSegments :
        #     if 'edge_357' in se.names :
        #         print ('Here is two : ', se.names, [i.names for i in se.links[0]], [i.names for i in se.links[1]], '\n') 
        
        listOfSegments, copiesNumber = merge_contigs(listOfSegments, copiesNumber, multiplicities, names, verbose = verbose)
        
        #stats_on_thresholds(listOfSegments, names, interactionMatrix, copiesNumber)
        #stats_on_thresholds(listOfSegments, names, lrInteractionMatrix, copiesNumber)
        # while True:
        #     r = 0
        
        #print('end of merge_contigs : ', [i.names for i in listOfSegments[names['262']].links[0]])

        # once all the contigs have been duplicated and merged, unfreeze everything so the cycle can start again
        for j in listOfSegments:
            j.unfreeze()
        
        print(str((i+1) / steps * 100) + "% of solving ambiguities done")
        
        if debugDir != '' :
            io.export_to_GFA(listOfSegments, exportFile = debugDir.strip('/')+'/'+'debug_gfa_step'+str(i)+'.gfa')
            f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
            f.write('Finished step '+ str(i)+ ' \n\n\n')
            f.close()
            
    #finish by breaking up long chimeras that can form sometimes
    
    if interactionMatrix.count_nonzero() > 0 :
        listOfSegments =  break_up_chimeras(listOfSegments, names, interactionMatrix, 100000)
        
    return listOfSegments, copiesNumber #return copiesNumber in case you want to run solve_ambiguities several times in a row



