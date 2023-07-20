#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:31:09 2021

@author: rfaure
"""
import numpy as np

def interactions_with_neighbors(
    segment,
    endOfSegment,
    candidateSegments,
    listOfNeighborEnds,
    listOfSegments,
    interactionMatrix,
    names, 
    copiesnumber,
    verbose = False,
    debugDir = '',
):
    
    
    for candidate in candidateSegments :
        if candidate.ID == segment.ID : #small loop, don't solve that !
            return [-1]
    
    ##first compute all contigs common to all candidates, to take them out

    commonContigs = []
    if len(candidateSegments) > 1 :
        commonContigs = compute_commonContigs(segment, candidateSegments, listOfNeighborEnds, min(1000000, 3*max([i.length for i in candidateSegments])))
    
    if debugDir != '':
        f = open(debugDir.strip('/')+'/'+'debug_log.txt', 'a')
        #f.write('common contigs: ' + str(commonContigs)+'\n')
        f.close()
    

    ##Now compute the score of each candidates    

    # bestsignature decides what contig is most characterisic of the segment
    bestSignature = np.min([copiesnumber[x] for x in segment.names])
    # return for each supercontig its absolute score and its relative score (wihtout the common parts)
    
    relativeScores = []
    returnRelativeScore = True
    
    absoluteScores = []

    for c in candidateSegments:
                
        absoluteScore, relativeScore, depthHere = c.interaction_with_contigs(segment, interactionMatrix, names, copiesnumber, commonContigs, bestSignature)
            
        if all([i in commonContigs for i in c.names]) :
            returnRelativeScore = False #if a contig is entirely contained in commonContigs, don't say that his score is 0, that would lead to errors
            # if verbose :
            #     print("common contigs between ", [i.names for i in candidateSegments], " : ", commonContigs)
    
        relativeScores.append(relativeScore)
        absoluteScores.append(absoluteScore)

    # if 'edge_27' in segment.names:    
    #     print('At contig ', segment.names, ' choosing between ',  [i.names for i in candidateSegments], ' and the result is ', relativeScores, absoluteScores)
    #     #print('Best signature : ', bestSignature, ' and the signatures are : ', [copiesnumber[x] for x in segment.names])
    #     print('Common contigs : ', commonContigs, '\n')
    
    if returnRelativeScore :
        return relativeScores
    else :
        return [-1]

#input : a segment, a end of segment a list of neighboring contig, a depth (length in nucleotides) beyond which we do not need to look
#output : a list of contig's names that are present in the prolongation of only one of the neighbor
def compute_commonContigs(segment, candidateSegments, listOfNeighborEnds, depth) :
    
    listOfTouchedContig = []
    
    for n in range(len(candidateSegments)) :
        vicinityContigs = set()
        propagate_vicinity(candidateSegments[n], listOfNeighborEnds[n], vicinityContigs, 0, depth, 0,30)
        listOfTouchedContig += [set()]
        for comm in vicinityContigs :
            listOfTouchedContig[-1].add(comm.split('$:')[0])
        
    commonContigs = set.intersection(*listOfTouchedContig)
    
    for name in segment.names :
        commonContigs.add(name) #to ensure no interaction is computed with itself

    
    return commonContigs

#a recursive function to find all contigs within the limitDepth of the end of a contig
def propagate_vicinity(segment, endOfSegment, vicinityContigs, depth, limitDepth, recursionDepth, recursionLimit) :
    
    if depth > limitDepth or recursionDepth >= recursionLimit :
        return 0
    else :
        for name in segment.names :
            vicinityContigs.add(name+"$:"+str(endOfSegment)+"$:"+str(segment.ID))
        for n, neighbor in enumerate(segment.links[1-endOfSegment]):
            if neighbor.names[0]+"$:"+str(segment.otherEndOfLinks[1-endOfSegment][n])+"$:"+str(neighbor.ID) not in vicinityContigs :
                propagate_vicinity(neighbor, segment.otherEndOfLinks[1-endOfSegment][n], vicinityContigs, depth+segment.length , limitDepth, 1+recursionDepth, recursionLimit)
        return 0
