#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from transform_gfa import load_gfa
"""
Created on Tue Apr 28 16:22:21 2020

@author: zaltabar

This file is for function that do not have anymore utility in the master code, but that took time to program
"""

#input : a gaf file (outputted by graphaligner)
#output : an interaction matrix with two values : 100 if the contigs are next to each other in some reads, 0 elsewhise
def longReads_interactionsMatrix(gafFile, names, segments, similarity_threshold = 0, whole_mapping = False):
     
    print('Building interaction matrix from the gaf file')
    f = open(gafFile, 'r')
    
    interactionMatrix = sparse.dok_matrix((len(segments), len(segments)))
    repeats = [0]*len(segments)
    
    allLinks = set()
    noNames = set() #list of names that are not found in names
   
    #these two values are in case some contigs of the gafs are just subnames of names
    namesPlus = names.copy()
    contigsPlus = {}
    for i in names :
        contigsPlus[i] = i
    
    for line in f :
        
        ls = line.split('\t')
        if ls[5].count('>') + ls[5].count('<') > 1 :
            
        
            if (not 'id:f' in ls[-2]) or (float(ls[-2].split(':')[-1]) > similarity_threshold) :
                
                if not whole_mapping or (float(ls[2]) == 0 and float(ls[1]) == float(ls[3])) :
                    
                    contigs = re.split('[><]' , ls[5])
                    orientations = "".join(re.findall("[<>]", ls[5]))
                    del contigs[0] #because the first element is always ''

                    for c1 in range(len(contigs)-1) :
                        for c2 in range(c1+1, len(contigs)):
                            
                            if contigs[c1] not in noNames and contigs[c2] not in noNames :
                            
                                #try to look for the name if it is not in names
                                if contigs[c1] not in namesPlus :
                                    found = False
                                    for i in names.keys() :
                                        if contigs[c1] in i :
                                            namesPlus[contigs[c1]] = names[i]
                                            contigsPlus[contigs[c1]] = i
                                            found = True
                                    if not found :
                                        noNames.add(contigs[c1]) #if you do not find it once, the cause is lost
                                            
                                if contigs[c2] not in namesPlus :
                                    found = False
                                    for i in names.keys() :
                                        if contigs[c2] in i :
                                            namesPlus[contigs[c2]] = names[i]
                                            contigsPlus[contigs[c2]] = i
                                            found = True
                                    if not found :
                                        noNames.add(contigs[c2])
                                
                                #if the two contigs exist, go ahead
                                if contigs[c1] in namesPlus and contigs[c2] in namesPlus :

                                    if namesPlus[contigs[c1]] != namesPlus[contigs[c2]] :
                                        interactionMatrix[namesPlus[contigs[c1]], namesPlus[contigs[c2]]] += 1
                                        if c2 == c1 +1 :
                                            if contigsPlus[contigs[c1]] == contigs[c1] and contigsPlus[contigs[c2]] == contigs[c2] :

                                                allLinks.add((contigsPlus[contigs[c1]], orientations[c1] == '>', contigsPlus[contigs[c2]], orientations[c2] == '<'))
                                                
                                            elif contigsPlus[contigs[c1]] == contigs[c1] :
                                                allLinks.add((contigsPlus[contigs[c1]], orientations[c1] == '>', contigsPlus[contigs[c2]], True))
                                                allLinks.add((contigsPlus[contigs[c1]], orientations[c1] == '>', contigsPlus[contigs[c2]], False))
                                            elif  contigsPlus[contigs[c2]] == contigs[c2] :
                                                allLinks.add((contigsPlus[contigs[c1]], True, contigsPlus[contigs[c2]], orientations[c2] == '<'))
                                                allLinks.add((contigsPlus[contigs[c1]], False, contigsPlus[contigs[c2]], orientations[c2] == '<'))
                                            else :
                                                allLinks.add((contigsPlus[contigs[c1]], True, contigsPlus[contigs[c2]], True))
                                                allLinks.add((contigsPlus[contigs[c1]], False, contigsPlus[contigs[c2]], True))
                                                allLinks.add((contigsPlus[contigs[c1]], True, contigsPlus[contigs[c2]], False))
                                                allLinks.add((contigsPlus[contigs[c1]], False, contigsPlus[contigs[c2]], False))
                                    else :
                                        repeats[namesPlus[contigs[c1]]] = max(contigs.count(contigs[c1]), repeats[namesPlus[contigs[c1]]])
                    
    return interactionMatrix, repeats, allLinks

def how_far_away_are_those_contigs(contig1, contig2, links, infContigs):
    connectedToContig1 = [contig1*2, contig1*2+1]
    
    distanceToContig1 = [0,0]
    
    pathToContig1 = [[],[]]
    
    pastLengthOfconnexion1 = 1
    
    while pastLengthOfconnexion1 < len(connectedToContig1) :
        
        pastLengthOfconnexion1 = len(connectedToContig1)
        newConnectedToContig1 = connectedToContig1[:]
        
        #we spread one step from contig1
        for cc, connectedContig in enumerate(connectedToContig1) : #that can be optimized if we need speed, since we go through the first contigs many times
            
            if cc < 2 or cc%2 == 1 : #that's to move only one way in the contig graph
                for newContig in links[connectedContig] :
                            
                    newdist = distanceToContig1[cc]+infContigs[int(connectedContig/2)][1]
                    if newContig in newConnectedToContig1 :
                        
                        oldpos = newConnectedToContig1.index(newContig)
                        
                        if newdist < distanceToContig1[oldpos] :
                            distanceToContig1[oldpos] = newdist
                            pathToContig1[oldpos] = pathToContig1[cc]+[cc]
                            
                            distanceToContig1[oldpos+1-2*(oldpos%2)] = newdist
                            pathToContig1[oldpos+1-2*(oldpos%2)] = pathToContig1[cc]+[cc]
                    else :
                        #we add both end of the contig to the connected contig
                        newConnectedToContig1 += [newContig]
                        distanceToContig1 += [newdist]
                        pathToContig1 += [pathToContig1[cc]+[cc]]
                        
                        newConnectedToContig1 += [newContig + 1 - 2*(newContig%2)]
                        distanceToContig1 += [newdist]
                        pathToContig1 += [pathToContig1[cc]+[connectedContig]]
        connectedToContig1 = newConnectedToContig1[:]
        
    #now we see where contig2 stands with respect to contig1
    if contig2*2 in connectedToContig1 :
        print(connectedToContig1)
        indexContig2 = connectedToContig1.index(contig2*2)
        return distanceToContig1[indexContig2]-infContigs[contig1][1], pathToContig1[indexContig2][1:] #because the length of contig 1 is always taken into account in the path length while it shouldn't
    else :
        print(connectedToContig1)
        return -1 #meaning contig1 and contig2 are not connected in this graph

def detect_fishy_links(links, confirmationOfLinks, coverage):

    # we're going to detect, when there is an ambiguity, if one path looks unlikely
    badlinks = [] * len(links)
    for i in coverage:  # to ensure we don't get absurdly high multiplicative factors
        if i < 0.01:
            i = 0.01

    for endOfContig in range(len(links)):

        if len(links[endOfContig]) > 1:

            weightedConfirmation = [-1] * len(links[endOfContig])
            for i in range(len(links[endOfContig])):

                if (
                    coverage[int(endOfContig / 2)] > 0.01
                    and coverage[int(links[endOfContig][i] / 2)] > 0.01
                ):
                    weightedConfirmation[i] = (
                        confirmationOfLinks[endOfContig][i]
                        / coverage[int(links[endOfContig][i] / 2)]
                    )

            maximum = np.max(weightedConfirmation)
            reliable = -1 not in weightedConfirmation

            if reliable:
                for i in range(len(weightedConfirmation)):
                    if weightedConfirmation[i] < 0.1 * maximum:

                        badlinks += [[endOfContig, i]]
                        print(
                            "There is a suspect link here : ",
                            endOfContig,
                            links[endOfContig],
                            confirmationOfLinks[endOfContig],
                            weightedConfirmation,
                        )
    return badlinks

def HiC_vs_GFA(hiccontacts, links, fragment_list):

    confirmationOfLinks = [[0 for i in j] for j in links]  # a list of list of 0 of the same dimensions as links

    for contact in hiccontacts:
        contig1 = fragment_list[contact[0]][0]
        contig2 = fragment_list[contact[1]][0]

        for j in range(len(links[contig1 * 2])):
            if (
                links[contig1 * 2][j] == contig2 * 2
                or links[contig1 * 2][j] == contig2 * 2 + 1
            ):
                confirmationOfLinks[contig1 * 2][j] += contact[2]
                for i in range(len(links[links[contig1 * 2][j]])):
                    if links[links[contig1 * 2][j]][i] == contig1 * 2:
                        confirmationOfLinks[links[contig1 * 2][j]][i] += contact[2]

        for j in range(len(links[contig1 * 2 + 1])):
            if (
                links[contig1 * 2 + 1][j] == contig2 * 2
                or links[contig1 * 2 + 1][j] == contig2 * 2 + 1
            ):
                confirmationOfLinks[contig1 * 2 + 1][j] += contact[2]
                for i in range(len(links[links[contig1 * 2 + 1][j]])):
                    if links[links[contig1 * 2 + 1][j]][i] == contig1 * 2 + 1:
                        confirmationOfLinks[links[contig1 * 2 + 1][j]][i] += contact[2]

    return confirmationOfLinks

# same as above but taking into account contigs that are two connexions away
def HiC_vs_GFAtwo(hiccontacts, links, fragment_list, coverage):  
    
    confirmationOfLinks = [
        [0 for i in j] for j in links
    ]  # a list of list of 0 of the same dimensions as links
    weightedconfirmationOfLinks = [[0 for i in j] for j in links]

    for contact in hiccontacts:

        contig1 = fragment_list[contact[0]][0]
        contig2 = fragment_list[contact[1]][0]

        for j, neighbor in enumerate(links[contig1 * 2]):
            # direct neighbor
            if (
                links[contig1 * 2][j] == contig2 * 2
                or links[contig1 * 2][j] == contig2 * 2 + 1
            ):
                confirmationOfLinks[contig1 * 2][j] += contact[2]
                weightedconfirmationOfLinks[contig1 * 2][j] += (
                    contact[2] / coverage[contig1] / coverage[contig2]
                )

                for i in range(len(links[links[contig1 * 2][j]])):
                    if links[links[contig1 * 2][j]][i] == contig1 * 2:
                        confirmationOfLinks[links[contig1 * 2][j]][i] += contact[2]
                        weightedconfirmationOfLinks[links[contig1 * 2][j]][i] += (
                            contact[2] / coverage[contig1] / coverage[contig2]
                        )
            # two connexions away
            for c in range(
                len(links[neighbor + 1 - 2 * neighbor % 2])
            ):  # we take the other end of the neighbor contig
                if (
                    links[neighbor + 1 - 2 * neighbor % 2][c] == contig2 * 2
                    or links[neighbor + 1 - 2 * neighbor % 2][c] == contig2 * 2 + 1
                ):
                    confirmationOfLinks[contig1 * 2][j] += contact[2]
                    weightedconfirmationOfLinks[contig1 * 2][j] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    confirmationOfLinks[neighbor + 1 - 2 * neighbor % 2][c] += contact[
                        2
                    ]
                    weightedconfirmationOfLinks[neighbor + 1 - 2 * neighbor % 2][c] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    for i in range(
                        len(links[links[neighbor + 1 - 2 * neighbor % 2][c]])
                    ):
                        if (
                            links[links[neighbor + 1 - 2 * neighbor % 2][c]][i]
                            == neighbor + 1 - 2 * neighbor % 2
                        ):
                            confirmationOfLinks[
                                links[neighbor + 1 - 2 * neighbor % 2][c]
                            ][i] += contact[2]
                            weightedconfirmationOfLinks[
                                links[neighbor + 1 - 2 * neighbor % 2][c]
                            ][i] += (contact[2] / coverage[contig1] / coverage[contig2])
                    for i in range(len(links[neighbor])):
                        if links[neighbor][i] == contig1 * 2:
                            confirmationOfLinks[neighbor][i] += contact[2]
                            weightedconfirmationOfLinks[neighbor][i] += (
                                contact[2] / coverage[contig1] / coverage[contig2]
                            )

        for j, neighbor in enumerate(links[contig1 * 2 + 1]):
            for j in range(len(links[contig1 * 2 + 1])):
                if (
                    links[contig1 * 2 + 1][j] == contig2 * 2
                    or links[contig1 * 2 + 1][j] == contig2 * 2 + 1
                ):
                    confirmationOfLinks[contig1 * 2 + 1][j] += contact[2]
                    weightedconfirmationOfLinks[contig1 * 2 + 1][j] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    for i in range(len(links[links[contig1 * 2 + 1][j]])):
                        if links[links[contig1 * 2 + 1][j]][i] == contig1 * 2 + 1:
                            confirmationOfLinks[links[contig1 * 2 + 1][j]][
                                i
                            ] += contact[2]
                            weightedconfirmationOfLinks[links[contig1 * 2 + 1][j]][
                                i
                            ] += (contact[2] / coverage[contig1] / coverage[contig2])

            # two connexions away
            otherEndOfNeighbor = neighbor + 1 - 2 * (neighbor % 2)

            for c in range(
                len(links[otherEndOfNeighbor])
            ):  # we take the other end of the neighbor contig
                if (
                    links[otherEndOfNeighbor][c] == contig2 * 2
                    or links[otherEndOfNeighbor][c] == contig2 * 2 + 1
                ):
                    confirmationOfLinks[contig1 * 2 + 1][j] += contact[2]
                    weightedconfirmationOfLinks[contig1 * 2 + 1][j] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    confirmationOfLinks[otherEndOfNeighbor][c] += contact[2]
                    weightedconfirmationOfLinks[otherEndOfNeighbor][c] += (
                        contact[2] / coverage[contig1] / coverage[contig2]
                    )

                    for i in range(len(links[links[otherEndOfNeighbor][c]])):

                        if links[links[otherEndOfNeighbor][c]][i] == otherEndOfNeighbor:
                            confirmationOfLinks[links[otherEndOfNeighbor][c]][
                                i
                            ] += contact[2]
                            weightedconfirmationOfLinks[links[otherEndOfNeighbor][c]][
                                i
                            ] += (contact[2] / coverage[contig1] / coverage[contig2])
                    for i in range(len(links[neighbor])):
                        if links[neighbor][i] == contig1 * 2:
                            confirmationOfLinks[neighbor][i] += contact[2]
                            weightedconfirmationOfLinks[neighbor][i] += (
                                contact[2] / coverage[contig1] / coverage[contig2]
                            )
                    # if verif != 2 :
                    #   print ('il y a un probleme')

    return confirmationOfLinks, weightedconfirmationOfLinks

# a function to delete small contigs made of repeated sequences that have no HiC contacts but tons of links
def crush_small_contigs(segments, interactionMatrix) :
    
    #first list all segments to be deleted
    for segment in segments :
        if segment.length < 5000 and segment in segment.links[0] : #this characterize a small repeated sequence
            if segment.HiCcoverage == 0 : # this segment is blind, can't do anything but crush it
                segment.ID = -segment.ID #segments marked with negative IDs should be deleted
                print(segment.names, ' is a small contig')
                
    print('Listed all segments to delete')
            
    #then delete all the links going toward unwanted segments (we don't want to reroute useless links)
    for segment in segments :
        for endOfSegment in range(2) :
            for n in range(len(segment.links[endOfSegment])-1, -1, -1) :
                if segment.links[endOfSegment][n].ID < 0 :
                    del segment.links[endOfSegment][n]
                    del segment.otherEndOfLinks[endOfSegment][n]
                    del segment.CIGARs[endOfSegment][n]                                        
    print('Deleted all links towards bad segments')
    
    #reroute all links going in and out of crushed segments
    # for segment in segments :
    #     if segment.ID < 0 :
    
    #         print('Crushing segment ', segment.names)
            
    #         #then make this segment disappear, by linking segments to the left of this contig to segment to the right :
    #         for l, leftneighbor in enumerate(segment.links[0]) :
    #             if leftneighbor.ID != segment.ID :
    #                 #print('Ln before : ', leftneighbor.ID, [se.ID for se in leftneighbor.links[segment.otherEndOfLinks[0][l]]])
    #                 leftneighbor.add_a_bunch_of_end_of_links(segment.otherEndOfLinks[0][l], segment.links[1], segment.otherEndOfLinks[1], ['0M' for i in segment.links[1]])
                    
    #         for r, rightneighbor in enumerate(segment.links[1]):
    #             if rightneighbor.ID != segment.ID :
    #                 #print('Rn before : ', [se.ID for se in rightneighbor.links[segment.otherEndOfLinks[1][l]]])
    #                 rightneighbor.add_a_bunch_of_end_of_links(segment.otherEndOfLinks[1][r], segment.links[0], segment.otherEndOfLinks[0], ['0M' for i in segment.links[0]])
    
    # print('Rerouted all links around segments to delete')
    #delete all segments that should be
    for se in range(len(segments)-1, -1, -1) :
        if segments[se].ID < 0 :
            del segments[se]
            
    print('There are ', len(segments), ' segments left')
                     
                
#function to freeze and sometimes solve small o-loops (works only if long reads are there)
def solve_small_loops(listOfSegments, names, repeats, lr_links, check_links) :
    
    for se in range(len(listOfSegments)) :
        
        segment = listOfSegments[se]
        if segment in segment.links[0] : #this is a small loop of length 0
            
            if segment.links[0].count(segment) == 1 : #this is a o-loop, let's flatten it
                
                if repeats != [] :
                    replications = 0
                    for contig in segment.names :
                        replications = max(replications, repeats[names[contig]])
                        
                    segment.flatten(replications)
                    #print('In solve_small_loops, flattening ', segment.names, segment.insideCIGARs)
                
                else :
                    segment.freeze(0)
                    segment.freeze(1)

        toRemove = []
        for n, neighbor in enumerate(segment.links[0]) : #trying to detect o-loops of length 1
            endOfLink = segment.otherEndOfLinks[0][n]
            index = s.find_this_link(neighbor, 1-endOfLink, segment.links[1], segment.otherEndOfLinks[1], warning = False) #returns -1 if it does not find anything
            
            if index != -1 : # then this is a o-loop of length 1
                    
                    if lr_links != [] :
                        cA0 = segment.names[0]
                        oA0 = (0 == segment.orientations[0])
                        cB0 = neighbor.names[-endOfLink]
                        oB0 = (neighbor.orientations[-endOfLink] == endOfLink)
                        
                        if not (cA0, oA0, cB0, oB0) in lr_links and not (cB0, oB0, cA0, oA0) in lr_links :
                            if (segment, 0, neighbor, endOfLink) not in toRemove and (neighbor, endOfLink, segment, 0) not in toRemove :
                                toRemove += [(segment, 0, neighbor, endOfLink)]
                        
                        endOfLink2 = segment.otherEndOfLinks[1][index]
                        cA1 = segment.names[-1]
                        oA1 = (1 == segment.orientations[-1])
                        cB1 = neighbor.names[-endOfLink2]
                        oB1 = (neighbor.orientations[-endOfLink2] == endOfLink2)
                        if not (cA1, oA1, cB1, oB1) in lr_links and not (cB1, oB1, cA1, oA1) in lr_links:
                            if (segment, 1, neighbor, endOfLink2) not in toRemove and  (neighbor, endOfLink2, segment, 1) not in toRemove:
                                toRemove += [(segment, 1, neighbor, endOfLink2)]
                        
        # print('3')
        # check_segments(listOfSegments)          
        
        for i in toRemove :
            # print('In o-loops : removing link from ', i[0].names, i[1], 'to ', i[2].names, i[3])
            # print('Links from ', i[0].names, i[1], ' : ', [j.names for j in i[0].links[i[1]]], i[0].otherEndOfLinks[i[1]])
            # print('Links from ', i[2].names, i[3], ' : ', [j.names for j in i[2].links[i[3]]], i[2].otherEndOfLinks[i[3]])
            i[0].remove_end_of_link(i[1], i[2], i[3])
            i[2].remove_end_of_link(i[3], i[0], i[1])
            
        # print('4')
        # check_segments(listOfSegments) 
                
                
            
def solve_l_loops(segments, lr_links): #l-loops occur when one end of a contig is in contact with both end of another contig or when a contig is linked to itself at one end
    
    for segment in segments :
        for endOfSegment in range(2) :
            toRemove = []
                    
            for n in range(len(segment.links[endOfSegment])-1) :
                    
                if segment.links[endOfSegment][n].ID == segment.links[endOfSegment][n+1].ID and  segment.otherEndOfLinks[endOfSegment][n] == segment.otherEndOfLinks[endOfSegment][n+1]: #the two links going toward the same contig are next to each other because links are sorted
                    
                    #here we have a l-loop
                    neighbor = segment.links[endOfSegment][n]
                    lenToRemove = len(toRemove)
                    
                    if lr_links != [] :
                        #let's check if the two links of the l-loop are confirmed by long reads
                        
                        if neighbor.ID == segment.ID : #here we have a l-loop of length 0
                            cA = segment.names[-endOfSegment]
                            oA = (endOfSegment == segment.orientations[-endOfSegment])
                            if not (cA, oA, cA, not oA) in lr_links :
                                toRemove += [(segment, endOfSegment, segment, endOfSegment)]
                            
                        else : #l-loop of length 1
                        
                            cA0 = segment.names[-endOfSegment]
                            oA0 = (endOfSegment == segment.orientations[-endOfSegment])
                            cB0 = neighbor.names[-segment.otherEndOfLinks[endOfSegment][n]]
                            oB0 = (neighbor.orientations[-segment.otherEndOfLinks[endOfSegment][n]] == segment.otherEndOfLinks[endOfSegment][n])
                            if not (cA0, oA0, cB0, oB0) in lr_links and not (cB0, not oB0, cA0, not oA0) in lr_links :
                                toRemove += [(segment, endOfSegment, neighbor, int(oB0))]
                                    
                            cA1 = segment.names[-endOfSegment]
                            oA1 = (endOfSegment == segment.orientations[-endOfSegment])
                            cB1 = neighbor.names[-segment.otherEndOfLinks[endOfSegment][n+1]]
                            oB1 = (neighbor.orientations[-segment.otherEndOfLinks[endOfSegment][n+1]] == segment.otherEndOfLinks[endOfSegment][n+1])
                            if not (cA1, oA1, cB1, oB1) in lr_links and not (cB1, not oB1, cA1, not oA1) in lr_links:
                                toRemove += [(segment, endOfSegment, neighbor, int(oB1))]
                    
                    if len(toRemove) == lenToRemove : #means that no links could be taken out
                        segment.freeze(endOfSegment)
                        
            for i in toRemove :
                i[0].remove_end_of_link(i[1], i[2], i[3])
                i[2].remove_end_of_link(i[3], i[0], i[1])

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





















