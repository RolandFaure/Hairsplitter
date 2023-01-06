#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 12:37:18 2021

@author: rfaure
"""

import re #to find all numbers in a mixed number/letters string (such as 31M1D4M), to split on several characters (<> in longReads_interactionMatrix)
from collections import Counter #to count the elements in a list quickly
from scipy import sparse #to handle sparse matrices
from determine_multiplicity import determine_multiplicity
from finish_untangling import merge_adjacent_contigs
from input_output import read_GAF
from input_output import read_TSV
import time
from copy import deepcopy

#from segment import find_this_link

import segment as sg

#Master function of the file
#Input : initial gfa (as a list of segments), a GAF file with long reads mapped to the segments, names (which is an index numbering the contig), multiplicity : the pre-computed ploidy of each contig (as numbered in names), exhaustive : if you want to delete all links not found in the .gaf
#Output : new gfa (as a list of segments) corrected with long reads, and modified copiesnumber (taking into account contigs that have been duplicated)
def bridge_with_long_reads(segments, names, copiesnumber, gafFile, supported_links2, multiplicities, exhaustive):
    
    supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file
    
    lines = []
    if '.gaf' in gafFile :
        print("Reading the gaf file...")
        read_GAF(gafFile, 0, 0, lines)
        print("Finished going through the gaf file.")
    elif '.tsv' in gafFile :
        print("Reading the tsv file...")
        read_TSV(gafFile, names, lines)
        print("Finished going through the tsv file.")
    else :
        print("ERROR: input format of mapped read not recognized. It should be .gfa or .gpa")
        sys.exit()
    
    #determine an approximate list of contigs that look haploid
    haploidContigs, haploidContigsNames = determine_haploid_contigs(lines, segments, names)
    #print(haploidContigsNames)
    sure_haploids = False
    
    #inventoriate all bridges in the list bridges : sequence of contigs found in the GAF containing at least one haploid contig. Do fill in the longContigs list while you're at it
    longContigs = [] #long contigs are contigs too long to be traversed by a long read in the gaf
    longContigs = [True for i in range(len(names))] #then all contigs that are in the middle of a read will be marked as False
    bridges = [[[],[]] for i in range(len(haploidContigs))] #bridges is a list inventoring at index haploidCOntigsNames[seg.names[0]] all the links left and right of the contig, supported by the gaf
    minimum_supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #minimum_supported links is the list of all links between different contigs found at least once in the gaf file
    inventoriate_bridges(lines, bridges, minimum_supported_links, haploidContigsNames, longContigs, names, segments) 
    
    # print("Bridge 1756 : ", bridges[haploidContigsNames['1756']])
    
    #if exhaustive, delete all links not found in the .gaf
    if exhaustive : 
        for s in segments :
            for end in range(2):
                index = 0
                while index < len(s.links[end]) :
                    if minimum_supported_links[2*names[s.names[0]]+end, 2*names[s.links[end][index].names[0]]+s.otherEndOfLinks[end][index]] == 0 :
                        sg.delete_link(s, end, s.links[end][index], s.otherEndOfLinks[end][index], warning = True)
                    else :
                        index += 1
    
    #now, from all the bridges, build consensus bridges
    consensus_bridges = [['',''] for i in range(len(haploidContigs))] #consensus bridge is essentially the same as bridges, except there is only one bridge left at each side for each contig
    print("Building consensus bridges from all the long reads")

    haploidContigs, haploidContigsNames, consensus_bridges = build_consensus_bridges(consensus_bridges, bridges, names, haploidContigs, haploidContigsNames)
    print("Done building consensus bridges                 ")

    # print("contig 1756 in haploid : ", "1756" in haploidContigsNames)
    # print("contig 93 consensus : ", consensus_bridges[haploidContigsNames["93"]])
    # print("contig 11 consensus : ", consensus_bridges[haploidContigsNames["11"]])
    # print("contig 2601 consensus : ", consensus_bridges[haploidContigsNames["2601"]])
    bridges = []
    
    # print("contig 121 in haploid : ", "121" in haploidContigsNames)

        
    print("Now we will determine through an iterative process what contigs of the assembly are present only once in the final genome")

    while not sure_haploids : #knowing what contigs are really haploid may take several iterations
        
        leng = len(haploidContigs)
        
       # print("consensus of 20: ", consensus_bridges[haploidContigsNames['20']])

        #consensus bridges overlap two by two (e.g. >2>3>4 right of 1 overlaps with <3<2<1 left of 4), so merge them to have a set of non-overlapping consensus bridges
        non_overlapping_bridges = [['',''] for i in range(len(haploidContigs))] 

        haploidContigs, haploidContigsNames, consensus_bridges, sure_haploids = merge_bridges(non_overlapping_bridges, consensus_bridges, haploidContigsNames, haploidContigs, longContigs, names, multiplicities)

        #print("nonoverlap consensus of 20: ", non_overlapping_bridges[haploidContigsNames['20']])
        #if this last phase of merge_contig detected no inconstistencies, sure_haploids=True and the program moves on
        
        if not sure_haploids :
            print("Out of ", leng, " supposed single-copy contigs, ", leng-len(haploidContigs), " were not actually haploid. Recomputing until all the single-copy contigs are robust")
        
    # print("contig 154 in haploid : ", "154" in haploidContigsNames)

    #from the consensus bridges, mark all links that are supported by the long reads
    supported_links = sparse.lil_matrix((len(names)*2, len(names)*2)) #supported links is the list of the links between different contigs found in the gaf file, and in how many different consensus
    compute_supported_links(supported_links, non_overlapping_bridges, haploidContigsNames, haploidContigs, longContigs, names, non_overlapping=True)
    
    supported_links = supported_links.tocoo()
    supported_links2 = supported_links2.todok()
    #now the haploid contigs are determined with confidence
    
    for r, c, m in zip(supported_links.row, supported_links.col, supported_links.data):
        supported_links2 [r,c] = max(supported_links2[r,c], m)  
        supported_links2 [c,r] = max(supported_links2[c,r], m)  
    #print("Link is supported after maxing with strength ", supported_links2[2*names['127'], 2*names['112']])
    # print("Link is supported after maxing with strength ", supported_links2[2*names['2505'], 2*names['2489']])
         
    # print("Contig 15 has for consensus : ", consensus_bridges[haploidContigsNames["15"]])

    #now actually unzip the graph using the instructions in non_overlapping_bridges
    

    print("Let's move on to actually untangling the graph")
    unzip_graph_with_bridges(segments, non_overlapping_bridges, copiesnumber, haploidContigs, haploidContigsNames, names, supported_links2.copy(), minimum_supported_links, multiplicities.copy(), longContigs)
    
    #now remove the tips that apparently came from an overestimation in the multiplicity
    print("Now we correct the last quirks by looking a posteriori at the graph               ")
    merge_adjacent_contigs(segments)
    trim_tips(segments, multiplicities, names, haploidContigsNames, supported_links2, exhaustive)
    
        
    #print(non_overlapping_bridges)
    
    return segments


#input : all aligned long reads
#output : A list of "haploid" contigs, i.e. contigs that have at most one possible other contig right and left AND THAT ARE LONG ENOUGH, BECAUSE SHORT CONTIGS ARE USELESS
def determine_haploid_contigs(lines, segments, names) :
    
    #haploidContigsIdx = set([i for i in range(len(segments))]) #list of the idx of all haploid contigs : we'll whittle it down
    neighborLeftRight = [({}, {}) for i in range(len(segments))] #list of neighbor contigs left and right of each contig : if a contig has only one neighbor left and right we'll say it's haploid
    
    for line in lines :
        
        contigs = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del contigs[0] #because the first element is always ''
    
        for c, contig in enumerate(contigs) :
            
            #if names[contig] in haploidContigsIdx :
                
                orientation = '><'.index(orientations[c])
                
                if c > 0 :
                    
                    if contigs[c-1] not in neighborLeftRight[names[contig]][orientation] :
                        neighborLeftRight[names[contig]][orientation][contigs[c-1]] = 1
                    else :
                        neighborLeftRight[names[contig]][orientation][contigs[c-1]] += 1
                    # if len (neighborLeftRight[names[contig]][orientation].keys()) > 1 :
                    #     #haploidContigsIdx.discard(names[contig])
      
                if c < len(contigs) -1 :
                    
                    if contigs[c+1] not in neighborLeftRight[names[contig]][1-orientation] :
                        neighborLeftRight[names[contig]][1-orientation][contigs[c+1]] = 1
                    else :
                        neighborLeftRight[names[contig]][1-orientation][contigs[c+1]] += 1
                        
                    # if len (neighborLeftRight[names[contig]][1-orientation].keys()) > 1 :
                    #     #haploidContigsIdx.discard(names[contig])

    
    haploidContigs = []
    for se, nei in enumerate(neighborLeftRight) :
        if segments[se].length > 100 : #that's because too short contigs cause trouble (in great part because erroneous contigs are often very short)
            #if (len(nei[0]) ==0 or max(nei[0].values()) > 0.9 * sum(nei[0].values())) and (len(nei[1])==0 or max(nei[1].values()) > 0.9 * sum(nei[1].values())) :
            if len(segments[se].links[0]) <= 1 and len(segments[se].links[1]) <= 1 :
                
                haploidContigs += [segments[se]]
        
       
    haploidContigs.sort(key= lambda x: x.length, reverse = True)
    
    haploidContigsNames = {} #contains the index of each contig (identified by its name) in the haploidContigs list
    index = 0
    for s in haploidContigs :
        haploidContigsNames[s.names[0]] = index
        index += 1
        
    return haploidContigs, haploidContigsNames
        
#input : a list of alignments of a gaf file
#output : the completed bridges list, with for each haploid contig a list of what was found left and right of the contig. 
def inventoriate_bridges(lines, bridges, minimum_supported_links, haploidContigsNames, longContigs, names, segments) :
    
    
    for l, line in enumerate(lines) :      
                    
        if (l+1) % 1000 == 0 :
            print("Inventoried ", l+1, " long reads over ", len(lines), end = '\r')
        
        contigs = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del contigs[0] #because the first element is always ''
        
        #first go through the alignment to make sure it is possible on the gfa
        possible = True
        for c, contig in enumerate(contigs) :
            if c>0 :
                or1 = '<>'.index(orientations[c-1])
                or2 = '><'.index(orientations[c])
                #check if the link actually exists (it should, if the aligner did its job correctly, but apparently sometimes SPAligner behaves strangely)
                if -1 == sg.find_this_link(segments[names[contig]], or2, segments[names[contigs[c-1]]].links[or1], segments[names[contigs[c-1]]].otherEndOfLinks[or1]) :
                    print ("WARNING: discrepancy between what's found in the alignment files and the inputted GFA graph. Link ", contigs[c-1:c+1], orientations[c-1:c+1], " not found in the gfa")
                    possible = False
                    
        #then, only inventoriate the bridge if it is possible with respect to the graph     
        if  possible :
    
            for c, contig in enumerate(contigs) :
                
                
                if c>0 :
                    
                    or1 = '<>'.index(orientations[c-1])
                    or2 = '><'.index(orientations[c])
                    
                    minimum_supported_links[2*names[contigs[c-1]] + or1, 2*names[contigs[c]] + or2] = 1
                    minimum_supported_links[2*names[contigs[c]] + or2, 2*names[contigs[c-1]] + or1] = 1
        
                if c > 0 and c < len(contigs) - 1 :
                    longContigs[names[contig]] = False
                                            
                if contig in haploidContigsNames :
                    
                    
                    if orientations[c] == ">" :
                        r = 0
                        #first look at what contigs are left of the contig of interest
                        bridges[haploidContigsNames[contig]][1] +=  [""]
                        for c2 in range(c+1, len(contigs)) :
                            
                            bridges[haploidContigsNames[contig]][1][-1] += orientations[c2] + contigs[c2]
                            
                            # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                            #     break
                            
                        #then look at what's left of the contig of interest (so mirror the orientations)
                        bridges[haploidContigsNames[contig]][0] +=  [""]
                        for c2 in range(c-1 , -1, -1) :
                            
                            if orientations[c2] == '>' :
                                bridges[haploidContigsNames[contig]][0][-1] += '<' + contigs[c2]
                            else :
                                bridges[haploidContigsNames[contig]][0][-1] += '>' + contigs[c2]
                                
                            # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                            #     break
                                
                    else :
                        
                        #first look at what contigs are left of the contig of interest
                        bridges[haploidContigsNames[contig]][0] +=  [""]
                        for c2 in range(c+1, len(contigs)) :
                            bridges[haploidContigsNames[contig]][0][-1] += orientations[c2] + contigs[c2]
                            # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                            #     break
                            
                        #then look at what's left of the contig of interest (so mirror the orientations)
                        bridges[haploidContigsNames[contig]][1] +=  [""]
                        for c2 in range(c-1 , -1, -1 ) :
                            
                            if orientations[c2] == '>' :
                                bridges[haploidContigsNames[contig]][1][-1] += '<' + contigs[c2]
                            else :
                                bridges[haploidContigsNames[contig]][1][-1] += '>' + contigs[c2]
                                
                            # if contigs[c2] in haploidContigsNames : #you can stop the bridge, you've reached the other side
                            #     break
                                                     
#input : list of bridges for each haploid contig
#output : completed consensus_bridges, where there is max one bridge at each end of contig
def build_consensus_bridges(consensus_bridges, bridges, names, haploidContigs, haploidContigsNames):
    
    not_actually_haploid = [] #a list of not actually haploid contigs : at this point, we will rule out as haploid contigs those that consensus back on themselves
    
    #print(allcontigs[0][1])
    
    for c in range(len(bridges)) :
                
        if (c)%100 == 0 :
            print("consensused ", c, " bridges out of ", len(consensus_bridges), end='\r')

            
        localContigs = [ [ re.split('[><]' , bridges[c][j][k])[1:] for k in range(len(bridges[c][j])) ] for j in range(2)]
        localOrientations = [ [ "".join(re.findall("[<>]", bridges[c][j][k])) for k in range(len(bridges[c][j])) ] for j in range(2)]
        
        for end in range(2) : #look left and right of each contig
        
            kept_reads = [i for i in range(len(bridges[c][end])) if bridges[c][end][i] != ""] #list keeping in mind all the reads still in the consensus
            pos = 0
            consensus = ''
            

            while len(kept_reads) > 0 :
                

                candidate2 = [localContigs[end][i][pos] for i in kept_reads]
                candidate1 = [localOrientations[end][i][pos] for i in kept_reads]
                
                
                cons1 = Counter(candidate1).most_common(1)[0]
                cons2 = Counter(candidate2).most_common(1)[0]
                
                # if ('11' in haploidContigs[c].names) :
                #     print("Consensuing : ", Counter(candidate2))
                                
                if cons1[1]>0 and (cons1[1] > 0.6*len(candidate1) or (cons1[1]>=2 and cons1[1] == len(candidate1)-1)) and (cons2[1] > 0.6*len(candidate1) or (cons2[1]>=2 and cons2[1] == len(candidate1)-1)) : #consensus is defined there as a 60% of the propositions agreeing or only 1 proposition disagreeing
                
                    if cons1[0] == '*' :
                        break
                    
                    consensus_bridges[c][end] = consensus_bridges[c][end] + cons1[0] + cons2[0]
                    
                    # if cons2[0] == haploidContigs[c].names[0] and len(consensus_bridges[c][end]) < len(names): #if the contig loops on itself, it is not haploid except if we're on a circular chromosome
                    #     #not_actually_haploid += [c]
                    #     #break
                    #     pass
                    
                    
                    #inventoriate this bridge in supported_links
                    
                    #print('Adding from ', consensus_bridges[c][end], " c: ", c)
                    
                    #update the list of reads still there to build the consensu
                    new_kept = []
                    for r in kept_reads:
                        if len(localContigs[end][r]) > pos+1 :
                            if localOrientations[end][r][pos] == cons1[0] and localContigs[end][r][pos] == cons2[0] :
                                new_kept += [r]
                        
                    kept_reads = new_kept
                    
                    # if cons2[0] in haploidContigsNames : #just go up to the next haploid contig
                    #     break
               
                else : #if there is no consensus let's stop
                    # if "11" in haploidContigs[c].names :
                    #     print("There are no consensus there: ", candidate1, candidate2)
                    break
                
                pos += 1
        bridges[c] = [] #to free memory
        
    new_consensus_bridges = []
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    not_actually_haploid = list(set(not_actually_haploid))
    not_actually_haploid.sort()
        
    index = 0
    indexNot = 0
    
    for i in range(len(haploidContigs)) :
        
        if indexNot<len(not_actually_haploid) and i == not_actually_haploid[indexNot] :
            indexNot += 1
        else :
            new_consensus_bridges += [consensus_bridges[i]]
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].names[0]] = index
            index += 1
                    
    return reliable_haploid_contigs, reliable_haploid_contigsNames, new_consensus_bridges
                
#input: list of consensus bridges. non_overlapping is an argument that is True when the consensus bridges are non-overlapping
#output: filled supported_links matrix
def compute_supported_links(supported_links, consensus_bridges, haploidContigsNames,haploidContigs, longContigs, names, non_overlapping = False) :
    
    symmetricSupport = deepcopy(supported_links)
    for c in range(len(consensus_bridges)) :
        
        for end in range(2) :
            
            
            contigs = re.split('[><]' , consensus_bridges[c][end])
            orientations = "".join(re.findall("[<>]", consensus_bridges[c][end]))
            del contigs[0] #because the first element is always ''
            
            previousContig = names[haploidContigs[c].names[0]] # variable containing the index of the last contig of the consensus (useful to fill supported_links)
            previousEnd = end #variable containing the free end of the last contig of the consensus (useful to fill supported_links)
            
            firstHapIdx = -1 #to see if there is another end to the bridge
            for co, contig in enumerate(contigs):
                if contig in haploidContigsNames :
                    firstHapIdx = co
                    break
            
            if len(contigs) > 0 :
                if longContigs[names[contigs[-1]]] : #it's a long contig, a bridge will be built
                    firstHapIdx = len(contigs)-1
                
     
            # if '7200' in contigs :
            #     print("Bridge from 433 will look like : ", contigs[:firstHapIdx+1], " ", firstHapIdx)
     
            for co in range(firstHapIdx+1):

                contig = contigs[co]
                # if c == haploidContigsNames['83467'] or c == haploidContigsNames['55374'] :
             

                current_end = 0
                if orientations[co] == "<" :
                    current_end=1
                    
                supported_links[2*names[contig]+current_end, 2*previousContig+previousEnd] += 1  #each link is present twice in supported_links, only add to one cell of supported_links, the symetric one can be added from the other end of the bridge
                symmetricSupport[2*previousContig+previousEnd, 2*names[contig]+current_end] += 1
                # if (contig == "7200" and names['4796'] == previousContig) or (contig == "4796" and names['7200'] == previousContig) :
                #     print("Increasing the link : ", haploidContigs[c].names[0], ",", consensus_bridges[c][end])
                if non_overlapping :
                    symmetricSupport[2*names[contig]+current_end, 2*previousContig+previousEnd] += 1
                    supported_links[2*previousContig+previousEnd, 2*names[contig]+current_end] += 1
                      
                previousContig = names[contig]
                previousEnd = 1-current_end
                
                if contig in haploidContigsNames :
                    break
      
    if not non_overlapping :                   
        supported_links = supported_links.todok()
        symmetricSupport = symmetricSupport.tocoo()
        #now the haploid contigs are determined with confidence
        
        for r, c, m in zip(symmetricSupport.row, symmetricSupport.col, symmetricSupport.data):
            supported_links[r,c] = max(supported_links[r,c], m)  
            
    # print("Strength of link : ", supported_links[2*names['4796'], 2*names['7200']])

                                                        
#input: list of all consensus bridges for all haploid contigs, plus a list of the haploid contigs
#output: a list of non-overlapping bridges left and right of each contig, with only full bridges (connecting 2 haploid contigs)
def merge_bridges(non_overlapping_bridges, consensus_bridges, haploidContigsNames, haploidContigs, longContigs, names, multiplicities) :
    
    sure_haploids = True
    not_actually_haploid = [] #a list of not actually haploid contigs among the haploidContigs
    
    #create dict firstHapIndices, indicating for each end the index of first haploid contig, -1 if there is none
    firstHapIndices = {}
    for c in range(len(consensus_bridges)) :
        for end in range(2) :
            
            contigs = re.split('[><]' , consensus_bridges[c][end])
            orientations = "".join(re.findall("[<>]", consensus_bridges[c][end]))
            del contigs[0] #because the first element is always ''
            
            firstHapIndices[2*c+end] = -1
            for co in range(len(contigs)):
                if contigs[co] in haploidContigsNames :
                    firstHapIndices[2*c+end] = co
                    break
            
    
    for c in range(len(consensus_bridges)) :
        
        for end in range(2) :
            
            if consensus_bridges[c][end] != '' :
            
                contigs = re.split('[><]' , consensus_bridges[c][end])
                orientations = "".join(re.findall("[<>]", consensus_bridges[c][end]))
                del contigs[0] #because the first element is always ''

                firstHapIdx = firstHapIndices[2*c+end]
                
                if firstHapIdx >= 0 : #in case there is a full bridge
                    #check if the two bridges are coherent
                    coherent = False
                    otherEnd = 0
                    if orientations[firstHapIdx] == "<":
                        otherEnd = 1
                        
                        # print ("here is contigs [0] : ", haploidContigs[c].names[0], " -- ", re.split('[><]' , consensus_bridges[haploidContigsNames[contigs[-1]]][otherEnd]))
                    
                    symmetrical = re.split('[><]' , consensus_bridges[haploidContigsNames[contigs[firstHapIdx]]][otherEnd])
                    
                    if len(symmetrical)>firstHapIdx+1 and haploidContigs[c].names[0] == symmetrical[firstHapIdx+1] : #+1 because the first element of symmetrical is always ''
                        coherent = True
                        
                    #if the symmetrical bridge does not lead anywhere, take this one as the truth
                    # if firstHapIndices[2*haploidContigsNames[contigs[firstHapIdx]]+otherEnd] <= -1 :
                    #     non_overlapping_bridges[c][otherEnd] = ''.join([orientations[i]+contigs[i] for i in range(firstHapIdx+1)])
                        
                    #     else : #if it's not the first time, there is a problem with the ploidy of the contigs
                    #         non_overlapping_bridges[c][otherEnd] = ''
                                                                        
                        
                                                
                    if not coherent :
                        
                        #check if the incoherence is not due to being on a very weakly supported (and probably false) path
                        
                        if len(symmetrical)<=firstHapIdx+1 : #means the symmetrical did not manage to consensus
                            # if '127' == contigs[firstHapIdx] :
                            #     print("two incoherent bridges between ", haploidContigs[c].names[0] , " and ", haploidContigs[haploidContigsNames[contigs[firstHapIdx]]].names[0], " : ", consensus_bridges[c][end], " vs ", consensus_bridges[haploidContigsNames[contigs[firstHapIdx]]][otherEnd])
                            not_actually_haploid += [haploidContigsNames[contigs[firstHapIdx]]]
                            sure_haploids = False
                            
                            # print("Contig ", contigs[firstHapIdx], " does not look haploid, seen from contig ",  haploidContigs[c].names[0], " : ", contigs[:firstHapIdx+1], ", whose bridge is ", symmetrical[:firstHapIdx+2])

                        # else : #if the path is really minoritary, the contig is probably not haploid, or at least you should not count on it
                        #     not_actually_haploid += [c]
                        #     sure_haploids = False
                        #     if haploidContigs[c].names[0] == '7' :
                        #         print("Contig ", haploidContigs[c].names[0], " does not look haploid, here is its bridge to the next contig ", contigs[firstHapIdx] , " : ", contigs[:firstHapIdx+1], ", whose bridge is ", symmetrical[:firstHapIdx+2])
                        #     #print("Haploid contigs : ", haploidContigsNames)
                            
                        
                        
                    else :
                        if c < haploidContigsNames[contigs[firstHapIdx]] : #keep only one of the two bridges
                            non_overlapping_bridges[c][end] = ''.join([orientations[i]+contigs[i] for i in range(firstHapIdx+1)])
                
                elif longContigs[names[contigs[-1]]] : #if it ends on a  multiploid long contig, extracting could still be done (though not for the last, long contig). That's because we know there is no bridge going the other direction through the long contig
                    non_overlapping_bridges[c][end] = consensus_bridges[c][end]
                
                else : #a half-bridge cannot be unambiguously extacted from the graph
                    non_overlapping_bridges[c][end] = ''
    
    #now update haploidCOntig and consensus_bridges
        
    new_consensus_bridges = []
    reliable_haploid_contigs = []
    reliable_haploid_contigsNames = {}
    not_actually_haploid = list(set(not_actually_haploid))
    not_actually_haploid.sort()
        
    index = 0
    indexNot = 0
    
    for i in range(len(haploidContigs)) :
        
        if indexNot<len(not_actually_haploid) and i == not_actually_haploid[indexNot] :
            indexNot += 1
            # if multiplicities[names[haploidContigs[i].names[0]]] == 1 :
            #     multiplicities[names[haploidContigs[i].names[0]]] += 1
        else :
            new_consensus_bridges += [consensus_bridges[i]]
            reliable_haploid_contigs += [haploidContigs[i]]
            reliable_haploid_contigsNames[haploidContigs[i].names[0]] = index
            index += 1
                    
    return reliable_haploid_contigs, reliable_haploid_contigsNames, new_consensus_bridges, sure_haploids #return the updated list of haploid contigs
          
                
#input : a list of segments and the non_overlapping_bridges
#output: a list of segment where all the bridges have been built          
def unzip_graph_with_bridges(segments, non_overlapping_bridges, copiesnumber, haploidContigs, haploidContigsNames, names, supported_links, minimum_supported_links, multiplicities, longContigs) :  
    
    #compute the minimum multiplicity of each contig, so that there are enough contigs to build all bridges, even when the depth suggests otherwise
    #first check with supported links
    for s, seg in enumerate(segments) :
        minLeft = 0
        minRight = 0
        for n, neighbor in enumerate(seg.links[0]) :
            minLeft += supported_links[2*names[seg.names[0]], 2*names[neighbor.names[0]]+seg.otherEndOfLinks[0][n]]
        for n, neighbor in enumerate(seg.links[1]) :
            minRight += supported_links[2*names[seg.names[0]]+1 , 2*names[neighbor.names[0]]+seg.otherEndOfLinks[1][n]]
        
        multiplicities[s] = max( min (minLeft, minRight), multiplicities[s])
        
    #second check with non_overlapping_bridges
    minimum_multiplicity = [0 for i in range(len(multiplicities))] 
    for c in range(len(non_overlapping_bridges)) :
        multiplicities[names[haploidContigs[c].names[0]]] += 1
        for end in range(2) :
            contigs = re.split('[><]' , non_overlapping_bridges[c][end])
            del contigs[0] #because the first element is always ''
            for contig in contigs :
                minimum_multiplicity[names[contig]] += 1
    for s in range(len(segments)) :
        multiplicities[s] = max(minimum_multiplicity[s], multiplicities[s])
        
    # print("Multiplicities of contig 127 : ", multiplicities[names['127']])
    
            
    alreadyDuplicated = [-1 for i in range(len(names))] #list useful for duplicating long contigs : only duplicate them from one side, the one that is in this list 
            
    l = len(segments)
    for se in range(l) :
        
        if (se)%100 == 0 :
            print("Processed ", se, " contigs out of ", l, ", while untangling with long reads", end = '\r')
        s = segments[se]
         
        if s.names[0] in haploidContigsNames :
            
                        
            for end in range(2) :
                
                
                if len(non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]) > 0  :#and ('358' in non_overlapping_bridges[haploidContigsNames[s.names[0]]][end] or '339' in non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]): #means there is a bridge to be built there
                                                  
                    contigs = re.split('[><]' , non_overlapping_bridges[haploidContigsNames[s.names[0]]][end])
                    orientations = "<>"[end] + "".join(re.findall("[<>]", non_overlapping_bridges[haploidContigsNames[s.names[0]]][end]))
                    del contigs[0] #because the first element is always ''
                    contigs = [s.names[0]] + contigs
                    contigsID = [segments[names[c]].ID for c in contigs]
                    newContigsIndices = [names[contigs[0]]] #list keeping track of the indices of the contigs on the path
                    oldContigsIndices = [names[i] for i in contigs]
                        
                    # if '6897' in contigs :
                    #     print("Unzipping contig ", s.names[0], " with bridge : ",contigs)
                    
                    haploidCoverage = (segments[names[contigs[0]]].depths[0]+segments[names[contigs[-1]]].depths[0])/2 #computation of the haploid coverage at this point in the assembly
                    
                    nextEnd = 0
                    if orientations[1] == '<' :
                        nextEnd = 1
                    CIGAR = segments[names[contigs[0]]].CIGARs[end][sg.find_this_link(segments[oldContigsIndices[1]], nextEnd, segments[names[contigs[0]]].links[end], segments[names[contigs[0]]].otherEndOfLinks[end])]
                    nextCIGAR = '';

                    
                    #take care of the first contig, nobody is going to do it elsewhise
                    # contig = segments[names[contigs[0]]]
                    # idx = 0
                    # if len(contig.links[end]) > 0 : #delete the good link, since it will be reestablished later (we don't want a double)
                    #     #print('here edge 291 neighbors are ', [i.names for i in segments[names['edge_291']].links[0]], " ", contigs)
                    #     while len(contig.links[end]) > idx :
                    #         neighbor = contig.links[end][idx]
                    #         if ((neighbor.names[0] != contigs[1] and supported_links[names[contigs[0]]*2+end , names[contigs[1]]*2+nextEnd] == 0) \
                    #             or (neighbor.names[0] == contigs[1] and neighbor.names[0] in haploidContigsNames)) \
                    #             and len(neighbor.links[contig.otherEndOfLinks[end][idx]]) > 1: #make sure not to create dead ends
                    #             success = sg.delete_link(contig, end, neighbor, contig.otherEndOfLinks[end][idx])
                    #         else :
                    #             idx += 1
                    

                        #print('now edge 291 neighbors are ', [i.names for i in segments[names['edge_291']].links[0]], " ", contigs)

                    #then take care of all other contigs
                    for c in range(1, len(contigs)) :
                        
                        contig = segments[names[contigs[c]]]
                        
                        
                        #multiplicity = max([1, round(contig.depths[0]/haploidCoverage), multiplicities[names[contigs[c]]]]) #the multiplicity is inferred from the coverage of the contig compared to the coverage of the two haploid contigs at the extremity of the bridge
                        multiplicity = multiplicities[names[contigs[c]]]
                        multiplicities[names[contigs[c]]] -= 1
                        # if contigs[c] == '7200' :
                        #     print("Link supported with strength : ", supported_links[2*names["7200"], 2*names["4796"]])
                            # print("Using once contig 7200 in overlap ", non_overlapping_bridges[haploidContigsNames[s.names[0]]][end], ", now at : ", multiplicities[names[contigs[c]]])
                        #neighborMultiplicity = max([1, round(segments[oldContigsIndices[c-1]].depths[0]/haploidCoverage), multiplicities[names[contigs[c-1]]]])
                        neighborMultiplicity = multiplicities[names[contigs[c-1]]]

                        #print(alreadyDuplicated[names['edge_291']])
                        
                        end1, end0 = 1, 1
                        if orientations[c] == '>' :
                            end1 = 0
                        if orientations[c-1] == '<' :
                            end0 = 0
                            
                        #remember the CIGAR before deleting all the links
                        if c < len(contigs)-1 :
                            nextEnd = 0
                            if orientations[c+1] == '<' :
                                nextEnd = 1
                            if sg.find_this_link(segments[oldContigsIndices[c+1]], nextEnd, contig.links[1-end1], contig.otherEndOfLinks[1-end1]) != -1 :
                                nextCIGAR = contig.CIGARs[1-end1][sg.find_this_link(segments[oldContigsIndices[c+1]], nextEnd, contig.links[1-end1], contig.otherEndOfLinks[1-end1])]
                            else :
                                print("Debug WARNING, ", contigs, " : looking for ", segments[oldContigsIndices[c+1]].names, " ", nextEnd, " from ", contig.names, " among ", [i.names for i in contig.links[1-end1]], " ", contig.otherEndOfLinks[1-end1], " ", s.names[0], " ",non_overlapping_bridges[haploidContigsNames[s.names[0]]][end])
                            
                        #print("Link supported with strength ", supported_links[names[contigs[c]*2+end1 , names[contigs[c-1]]*2+end0])
                        
                        # if s.names[0] == '55374' and end == 0 and contigs[c] == '2505' :
                        #     print("multiplicity of 2505 : " , multiplicity)
                                                    
                        if multiplicity > 1 and c < len(contigs)-1 : #if multiplicity>1, the contig should be duplicated 
                        
                            newSegment = sg.Segment(contig.names, contig.orientations, contig.lengths, contig.insideCIGARs, HiCcoverage = contig.HiCcoverage, readCoverage = [i/multiplicity for i in contig.depths])
                            segments.append(newSegment)
                            newContigsIndices += [len(segments) - 1]
                            
                            # if '4796' in newSegment.names : 
                            #     print ("Hrer I duplicate, ", contigs)
                                # print (supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0])
                            
                            for n in newSegment.names :
                                copiesnumber[n] += 1

                            #add the link to form the new bridge
                            sg.add_link(segments[-1], end1, segments[newContigsIndices[c-1]], end0, CIGAR)

                            #delete the old link if and only if it was only supported by one path only
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            supported_links[names[contigs[c-1]]*2+end0, names[contigs[c]]*2+end1] -= 1
                            minimum_supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            minimum_supported_links[names[contigs[c-1]]*2+end0, names[contigs[c]]*2+end1] -= 1
                                
                            if supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] == 0 and len(segments[oldContigsIndices[c-1]].links[end0])>1\
                                and sg.find_this_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]].links[end0], segments[oldContigsIndices[c-1]].otherEndOfLinks[end0], warning=False) != -1\
                                and len(segments[oldContigsIndices[c]].links[end1]) > 1: #be sure not to create dead ends
                                # if "130" in segments[oldContigsIndices[c-1]].names : 
                                #     print("remove links left of ", segments[oldContigsIndices[c]].names, " : ", segments[oldContigsIndices[c-1]].names)
                              sg.delete_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]], end0, warning = True) #though it is very hard to be sure, it is not impossible at that point that we actually delete a link that is present in another bridge, so don't warn
                                
                            #since contig has been duplicated, lower its depth
                            contig.divide_depths(multiplicity/(multiplicity-1))

                            
                            # if "3682" in newSegment.names :
                            #     print("New 3682 : ", non_overlapping_bridges[haploidContigsNames[s.names[0]]][end])
                            #     print([s.names for s in segments[-1].links[0]])
                            #     print([s.names for s in segments[-1].links[1]])
                                    

                        elif c == len(contigs)-1 : #we're at the end of the bridge
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            minimum_supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            
                            if supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] == 0 and len(segments[oldContigsIndices[c-1]].links[end0])>1 \
                                and sg.find_this_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]].links[end0], segments[oldContigsIndices[c-1]].otherEndOfLinks[end0], warning=False) != -1\
                                and len(segments[oldContigsIndices[c-1]].links[end0]) > 1: #be sure not to create dead ends
                                sg.delete_link(segments[oldContigsIndices[c]], end1, segments[oldContigsIndices[c-1]], end0, warning = True) 
  
                            if segments[newContigsIndices[c-1]] not in contig.links[end1] : #the link could have been already kept to prevent dead ends
                                sg.add_link(contig, end1, segments[newContigsIndices[c-1]], end0, CIGAR)
                            
                                 
                        else : #i.e. contig has multiplicity 1 and is not a longContig and is not the endpoint of the bridge
                            

                            # if the contig has multiple incoming reads left, get rid of all of them (the good one is reestablished later). Do the same at the other end of the contig
                            supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1
                            minimum_supported_links[names[contigs[c]]*2+end1 , names[contigs[c-1]]*2+end0] -= 1       
                            
                            if len(contig.links[end1]) > 0 :
                                for n, neighbor in enumerate(contig.links[end1]) :
                                    if neighbor.ID != segments[newContigsIndices[c-1]].ID:
                                        #then delete the links, except if they are still supported or create dead ends
                                        if supported_links[names[contigs[c]]*2+end1 , names[neighbor.names[0]]*2+contig.otherEndOfLinks[end1][n]] <= 0 \
                                            and minimum_supported_links[names[contigs[c]]*2+end1 , names[neighbor.names[0]]*2+contig.otherEndOfLinks[end1][n]] <= 0\
                                            and len(neighbor.links[contig.otherEndOfLinks[end1][n]]) > 1 :
                                            sg.delete_link(contig, end1, neighbor, contig.otherEndOfLinks[end1][n])
                                        
                            # print("RIght of 138 there is ", [i.ID for i in segments[names['138']].links[1]], " ", [i for i in segments[names['119']].otherEndOfLinks[1]], " ", segments[names['119']].ID)
                            # print("RIght of 119 there is ", [i.ID for i in segments[names['119']].links[1]], " ", [i for i in segments[names['138']].otherEndOfLinks[1]], " ", segments[names['138']].ID)
                            if len(contig.links[1-end1]) > 0 and c < len(contigs)-1 : 
                                
                                nextEnd = 2
                                for n, neighbor in enumerate(contig.links[1-end1]) :
                                    if neighbor.names[0] == contigs[c+1] :
                                        nextEnd = contig.otherEndOfLinks[1-end1][n]
                                
                                if not (c == len(contigs)-2 and alreadyDuplicated[names[contigs[c+1]]] == 1-nextEnd and contigs[c+1] not in haploidContigs): #delete all the links right of the contig, the only good one will be reestablished later, except if looking at a long contig duplicated from the other end
                                    idx = 0
                                    while len(contig.links[1-end1]) > idx :
                                        neighbor = contig.links[1-end1][idx]
                                        if len(neighbor.links[contig.otherEndOfLinks[1-end1][idx]]) > 1 : #you don't want to create dead ends
                                            sg.delete_link(contig, 1-end1, neighbor, contig.otherEndOfLinks[1-end1][idx])
                                        else :
                                            idx += 1

                            
                            #now the link the contig to the contig right at its left, except if the link had already been kept to prevent making dead ends
                            
                            if segments[newContigsIndices[c-1]] not in contig.links[end1] :
                                sg.add_link(contig, end1, segments[newContigsIndices[c-1]], end0, CIGAR) 
                            
                            newContigsIndices += [oldContigsIndices[c]]
                            
                        CIGAR = nextCIGAR
                    
                        
                    #print('edge 291 neighbors are ', [i.names for i in segments[names['edge_291']].links[0]], " ", contigs)
                # if supported_links[names['2484']*2+1 , names['2505']*2] == 1:
                #     print(non_overlapping_bridges[haploidContigsNames[s.names[0]]][end], "  ", [i.names for i in segments[names['2484']].links[1]], " ", s.names[0], " ", end)
                #     while True :
                #         rien = 0

      
#input : the graph of segments as well as the list of multiplicities
#output : the graph where the dubious tips have been deleted (those where multiplicity has been overestimated). Segments that have no supported links left and right are disconnected (which they weren't before because of dead ends creation)
def trim_tips(segments, multiplicities, names, haploidContigsNames, supported_links, exhaustive):
    
    # detach contigs that are totally unsupported (sometimes misassembled contigs are connected to the graph but totally unsupported : they remain there because we did not want to create dead ends)
    if exhaustive :
        for s, seg in enumerate(segments):
            
            if len(seg.names) == 1 : #our only target here is misassembled contigs
                someSupported = False
                for end in range(2) :
                    
                    # if multiplicities[s] < len(seg.links[end]):
                    
                    #if some links are supported and some aren't, delete those who aren't
                    for n, neighbor in enumerate(seg.links[end]) :
                        
                        otherEnd = seg.otherEndOfLinks[end][n]
                        
                        if supported_links[2*names[seg.names[0]]+end, 2*names[neighbor.names[0]]+otherEnd] == 1 : #this is supported
                            someSupported = True
                        if len(neighbor.links[otherEnd]) == 1 : #let's not create extra dead ends
                            someSupported = True
                    
                    
                if not someSupported :
                    for end in range(2) :
                        for n, neighbor in enumerate(seg.links[end]) :
                            otherEnd = seg.otherEndOfLinks[end][n]
                            sg.delete_link(seg, end, neighbor, otherEnd)
                        
                        
    #now trim tips, i.e. delete suspicious dead ends
    toDelete = []
    for s, seg in enumerate(segments):
        
        for end in range(2) :
            
            if len(seg.links[1-end]) == 0 and len(seg.links[end]) == 1 and seg.length < 2000 : #if this contig is a very short dead end, delete it

                neighbor = seg.links[end][0]
                neighborEnd = seg.otherEndOfLinks[end][0]
                #print("Checking if ", seg.names, " is a tip")
                if any([extended_length(i, neighbor.otherEndOfLinks[neighborEnd][e], 10*seg.length, 30) for e,i in enumerate(seg.links[end][0].links[seg.otherEndOfLinks[end][0]])]) : #this means we're in a very short dead end
                                    
                    if all([i not in haploidContigsNames for i in seg.names]) : #then it means it's probably an error in the determination of the multiplicity
                    
                        sg.delete_link(seg, end, seg.links[end][0], seg.otherEndOfLinks[end][0])
                        toDelete += [s]
                 
            if len(seg.links[1-end]) == 0 and len(seg.links[end]) == 1 : #if this contig is a replica of a parallel contig, merge both together
            
                    neighbor = seg.links[end][0]
                    neighborEnd = seg.otherEndOfLinks[end][0]
                                    
                    for n2, neighbor2 in enumerate(neighbor.links[neighborEnd]) :
                        
                        otherEnd2 = neighbor.otherEndOfLinks[neighborEnd][n2]
                        if neighbor2.length > seg.length : #if we're not looking at the same contig
                        
                            #see if it's parallel to the dead end, i.e. has the same contigs in the same order
                            same = True
                            for co in range(len(seg.names)) :
                                nameParallel1 = neighbor2.names[ otherEnd2*(len(neighbor2.names)-1) + (-2*otherEnd2+1)*co ]
                                nameParallel2 = seg.names[ end*(len(seg.names)-1) + (-2*end+1)*co ]
                                if nameParallel1 != nameParallel2 :
                                    same = False
                                    break
                                
                            if same :
                                # print("Found two parallels contigs : ", seg.names, neighbor2.names)
                                neighbor2.multiply_end_depths(2, otherEnd2, len(seg.names))
                                sg.delete_link(seg, end, seg.links[end][0], seg.otherEndOfLinks[end][0])
                                toDelete += [s]
        
                        
    
    for i in toDelete[::-1]:
        del segments[i]

#a function returning True if you can go far (up to threshold) with neighbors of neighbors of neighbors of... 
#it returns False if it needs to recur deeper than thresholdContigs (even though it might be true)
def extended_length(segment, end, thresholdLength, thresholdContigs) :
    
    #print("Extended length called with threshold ", thresholdLength, " on segment , ", sg.names)
    
    if thresholdContigs == 0 :
        return False
    
    if segment.length > thresholdLength :
        return True
    
    #start by looking down the longest contig, it will be fastest
    longestContig = [i for i in range(len(segment.links[1-end]))]
    longestContig.sort(key= lambda x : segment.links[1-end][x].length, reverse = True)
    
    #print("Longest contigs : ", longestContig, [segment.links[1-end][i].length for i in longestContig])
    
    for n in longestContig[:min(len(longestContig), 2)] : #only explore the 2 most promising neighbors, beyond it's not worth it
        
        neighbor = segment.links[1-end][n]
        if extended_length(neighbor, segment.otherEndOfLinks[1-end][n], thresholdLength-segment.length, thresholdContigs-1) :
            return True
        
    return False
        
        
        
        
        
        
        
        
        
        
        
        
        