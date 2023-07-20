#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 13:02:36 2021

@author: rfaure
"""

#main function of the file : tries to estimate how many copies of each contig is actually present in the actual assembly, based on the topology of the graph and the coverage
#input : a gfa (as a list of segments), with mandatory coverage information ; names (to know at what index to put each contig)
#output : computed_multiplicity, which is a list containing the theoretical multiplicity of contig 'a' at position names[a], as well as updated supported_links, telling which links actually exist
def determine_multiplicity(segments, names, supported_links, reliable_coverage=True) :

    computed_multiplicity = [0 for i in range(len(segments))]
    
    #first compute the 'reference haploid coverage', the average coverage for all reads with exactly one neighbors at each end
    refCoverages = 0
    weightedNumberOfRefContigs = 1 #1 and not 0 to be sure not to divide by 0
    for s in segments :
        
        links = s.get_links()
        if len(links[0]) <= 1 and len(links[1]) <= 1 : #if the contig has less than 1 neighbor at each end we can suppose it is haploid
        
            weightedNumberOfRefContigs += s.length
            refCoverages += s.length * s.depths[0]

    #print(len(segments), refCoverages, weightedNumberOfRefContigs)
    refCoverage = refCoverages / weightedNumberOfRefContigs
    
    if (refCoverage == 1 or refCoverage == 0) and reliable_coverage :
        reliable_coverage = False
        
    if not reliable_coverage :
        refCoverage = 1
        
    #then inventoriate all haploid contigs
    for s in segments :
        
        links = s.get_links()
        if len(links[0]) <= 1 and len(links[1]) <= 1 and (round(s.depths[0] / refCoverage) <= 1 or refCoverage == 1): #if the contig has 1 neighbor at each end we can suppose it is haploid
        
            computed_multiplicity[names[s.names[0]]] = 1
                   
    #now infer greedily the coverage of contigs next to haploid ones

    s = 0 #s is the index of the segment we're looking at
    unchanged = 0
    while unchanged < len(segments) :
        
        if computed_multiplicity[s] == 0 :
            
            new_multiplicity = 0
            new_multiplicity1 = 0
            new_multiplicity2 = 0
            confidence = False
            if len(segments[s].links[0]) > 0 and all([computed_multiplicity[names[i.names[0]]] > 0 for i in segments[s].links[0]]) and all([len(i.links[segments[s].otherEndOfLinks[0][j]]) ==1 for j,i in enumerate(segments[s].links[0])]) :
                new_multiplicity1 = sum([computed_multiplicity[names[i.names[0]]] for i in segments[s].links[0]])
                        
            if len(segments[s].links[1]) > 0 and all([computed_multiplicity[names[i.names[0]]] > 0 for i in segments[s].links[1]]) and all([len(i.links[segments[s].otherEndOfLinks[1][j]]) ==1 for j,i in enumerate(segments[s].links[1])]) :
                new_multiplicity2 = sum([computed_multiplicity[names[i.names[0]]] for i in segments[s].links[1]])

            if new_multiplicity1 == new_multiplicity2 : #in this case the confidence is very high
                new_multiplicity = new_multiplicity1
                confidence = True
            else :
                new_multiplicity = max(new_multiplicity1, new_multiplicity2) #that includes the case where both are > 0, in which case we are not confident
                
            # if 'edge_224' in segments[s].names :
            #      print("Looking at 224 and ", new_multiplicity)
            if new_multiplicity > 0 and (segments[s].depths[0]/refCoverage > new_multiplicity/1.5 or confidence or refCoverage == 1):
                    computed_multiplicity[s] = new_multiplicity
                    #print("coucou, ", segments[s].names, " ", new_multiplicity)
                    unchanged = -1
            
            #also take care of the supported_links
            if new_multiplicity1 == new_multiplicity and new_multiplicity > 0:
                for n, neighbor in enumerate(segments[s].links[0]) :
                    supported_links[2*names[segments[s].names[0]], 2*names[neighbor.names[0]]+segments[s].otherEndOfLinks[0][n]] = computed_multiplicity[names[neighbor.names[0]]]
                    supported_links[2*names[neighbor.names[0]]+segments[s].otherEndOfLinks[0][n], 2*names[segments[s].names[0]]] = computed_multiplicity[names[neighbor.names[0]]]
            if new_multiplicity2 == new_multiplicity and new_multiplicity > 0:
                for n, neighbor in enumerate(segments[s].links[1]) :
                    supported_links[2*names[segments[s].names[0]]+1, 2*names[neighbor.names[0]]+segments[s].otherEndOfLinks[1][n]] = computed_multiplicity[names[neighbor.names[0]]]
                    supported_links[2*names[neighbor.names[0]]+segments[s].otherEndOfLinks[1][n], 2*names[segments[s].names[0]]+1] = computed_multiplicity[names[neighbor.names[0]]]
            
        else :
            #check if we can infer the multiplicity of a neighborign contig (possible if there is strictly 1 that is unknown)
            for end in range(2) :
                
                if all([len(i.links[segments[s].otherEndOfLinks[end][j]]) == 1 for j,i in enumerate(segments[s].links[end])])  :
                    number0 = 0
                    index0 = 0
                    for n, neighbor in enumerate(segments[s].links[end]) :
                        if computed_multiplicity[names[neighbor.names[0]]] == 0 :
                            number0 += 1
                            index0 = n
                    if number0 == 1 :
                        new_multiplicity = computed_multiplicity[s] - sum([computed_multiplicity[names[i.names[0]]] for i in segments[s].links[end]])
                         
                        if new_multiplicity > 0 and (segments[s].depths[0]/refCoverage >= new_multiplicity/1.5 or refCoverage == 1):
                            computed_multiplicity[names[segments[s].links[end][index0].names[0]]] = new_multiplicity
                            #print("salut, ", s+1, " ", names[segments[s].links[end][index0].names[0]]+1, " ", new_multiplicity)
                            unchanged = -1
                            
                            supported_links[2*names[segments[s].names[0]]+end, 2*names[segments[s].links[end][index0].names[0]]+segments[s].otherEndOfLinks[end][index0]] = new_multiplicity
                            supported_links[2*names[segments[s].links[end][index0].names[0]]+segments[s].otherEndOfLinks[end][index0], 2*names[segments[s].names[0]]+end] = new_multiplicity
                
        
        s = (s+1)%len(segments)
        unchanged += 1
            
    # print("Propagation finished")    
    
    #now the propagation has stopped, try to see if some multiplicities can be inferred from contig with known multiplicty
    found = False
    if refCoverage != 1 : #that would mean the coverage are not reliable
        for s, seg in enumerate(segments) :
            if computed_multiplicity[s] > 0 :
                for end in range(2) :
                    
                    if all([len(i.links[seg.otherEndOfLinks[end][j]]) ==1 for j,i in enumerate(seg.links[end])]) :
                        
                        covTot = sum([i.depths[0] for i in seg.links[end]])
                        tot = 0
                        
                        if covTot != 0 :
                            for n, neighbor in enumerate(seg.links[end]) :
                                if computed_multiplicity[names[neighbor.names[0]]] == 0 :
                                    computed_multiplicity[names[neighbor.names[0]]] = round(computed_multiplicity[s] * neighbor.depths[0]/covTot)
                                    if round(computed_multiplicity[s] * neighbor.depths[0]/covTot) > 0:
                                        found = True
                                        supported_links[2*names[seg.names[0]]+end, 2*names[neighbor.names[0]]+seg.otherEndOfLinks[end][n]] = round(computed_multiplicity[s] * neighbor.depths[0]/covTot)
                                        supported_links[2*names[neighbor.names[0]]+seg.otherEndOfLinks[end][n], 2*names[seg.names[0]]+end] = round(computed_multiplicity[s] * neighbor.depths[0]/covTot)
                                        #print("Inferring multiplicity of ", neighbor.names, " at ", round(computed_multiplicity[s] * neighbor.depths[0]/covTot) , ", from contig ", seg.names)
                                    propagate_multiplicity(computed_multiplicity, segments, names, names[neighbor.names[0]], supported_links, refCoverage)
        
    #if nothing else worked, start the propagation again by giving a multiplicity to the largest contig without multiplicity yet
    sortedContigs = [names[i.names[0]] for i in sorted(segments, key=lambda x:x.length, reverse=True)]
    sortedContigidx = len(sortedContigs)-1
    while sortedContigidx >= 0 :

        if computed_multiplicity[sortedContigidx] == 0 :
            
    
            minLeft = 0
            minRight = 0
            for n, neighbor in enumerate(segments[sortedContigidx].links[0]) :
                if len(neighbor.links[segments[sortedContigidx].otherEndOfLinks[0][n]]) == 1 :
                    minLeft += computed_multiplicity[names[neighbor.names[0]]]
            for n, neighbor in enumerate(segments[sortedContigidx].links[1]) :
                if len(neighbor.links[segments[sortedContigidx].otherEndOfLinks[1][n]]) == 1 :
                    minRight += computed_multiplicity[names[neighbor.names[0]]]
                
            computed_multiplicity[sortedContigidx] = max(1, minLeft, minRight) #choosing 1 there means that computed_multiplicity is a minimum multiplicity
            propagate_multiplicity(computed_multiplicity, segments, names, sortedContigidx, supported_links, refCoverage)
            
            #print("Finding new multiplicities, ", computed_multiplicity.count(0), " multiplicities left to infer")
            
        sortedContigidx -= 1
    

    return refCoverage, computed_multiplicity

#recursive function to propagate multiplicity :
# input : a contig with a fresh multiplicity
# output : list of multiplicities updated with every multiplicity that could be deduced from this
def propagate_multiplicity(multiplicities, segments, names, contigIdx, supported_links, refCoverage) :
    
    #print("called on ", segments[contigIdx].names[0], " ", multiplicities[contigIdx])
    
    for end in range(2):
        for n, neighbor in enumerate(segments[contigIdx].links[end]):
            
            neighborEnd = segments[contigIdx].otherEndOfLinks[end][n]
            if multiplicities[names[neighbor.names[0]]] == 0 :
                
                if all([multiplicities[names[i.names[0]]] > 0 for i in neighbor.links[neighborEnd]]) and all([len(neighbor.links[neighborEnd][i].links[neighbor.otherEndOfLinks[neighborEnd][i]]) == 1 for i in range(len(neighbor.links[neighborEnd]))]) :
                    
                    multiplicities[names[neighbor.names[0]]] = sum([multiplicities[names[i.names[0]]] for i in neighbor.links[neighborEnd]])
                    
                    for n2, neighbor2 in enumerate(neighbor.links[neighborEnd]) :
                        supported_links[2*names[neighbor.names[0]]+neighborEnd, 2*names[neighbor2.names[0]]+neighbor.otherEndOfLinks[neighborEnd][n2]] = multiplicities[names[neighbor2.names[0]]]
                        supported_links[2*names[neighbor2.names[0]]+neighbor.otherEndOfLinks[neighborEnd][n2], 2*names[neighbor.names[0]]+neighborEnd] = multiplicities[names[neighbor2.names[0]]]
                      
                    propagate_multiplicity(multiplicities, segments, names, names[neighbor.names[0]], supported_links, refCoverage) #recursive call : found a multiplicity, move on
                    
            elif multiplicities[names[neighbor.names[0]]] > 0 :
                
                if all([len(i.links[neighbor.otherEndOfLinks[neighborEnd][j]]) == 1 for j,i in enumerate(neighbor.links[neighborEnd])])  :
                
                    number0 = 0
                    index0 = 0
                    for n2, neighbor2 in enumerate(neighbor.links[neighborEnd]) :
                        if multiplicities[names[neighbor2.names[0]]] == 0 :
                            number0 += 1
                            index0 = n2
                    if number0 == 1 :
                        new_multiplicity = multiplicities[names[neighbor.names[0]]] - sum([multiplicities[names[i.names[0]]] for i in neighbor.links[neighborEnd]])
                         
                        if new_multiplicity > 0 and (neighbor.links[neighborEnd][index0].depths[0]/refCoverage >= new_multiplicity/1.5  or refCoverage==1):
                            multiplicities[names[neighbor.links[neighborEnd][index0].names[0]]] = new_multiplicity
                            # if '75' in neighbor.links[neighborEnd][index0].names :
                            #     print("propagated multiplicity2: found the multiplicity of contig ", neighbor.links[neighborEnd][index0].names[0], " : ", multiplicities[names[neighbor.links[neighborEnd][index0].names[0]]], ", from multiplicity of ", neighbor.names)
                            
                            supported_links[2*names[neighbor.names[0]]+neighborEnd, 2*names[neighbor.links[neighborEnd][index0].names[0]]+neighbor.otherEndOfLinks[neighborEnd][index0]] = new_multiplicity
                            supported_links[2*names[neighbor.links[neighborEnd][index0].names[0]]+neighbor.otherEndOfLinks[neighborEnd][index0], 2*names[neighbor.names[0]]+neighborEnd] = new_multiplicity
                            
                            propagate_multiplicity(multiplicities, segments, names, names[neighbor.links[neighborEnd][index0].names[0]], supported_links, refCoverage)

                
            else : #that's if multiplicity <0, weird
                print("debug WARNING: negative multiplicity")
    
    #now spread the multiplicity based on coverage, if applicable
    if refCoverage != 1 :
        seg = segments[contigIdx]
        for end in range(2) :
                            
            if all([len(i.links[seg.otherEndOfLinks[end][j]]) ==1 for j,i in enumerate(seg.links[end])]) :
                
                covTot = sum([i.depths[0] for i in seg.links[end]])
                tot = 0
                
                for n, neighbor in enumerate(seg.links[end]) :
                    if multiplicities[names[neighbor.names[0]]] == 0 :
                        new_multiplicity = max(min(round(multiplicities[contigIdx] * neighbor.depths[0]/covTot), multiplicities[contigIdx]-len(seg.links[end])+1),1)
                        multiplicities[names[neighbor.names[0]]] = new_multiplicity
                        # if '75' in neighbor.names :
                        #     print("Inferffing the multiplicity of contig ", neighbor.names, " from " , seg.names, " of multiplicity ", multiplicities[contigIdx], " the depth being ", neighbor.depths[0], " and covTot ",  covTot)
                        if round(multiplicities[contigIdx] * neighbor.depths[0]/covTot) > 0:
                            found = True
                            supported_links[2*names[seg.names[0]]+end, 2*names[neighbor.names[0]]+seg.otherEndOfLinks[end][n]] = new_multiplicity
                            supported_links[2*names[neighbor.names[0]]+seg.otherEndOfLinks[end][n], 2*names[seg.names[0]]+end] = new_multiplicity
                            #print("Inferring multiplicity of ", neighbor.names, " at ", round(computed_multiplicity[s] * neighbor.depths[0]/covTot) , ", from contig ", seg.names)
                        propagate_multiplicity(multiplicities, segments, names, names[neighbor.names[0]], supported_links, refCoverage)
        
        
        
        
        
        
        
        
        
        
        
        
        
        