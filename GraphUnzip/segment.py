#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import random

#a segment is a supercontig
class Segment:

    def __init__(self, segNamesOfContig, segOrientationOfContigs, segLengths, segInsideCIGARs = None, segLinks = [[],[]], segOtherEndOfLinks = [[],[]], segCIGARs = [[],[]], lock = False, HiCcoverage = 0, readCoverage = []):
        
        if len(segLinks[0]) != len(segOtherEndOfLinks[0]) or len(segLinks[1]) != len(segOtherEndOfLinks[1]) :
            print('ERROR in the links while initializing a segment')
            return 0
        
        if any(i!=0 and i!=1 for i in segOtherEndOfLinks[0]) or any(i!=0 and i!=1 for i in segOtherEndOfLinks[1]):
            print('ERROR in the links while initializing a segment')
            return 0
        
        if len(segNamesOfContig) != len(segOrientationOfContigs) :
            print('ERROR in initializing the orientations of contigs within a segment')
            return 0
        
        if segInsideCIGARs == None :
            segInsideCIGARs = ['*' for i in range(len(segNamesOfContig)-1)]
            
        if segCIGARs == [[],[]] :
            segCIGARs = [['*' for i in segLinks[0]],['*' for i in segLinks[1]]]
        
        self._id = random.random() #a random number identifying the segment
        self._HiCcoverage = HiCcoverage #to keep in mind how many HiC contacts this segment has in total
        
        #this group of attributes are linked arrays : element n in one corresponds with element n in the other. Therefore they shouldn't be modified independantly
        self._namesOfContigs = segNamesOfContig.copy() #names are strings with which sequences are described in the GFA
        self._orientationOfContigs = segOrientationOfContigs.copy() #1 being '+' orientation, 0 the '-' orientation
        self._lengths = segLengths.copy()
        self._insideCIGARs = segInsideCIGARs.copy()
        if readCoverage != [] :
            self._depths = readCoverage.copy() #to keep in mind the read coverage of the contigs
        else :
            self._depths = [1 for i in range(len(self._lengths))]
            
        self._copiesOfContigs = [-1]*len(segNamesOfContig) #this is used exclusively while exporting, to indicate which copy of which contig is in the segment (copy 0/ copy 1 / copy 2 ...)
        
        #this group of attribute are linked arrays : one should never be modified without the others 
        #They are sorted by ID of the neighbor : important for handling quickly big nodes
        lists_keyed = [[(segLinks[0][i], segOtherEndOfLinks[0][i], segCIGARs[0][i]) for i in range(len(segLinks[0]))], [(segLinks[1][i], segOtherEndOfLinks[1][i], segCIGARs[1][i]) for i in range(len(segLinks[1]))]]
        lists_keyed[0].sort(key = lambda x: x[0].ID)
        lists_keyed[1].sort(key = lambda x: x[0].ID)
        self._links = [[i[0] for i in lists_keyed[0]], [i[0] for i in lists_keyed[1]]] #two lists of segments with which the segment is linked, at the left end and at the right end
        self._otherEndOfLinks = [[i[1] for i in lists_keyed[0]], [i[1] for i in lists_keyed[1]]] #for each link, indicates the side of the other segment on which the link arrives
        self._CIGARs = [[i[2] for i in lists_keyed[0]], [i[2] for i in lists_keyed[1]]] #for each link, indicates the CIGAR string found in the GFA
        
        self._freezed = [False, False] #do not duplicate from one end if frozen at this end
        self._locked = lock #That is to duplicate a contig only once in each merge_contigs
        
        
    # getters
    
    def get_id(self):
        return self._id
    
    # def __hash__(self):
    #     return self._id
    
    def get_orientations(self):
        return self._orientationOfContigs
    
    def get_insideCIGARs(self):
        return self._insideCIGARs
    
    def get_links(self):
        return self._links
    
    def get_otherEndOfLinks(self):
        return self._otherEndOfLinks
    
    def get_CIGARs(self):
        return self._CIGARs
    
    def get_lengths(self):
        return self._lengths
    
    def get_length(self):
        return np.sum(self._lengths)
    
    def get_namesOfContigs(self):
        return self._namesOfContigs
    
    def get_copiesOfContigs(self):
        return self._copiesOfContigs
    
    def get_freezed(self):
        return self._freezed
    
    def get_locked(self):
        return self._locked
    
    def get_coverage(self):
        return self._HiCcoverage
        
    def get_depths(self):
        return self._depths
    
    def get_depth(self):
        sumdepth = 0
        sumlength = 1 #1 and not 0 to be sure not to divide by 0
        for i in range(len(self._depths)) :
            sumdepth += self._depths[i]*self._lengths[i]
            sumlength += self._lengths[i]
            
        return sumdepth / sumlength
    
    def full_name(self) :
        return '_'.join([self._namesOfContigs[i]+'-'+str(self._copiesOfContigs[i]) for i in range(len(self._namesOfContigs))])
    
    def print_complete(self):
        print(self._namesOfContigs, [s.names for s in self._links[0]], \
              [s.names for s in self._links[1]])
    
    # setters 
    def set_copiesNumber(self, copiesNumberForNow):
        for c, contig in enumerate(self._namesOfContigs) :
            if contig in copiesNumberForNow :
                self._copiesOfContigs[c] = copiesNumberForNow[contig]
                copiesNumberForNow[contig] += 1
            else :
                self._copiesOfContigs[c] = 0
                copiesNumberForNow[contig] = 1
    
    def set_coverage(self, newCoverage) :
        self._HiCcoverage = newCoverage
    
    def freeze(self, endOfSegment): 
        self._freezed[endOfSegment] = True
  
    def freezeNode(self, endOfSegment):
        self._freezed[endOfSegment] = True
        for n, neighbor in enumerate(self._links[endOfSegment]) :
            neighbor.freeze(self._otherEndOfLinks[endOfSegment][n])
        
    def unfreeze(self):
        self._freezed = [False, False]
        
    def set_locked(self, b):
        self._locked = b
        
    def lockNode(self, endOfSegment):
        self._locked = True
        for i in self._links[endOfSegment]:
            i.locked = True
            
    def set_id(self, newID) : #few cases where that is useful
        self._id = newID 
        
    def divide_depths(self, n) : #when duplicating a segment, you need to lower the coverage of all replicas
        for i in range (len(self._depths)) :
            self._depths[i] /= n
            
    def multiply_end_depths(self, n, end, numberOfContigs) : #this fucntion is useful when merging a wrongly duplicated dead end
        for co in range(numberOfContigs) :
           self._depths[end*(len(self._depths)-1) + (-2*end+1)*co ]
            
    def length1(self): #use this function to set length of a segment to 1 (instead of 0, mostly)
        self._lengths = [max(1, i) for i in self._lengths]
                    
    # properties
    
    ID = property(get_id, set_id)
    HiCcoverage = property(get_coverage, set_coverage)
    depths = property(get_depths)
    depth = property(get_depth)
    length = property(get_length)
    
    names = property(get_namesOfContigs)
    orientations = property(get_orientations)
    lengths = property(get_lengths)
    insideCIGARs = property(get_insideCIGARs)
    
    copiesnumber = property(get_copiesOfContigs)
    
    links = property(get_links)
    otherEndOfLinks = property(get_otherEndOfLinks)
    CIGARs = property(get_CIGARs)
    
    freezed = property(get_freezed) #segment is freezed if comparison of links was inconclusive
    locked = property(get_locked, set_locked) #segment is locked if it was already duplicated in merge_contigs and that it should not be for a second time

    #other functions that handle segments
    
    #function which goals is to return the intensity of Hi-C contacts between another segment and this one
    def interaction_with_contigs(self, segment, interactionMatrix, names, copiesnumber = None, commonContigs = set(), bestSignature = 1000):
        
        if copiesnumber == None :
            copiesnumber = [1 for i in interactionMatrix]
            
        absoluteScore = 0
        relativeScore = 0
        
        # orientation = -1 # if supercontig is directly linked to the candidates, then this variable tells us by which end
        
        # if segment in self.links[0]:
        #     orientation = 0
        # elif segment in self.links[1]:
        #     orientation = 1
        # else: 
        #     print('ERROR : trying to compute an interaction but the contigs do not touch each other')
        #     print('Looking for ', self._namesOfContigs, ' from ', segment.names)
        #     return 0, 0, 1
            
        depth = 1
        #first compute interactions with self
        for co, contig in enumerate(self.names) :
                 
            for c, contigInSegment in enumerate(segment.names):
                
                if contig not in commonContigs and copiesnumber[contigInSegment] <= bestSignature:
                    
                    absoluteScore += interactionMatrix[names[contig], names[contigInSegment]]
                    relativeScore += interactionMatrix[names[contig], names[contigInSegment]]
                else:
                    absoluteScore += interactionMatrix[names[contig], names[contigInSegment]]
                                        
        # #now compute the interaction with neighbors of self 
        # endOfSegment = 1-orientation
        # for neighbor in self.links[endOfSegment] :
        #     for co, contig in enumerate(neighbor.names) :
        #         for c, contigInSegment in enumerate(segment.names):
                
        #             if contig not in commonContigs and copiesnumber[contigInSegment] <= bestSignature:
                        
        #                 depth = 2
                        
        #                 absoluteScore += interactionMatrix[names[contigInSegment],names[contig]]
        #                 relativeScore += interactionMatrix[names[contigInSegment],names[contig]]
        #             else:
        #                 absoluteScore += interactionMatrix[names[contigInSegment],names[contig]]
                
            
        return absoluteScore, relativeScore, depth
    
    def add_link_from_GFA(self, GFAline, names, segments, leftOrRight) : #leftOrRight = 0 when the segment is at the beginning of a link (left of a GFA line), 1 otherwise
        
        l = GFAline.strip('\n').split('\t')
 
        if len(l) < 5 :
            print('ERROR : expected at least 5 fields in line ', GFAline)
        
        if l[0] != 'L':
            print('ERROR : trying to add a link from a GFA line that does not start with "L"')
        
        else :
            o1,o2 = -1, -1
            
            if l[2] == '-':
                o1 = 0
            elif l[2] == '+' :
                o1 = 1
                
            if l[4] == '-':
                o2 = 0
            elif l[4] == '+' :
                o2 = 1
                
            if o1 == -1 or o2 == -1 :
                print('ERROR while creating a link : orientations not properly given.')
                print('Problematic line : ', GFAline)   

            
            if leftOrRight == 0 and o1 == 0:

                #then comes a little test to see if the link has already been added (for example if their are several lines in the GFA describing the same link). An 'or' condition to allow it when an end of link is linked with itself
                if find_this_link(segments[names[l[3]]], 1-o2, self._links[0], self._otherEndOfLinks[0]) == -1 or (segments[names[l[3]]].ID == self._id and 1-o2 == 0):

                    index = index_at_which_new_link_should_be_inserted(segments[names[l[3]]], self._links[0], 1-o2 ,self._otherEndOfLinks[0])
                    self._links[0].insert(index, segments[names[l[3]]])
                    self._otherEndOfLinks[0].insert(index, 1-o2)
                    if len(l) > 5 :
                        self._CIGARs[0].insert(index, l[5])
                    else :
                        self._CIGARs[0].insert(index, '*')
                    
            elif leftOrRight == 0 and o1 == 1 :
                if find_this_link(segments[names[l[3]]], 1-o2, self._links[1], self._otherEndOfLinks[1]) == -1 or (segments[names[l[3]]].ID == self._id and 1-o2 == 1):
                    
                    index = index_at_which_new_link_should_be_inserted(segments[names[l[3]]], self._links[1],  1-o2 ,self._otherEndOfLinks[1])
                    
                    self._links[1].insert(index, segments[names[l[3]]])
                    self._otherEndOfLinks[1].insert(index, 1-o2)
                    if len(l) > 5 :
                        self._CIGARs[1].insert(index, l[5])
                    else :
                        self._CIGARs[1].insert(index, '*')
                
            elif leftOrRight == 1 and o2 == 1 :
                if find_this_link(segments[names[l[1]]], o1, self._links[0], self._otherEndOfLinks[0]) == -1 or (segments[names[l[1]]].ID == self._id and o1 == 0) :
                    
                    index = index_at_which_new_link_should_be_inserted(segments[names[l[1]]], self._links[0],  o1 ,self._otherEndOfLinks[0])
                    self._links[0].insert(index, segments[names[l[1]]])
                    self._otherEndOfLinks[0].insert(index, o1)
                    if len(l) > 5 :
                        self._CIGARs[0].insert(index, l[5])
                    else :
                        self._CIGARs[0].insert(index, '*')
                    
            elif leftOrRight == 1 and o2 == 0 :
                if find_this_link(segments[names[l[1]]], o1, self._links[1], self._otherEndOfLinks[1]) == -1 or (segments[names[l[1]]].ID == self._id and o1 == 1):
                    index = index_at_which_new_link_should_be_inserted(segments[names[l[1]]], self._links[1],  o1 ,self._otherEndOfLinks[1])
                    self._links[1].insert(index, segments[names[l[1]]])
                    self._otherEndOfLinks[1].insert(index, o1)
                    if len(l) > 5 :
                        self._CIGARs[1].insert(index, l[5])
                    else :
                        self._CIGARs[1].insert(index, '*')
            
            else :
                print('ERROR while trying to add a new link from the gfa : could not locate a correct name')
        
    #this adds the end of a links, but only on this segment, not on the other end
    def add_end_of_link(self, endOfSegment, segment2, endOfSegment2, CIGAR = '*'):
        
        #print('A', len(segment2.otherEndOfLinks[1]), len(segment2.links[1]), len(segment2.CIGARs[1]))
        #print(self._namesOfContigs, segment2.names)
        index = index_at_which_new_link_should_be_inserted(segment2, self._links[endOfSegment], endOfSegment2, self._otherEndOfLinks[endOfSegment])

        self._links[endOfSegment].insert(index, segment2)
        #print('B', len(segment2.otherEndOfLinks[1]), len(segment2.links[1]), len(segment2.CIGARs[1]))

        self._otherEndOfLinks[endOfSegment].insert(index, endOfSegment2)
        self._CIGARs[endOfSegment].insert(index, CIGAR)        

    #this function is useful for rerouting around big nodes. It adds end of links more efficiently than if done individually, and checks for doubles
    def add_a_bunch_of_end_of_links(self, endOfSegment, listOfSegmentsToAdd, listOfEndOfSegmentsToAdd, CIGARsToAdd) :
 
        #All the list of segments being sorted by ID, this is done on the basis of the merging of merge sort.
        newLinks = []
        newEndOfLinks = []
        newCIGARs = []

        indexOfSegmentToAdd = 0
        indexOfLinks = 0
        
        while indexOfSegmentToAdd < len(listOfSegmentsToAdd) and indexOfLinks < len(self._links[endOfSegment]) :
            
            if listOfSegmentsToAdd[indexOfSegmentToAdd].ID < self._links[endOfSegment][indexOfLinks].ID :
                newLinks.append(listOfSegmentsToAdd[indexOfSegmentToAdd])
                newEndOfLinks.append(listOfEndOfSegmentsToAdd[indexOfSegmentToAdd])
                newCIGARs.append( CIGARsToAdd[indexOfSegmentToAdd] )
                indexOfSegmentToAdd += 1
                
            elif listOfSegmentsToAdd[indexOfSegmentToAdd].ID > self._links[endOfSegment][indexOfLinks].ID :
                newLinks.append( self._links[endOfSegment][indexOfLinks] )
                newEndOfLinks.append( self._otherEndOfLinks[endOfSegment][indexOfLinks] )
                newCIGARs.append( self._CIGARs[endOfSegment][indexOfLinks] )
                indexOfLinks += 1
                
            #elsewhise the link goes toward the same segment, but maybe not the same end of this segment
            elif listOfEndOfSegmentsToAdd[indexOfSegmentToAdd] < self._otherEndOfLinks[endOfSegment][indexOfLinks] :
                newLinks.append( listOfSegmentsToAdd[indexOfSegmentToAdd] )
                newEndOfLinks.append( listOfEndOfSegmentsToAdd[indexOfSegmentToAdd] )
                newCIGARs.append( CIGARsToAdd[indexOfSegmentToAdd] )
                indexOfSegmentToAdd += 1
                
            elif listOfEndOfSegmentsToAdd[indexOfSegmentToAdd] > self._otherEndOfLinks[endOfSegment][indexOfLinks] :
                newLinks.append( self._links[endOfSegment][indexOfLinks] )
                newEndOfLinks.append( self._otherEndOfLinks[endOfSegment][indexOfLinks] )
                newCIGARs.append( self._CIGARs[endOfSegment][indexOfLinks] )
                indexOfLinks += 1
                
            else : #it means that there is twice the same link, so add it only once
                newLinks.append( self._links[endOfSegment][indexOfLinks] )
                newEndOfLinks.append( self._otherEndOfLinks[endOfSegment][indexOfLinks])
                newCIGARs.append( self._CIGARs[endOfSegment][indexOfLinks])
                indexOfLinks += 1
                indexOfSegmentToAdd += 1
                
        if indexOfSegmentToAdd >= len(listOfEndOfSegmentsToAdd) :
            newLinks += self._links[endOfSegment][indexOfLinks:]
            newEndOfLinks += self._otherEndOfLinks[endOfSegment][indexOfLinks:]
            newCIGARs += self._CIGARs[endOfSegment][indexOfLinks:]
            
        elif indexOfLinks >= len(self._links[endOfSegment]) :
            newLinks += listOfSegmentsToAdd[indexOfSegmentToAdd:]
            newEndOfLinks += listOfEndOfSegmentsToAdd[indexOfSegmentToAdd:]
            newCIGARs += CIGARsToAdd[indexOfSegmentToAdd:]
        
        self._links[endOfSegment] = newLinks
        self._otherEndOfLinks[endOfSegment] = newEndOfLinks
        self._CIGARs[endOfSegment] = newCIGARs

    def remove_end_of_link(self, endOfSegment, segmentToRemove, endOfSegmentToRemove = None, warning = True): #endOfSegmentToRemove is there in case there exists two links between self[endOfSegment] and segment to remove. Needed for extra security
        
        #first determine the index of the segment to remove
        #print('Removing ', segmentToRemove.names, endOfSegmentToRemove, ' from ', self._namesOfContigs)
        #print('Among these links :', [i.names for i in self._links[endOfSegment]], self._otherEndOfLinks[endOfSegment])
        index = find_this_link(segmentToRemove, endOfSegmentToRemove, self._links[endOfSegment], self._otherEndOfLinks[endOfSegment], warning = warning)
        #index = self._links[endOfSegment].index(segmentToRemove)
   
        #then remove the end of unwanted link in all attributes
        if index != -1 :
            del self._links[endOfSegment][index]
            del self._otherEndOfLinks[endOfSegment][index]
            del self._CIGARs[endOfSegment][index]
            return True
        elif index == -1 and warning:
             print('Trying unsuccesfully to remove ', segmentToRemove.names, ' from ', self._namesOfContigs)
             return False
        
    #returns two contigs, equal to this contig but split at axis, corresponding to the number of contigs left of the junction
    def break_contig(self, axis) :
        
        newSegment1 = Segment(self._namesOfContigs[:axis], self._orientationOfContigs[:axis], self._lengths[:axis], self._insideCIGARs[:axis-1], [self._links[0], []], [self._otherEndOfLinks[0], []], [self._CIGARs[0], []], readCoverage = self._depths[:axis])
        
        newSegment2 = Segment(self._namesOfContigs[axis:], self._orientationOfContigs[axis:], self._lengths[axis:], self._insideCIGARs[axis:], [[], self._links[1]], [[], self._otherEndOfLinks[1]], [[], self._CIGARs[1]], readCoverage = self._depths[axis:])
        
        return newSegment1, newSegment2
    
     
    #function to be used on small loops only
    def flatten(self, replicas) :
        if self not in self._links[0] :
            
            print('ERROR : in segment.flatten, trying to flatten something that is not a loop')
            
        else :
            
            for i in range(len(self._depths)) :
                self._depths[i] /= replicas+1
            
            newName = self._namesOfContigs.copy()
            newOrientations = self._orientationOfContigs.copy()
            newLengths = self._lengths.copy()
            newinsideCIGARs = self._insideCIGARs.copy()
            newCopies = self._copiesOfContigs.copy()
            newDepths = self._depths.copy()
            for i in range(replicas) :
                newName += self._namesOfContigs
                newOrientations += self._orientationOfContigs
                newLengths += self._lengths
                newCopies += self._copiesOfContigs
                newinsideCIGARs += [self._CIGARs[0][self._links[0].index(self)]] + self._insideCIGARs
                newDepths += self._depths
            
            self._namesOfContigs = newName
            self._orientationOfContigs = newOrientations
            self._lengths = newLengths
            self._copiesOfContigs = newCopies
            self._insideCIGARs = newinsideCIGARs
            self._depths = newDepths
            
           # print('In segment.flatten : ', self._namesOfContigs, self._insideCIGARs, [self._CIGARs[0][self._links[0].index(self)]])
            #print('Links before any removal, ', [i.names for i in self._links[0]], '\n')
            self.remove_end_of_link(0, self)
            self.remove_end_of_link(1, self)
  
    #a function to delete all the links of a segment (typically before deleting it)
    def cut_all_links(self) :
        
        for end in range(2) :
            
            for n, neighbor in enumerate(self._links[end]) :
                
                otherEnd = self._otherEndOfLinks[end][n]
                neighbor.remove_end_of_link(otherEnd, self, end)
        
        self._links = [[],[]]
        self._otherEndOfLinks = [[],[]]
        self._CIGARs = [[],[]]
        
#This function is OUTSIDE the class. It takes two segments and the end of the first segment which is linked to the second. It appends a merged contig to the listOfSegments, without modifying the two inputed segments
def merge_two_segments(segment1, endOfSegment1, segment2, listOfSegments):
    
    if segment1.links[endOfSegment1].count(segment2) > 1 : #this means a loop
        return 0
      
    # creating a new segment
    orientation1 = endOfSegment1*2-1
    endOfSegment2 = segment1.otherEndOfLinks[endOfSegment1][segment1.links[endOfSegment1].index(segment2)]
    orientation2 = -2*endOfSegment2+1
    
    orientationOfContigs1 = segment1.orientations
    orientationOfContigs2 = segment2.orientations
    
    if orientation1 == -1 : #then change the orientation of all the contigs within the segment, since the segment will be mirrored in the new supersegment
        orientationOfContigs1 = [1-i for i in orientationOfContigs1]
    if orientation2 == -1 :
        orientationOfContigs2 = [1-i for i in orientationOfContigs2]
        
    CIGAR = segment1.CIGARs[endOfSegment1][segment1.links[endOfSegment1].index(segment2)]
    
    newSegment = Segment(segment1.names[::orientation1] + segment2.names[::orientation2],\
                                orientationOfContigs1[::orientation1]+orientationOfContigs2[::orientation2],\
                                segment1.lengths[::orientation1]+segment2.lengths[::orientation2], \
                                segment1.insideCIGARs[::orientation1] + [CIGAR] + segment2.insideCIGARs[::orientation2],\
                                segLinks = [segment1.links[1-endOfSegment1], \
                                segment2.links[1-endOfSegment2]], \
                                segOtherEndOfLinks = [segment1.otherEndOfLinks[1-endOfSegment1], \
                                segment2.otherEndOfLinks[1-endOfSegment2]],\
                                segCIGARs = [segment1.CIGARs[1-endOfSegment1], \
                                segment2.CIGARs[1-endOfSegment2]],
                                lock = True,
                                HiCcoverage = segment1.HiCcoverage + segment2.HiCcoverage,
                                readCoverage = segment1.depths[::orientation1] + segment2.depths[::orientation2])
            
    listOfSegments.append(newSegment)
    
    #building the other end of links with the new segment
    self_loop_CIGAR = ''
    for n, neighbor in enumerate(newSegment.links[0]) :
        neighbor.add_end_of_link(newSegment.otherEndOfLinks[0][n], newSegment, 0, CIGAR = newSegment.CIGARs[0][n])
        if neighbor.ID == segment2.ID and newSegment.otherEndOfLinks[0][n] == 1 - endOfSegment2 : #check if the new contig should loop back on itself
            self_loop_CIGAR = newSegment.CIGARs[0][n]

    for n, neighbor in enumerate(newSegment.links[1]) :
        #print(len(newSegment.otherEndOfLinks[1]), len(newSegment.links[1]), len(newSegment.CIGARs[1]))
        neighbor.add_end_of_link(newSegment.otherEndOfLinks[1][n], newSegment, 1, CIGAR = newSegment.CIGARs[1][n])
        
    if self_loop_CIGAR != '' :
        add_link(newSegment, 0, newSegment, 1, self_loop_CIGAR)

#function creating a link between two ends of contigs, OUTSIDE of the class
def add_link(segment1, end1, segment2, end2, CIGAR = '*'):
    segment1.add_end_of_link(end1, segment2, end2, CIGAR)
    segment2.add_end_of_link(end2, segment1, end1, CIGAR)
    
def delete_link(segment1, end1, segment2, end2, warning = True) :
    success1 = segment1.remove_end_of_link(end1, segment2, end2, warning = warning)
    success2 = segment2.remove_end_of_link(end2, segment1, end1, warning = warning)
    return success1 and success2
           
def compute_copiesNumber(listOfSegments):
    cn = {}
    for s in listOfSegments :
        s.set_copiesNumber(cn)
        
    return cn

#returns the position of the link pointing towards segment and its endOfSegment
def find_this_link(segment, endOfSegment, listOfLinks, listOfEndsOfLinks, warning = False) :

    lo = 0
    hi = len(listOfLinks)

    while lo < hi:
        mid = (lo+hi)//2
        if segment.ID < listOfLinks[mid].ID:
            hi = mid
        elif segment.ID > listOfLinks[mid].ID:
            lo = mid+1
        else :
            #print('Found : ', endOfSegment , listOfEndsOfLinks[mid])
            if endOfSegment == None :
                return mid
                
            elif endOfSegment == listOfEndsOfLinks[mid] :
                return mid
                
            elif endOfSegment > listOfEndsOfLinks[mid] :
                mid += 1
                while mid < len(listOfLinks) and listOfLinks[mid].ID == segment.ID :
                    if endOfSegment == listOfEndsOfLinks[mid] :
                        return mid
                    mid += 1
                    
                break
                
            elif endOfSegment < listOfEndsOfLinks[mid] :
                mid -= 1
                while mid >= 0 and listOfLinks[mid].ID == segment.ID :
                    if endOfSegment == listOfEndsOfLinks[mid] :
                        return mid
                    mid -= 1
                break
    
    if not warning :
        return -1
        
    print('In find_this_link : did not find the link')
    #print([[listOfLinks[se].names, listOfEndsOfLinks[se]] for se in range(len(listOfLinks))])
    print('Did not find ', segment.names , endOfSegment, ' among ', [i.names for i in listOfLinks], listOfEndsOfLinks)
    
    return -1

#returns the index at which a segment should be inserted in a list sorted by ID : useful because links[0] and links[1] are kept sorted at all times
def index_at_which_new_link_should_be_inserted(segment, listOfSegments, endOfLink, listOfEndOfLinks) :
    lo = 0
    hi = len(listOfSegments)

    while lo < hi:
        mid = (lo+hi)//2
        if segment.ID < listOfSegments[mid].ID or (segment.ID == listOfSegments[mid].ID and endOfLink < listOfEndOfLinks[mid]):
            hi = mid
        else:
            lo = mid+1
            
    while lo < len(listOfSegments) and listOfSegments[lo].ID == segment.ID and endOfLink == listOfEndOfLinks[lo] :
        lo += 1
    return lo

def check_if_all_links_are_sorted(listOfSegments) :
    
    for segment in listOfSegments :
        for endOfSegment in range(2) :
            for n in range(len(segment.links[endOfSegment])-1) :
                if segment.links[endOfSegment][n].ID > segment.links[endOfSegment][n+1].ID :
                    print('Problem in the links of ', segment.names, ' : ', [s.ID for s in segment.links[endOfSegment]])
    
    for segment in listOfSegments :
        for endOfSegment in range(2) :
            for n, neighbor in enumerate(segment.links[endOfSegment]) :
                if segment not in neighbor.links[segment.otherEndOfLinks[endOfSegment][n]] :
                    print('Non-reciprocal links : ', segment.names, segment.ID, neighbor.names, neighbor.ID)
                
#funtion to delete links that are present twice in the graph (often because they are present twice in the gfa)
def delete_links_present_twice(segments):
    

    for segment in segments :
        toBeRemoved = []
        for endOfSegment in range(2) :
            
            for n1 in range(len(segment.links[endOfSegment])-1) :
                
                for n2 in range(n1+1, len(segment.links[endOfSegment])) :
                    
                    if segment.links[endOfSegment][n1].ID == segment.links[endOfSegment][n2].ID and segment.otherEndOfLinks[endOfSegment][n1] == segment.otherEndOfLinks[endOfSegment][n2] and segment.links[endOfSegment][n1].ID != segment.ID:
                        
                        segment.links[endOfSegment][n2].remove_end_of_link(segment.otherEndOfLinks[endOfSegment][n2], segment, endOfSegment)
                        toBeRemoved += [[endOfSegment, segment.links[endOfSegment][n2], segment.otherEndOfLinks[endOfSegment][n2]]]

        for r in toBeRemoved :
            segment.remove_end_of_link(r[0], r[1], r[2])


## A few lines to test the functions of the file

# s1 = Segment([0,1], [1,1], [1000])
# s2 = Segment([2], [1], [500])
# add_link(s1,1,s2,0,'60M')

# listOfSegments = [s1, s2]

# merge_two_segments(s1, 1, s2, listOfSegments)

# print([i.orientations for i in listOfSegments])

                    
                    
                    