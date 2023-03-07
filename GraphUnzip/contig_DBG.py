#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 12:37:18 2021

@author: rfaure
"""

from input_output import read_GAF
import segment as sg
import sys
import re
from copy import deepcopy

#define reverse complement function for the paths
def reverse_complement(path):
    reverse_comp = []
    for i in range (len(path)-1, -1, -1):
        orientation = ">"
        if path[i][1] == '>':
            orientation = '<'
        reverse_comp.append((path[i][0], orientation))
    return reverse_comp

#function to output the graph in gdf format
def output_graph(graph_nodes, graph_edges, filename):
    #output the graph in gdf format
    gdfFile = open(filename, "w")
    gdfFile.write("nodedef>name VARCHAR, label VARCHAR\n")
    for node in graph_nodes:
        name = ""
        for i in range(len(node)):
            # print("node : ", node)
            name += str(node[i][0])
            name += str(node[i][1])
        for end in range(2):
            gdfFile.write(name+"_"+str(end)+"," + name+"_"+str(end) + "\n")
    
    gdfFile.write("edgedef>node1 VARCHAR, node2 VARCHAR\n")
    
    for node in graph_nodes:
        name = ""
        for i in range(len(node)):
            # print("node : ", node)
            name += str(node[i][0])
            name += str(node[i][1])
        gdfFile.write(name + "_0," + name+"_1\n")

    for edge in graph_edges:
        name1 = ""
        name2 = ""
        # print("edge : ", edge)
        for i in range(len(edge[0][0])):
            name1 += str(edge[0][0][i][0]) #first kmer, the name
            name1 += str(edge[0][0][i][1])
        name1 += "_"+str(int(edge[0][1]))
        for i in range(len(edge[1][0])):
            name2 += str(edge[1][0][i][0])
            name2 += str(edge[1][0][i][1])
        name2 += "_"+str(int(edge[1][1]))
        gdfFile.write(name1 + "," + name2+"\n")

#supress all kmers that have an abundance lower than the limit   
def clean_graph(graph_nodes, graph_edges, graph_edges_dict, kmer_abundance, limit_abundance) :
    #delete from the graph all the kmers that have an abundance less than limit_abundance
    for kmer in kmer_abundance:
        if kmer_abundance[kmer] <= limit_abundance:
            graph_nodes.remove(kmer)
            del graph_edges_dict[(kmer, True)]
            del graph_edges_dict[(kmer, False)]
    
    toDelete = set()
    for edge in graph_edges:
        if edge[0][0] not in graph_nodes or edge[1][0] not in graph_nodes:
            toDelete.add(edge)

    for delete in toDelete:
        graph_edges.remove(delete)

    for node in graph_nodes:
        toDelete = set()
        for edge in graph_edges_dict[(node, True)]:
            if edge[0] not in graph_nodes:
                toDelete.add(edge)
        for delete in toDelete:
            graph_edges_dict[(node, True)].remove(delete)

        toDelete = set()
        for edge in graph_edges_dict[(node, False)]:
            if edge[0] not in graph_nodes:
                toDelete.add(edge)
        for delete in toDelete:
            graph_edges_dict[(node, False)].remove(delete)

#function to build a DBG graph from the paths
def build_graph(k, paths, graph_nodes, graph_edges, graph_edges_dict, kmer_abundance):
    #build a de Bruijn graph from the paths for k = 2
    for path in paths:

        previous_kmer = []
        previous_orientation = ''
        for seg_start in range(len(path)-k+1):
            #if the k-mer is not in the graph, add them
            kmer_foward = tuple(path[seg_start:seg_start+k])
            kmer = kmer_foward
            kmer_reverse = tuple(reverse_complement(path[seg_start:seg_start+k]))
            reverse = False
            if hash(kmer_foward) < hash(kmer_reverse):
                kmer = kmer_reverse
                reverse = True

            # print("kmer : ", kmer, " and orientation : ", orientation, " and path[seg_start:seg_start+k] : ", path[seg_start:seg_start+k], " and reverse_complement(path[seg_start:seg_start+k]) : ", reverse_complement(path[seg_start:seg_start+k]))

            if kmer not in graph_nodes:
                # print("kmer : ", kmer)
                if kmer_foward[0] == ("edge_31@3@0",">") or kmer_reverse[0] == ("edge_31@3@0",">")  :
                    print("for k = ", k,  " kmer : ", kmer, " and orientation : ", orientation)
                    print("from path : ", path)
                graph_nodes.add(kmer)
                graph_edges_dict[(kmer, True)] = set()
                graph_edges_dict[(kmer, False)] = set()
                kmer_abundance[kmer] = 0

            kmer_abundance[kmer] += 1

            orientation = False 
            if reverse :
                orientation = True
            
            #add the edge to the right
            if seg_start > 0 :
                # print("adding previous kmer : ", previous_kmer, " and orientation : ", previous_orientation, " to kmer : ", kmer, " and orientation : ", orientation)
                graph_edges.add(((previous_kmer, previous_orientation), (kmer, orientation)))
                graph_edges_dict[(previous_kmer, previous_orientation)].add((kmer, orientation))
                graph_edges_dict[(kmer, orientation)].add((previous_kmer, previous_orientation))

                # if kmer_foward[0] == ("edge_31@3@0",">") or kmer_reverse[0] == ("edge_31@3@0",">")  :
                #     print("for k = ", k,  " kmer : ", kmer, " and orientation : ", orientation)
                #     print("adding an edge from previous kmer : ", previous_kmer, " and orientation : ", previous_orientation)

            previous_kmer = kmer
            previous_orientation = not orientation
    
    return graph_nodes, graph_edges, graph_edges_dict, kmer_abundance

#function to build a graph from the given paths, but not DBG
def list_existing_links(k, paths, graph_nodes, graph_edges, graph_edges_dict, kmer_abundance):
    #build a de Bruijn graph from the paths for k = 2
    for path in paths:

        previous_kmer = []
        previous_orientation = ''
        for seg_start in range(len(path)-k+1):
            #if the k-mer is not in the graph, add them
            kmer_foward = tuple(path[seg_start:seg_start+k])
            kmer = kmer_foward
            kmer_reverse = tuple(reverse_complement(path[seg_start:seg_start+k]))
            reverse = False
            if hash(kmer_foward) < hash(kmer_reverse):
                kmer = kmer_reverse
                reverse = True

            # print("kmer : ", kmer, " and orientation : ", orientation, " and path[seg_start:seg_start+k] : ", path[seg_start:seg_start+k], " and reverse_complement(path[seg_start:seg_start+k]) : ", reverse_complement(path[seg_start:seg_start+k]))

            if kmer not in graph_nodes:
                graph_nodes.add(kmer)
                graph_edges_dict[(kmer, True)] = set()
                graph_edges_dict[(kmer, False)] = set()
                kmer_abundance[kmer] = 0

            kmer_abundance[kmer] += 1

            orientation = False 
            if reverse :
                orientation = True
            
            #add the edge to the right
            if seg_start > 0 :
                # print("adding previous kmer : ", previous_kmer, " and orientation : ", previous_orientation, " to kmer : ", kmer, " and orientation : ", orientation)
                graph_edges.add(((previous_kmer, previous_orientation), (kmer, orientation)))
                graph_edges_dict[(previous_kmer, previous_orientation)].add((kmer, orientation))
                graph_edges_dict[(kmer, orientation)].add((previous_kmer, previous_orientation))

                # if kmer_foward[0] == ("edge_31@3@0",">") or kmer_reverse[0] == ("edge_31@3@0",">")  :
                #     print("for k = ", k,  " kmer : ", kmer, " and orientation : ", orientation)
                #     print("adding an edge from previous kmer : ", previous_kmer, " and orientation : ", previous_orientation)

            previous_kmer = kmer
            previous_orientation = not orientation
    
    return graph_nodes, graph_edges, graph_edges_dict, kmer_abundance

#function that builds contigs from the graph
def build_contigs (graph_nodes, graph_edges, graph_edges_dict, k):
    contigs = []
    already_visited = set()
    for node in graph_nodes:
        #check if it can be the beginning of a contig
        if node not in already_visited:
            starting_node = node
            starting_end = False

            while len(graph_edges_dict[(starting_node, starting_end)]) == 1 and len(graph_edges_dict[list(graph_edges_dict[(starting_node, starting_end)])[0]]) == 1 :
                new_node = list(graph_edges_dict[(starting_node, starting_end)])[0][0]
                starting_end = not list(graph_edges_dict[(starting_node, starting_end)])[0][1]
                starting_node = new_node
                if starting_node == node : #we did a loop !
                    break
                
            #now we have the starting node, we can build the contig starting from it
            contig = []
            starting_end = not starting_end
            firstKmer = deepcopy(starting_node)
            if not starting_end:
                firstKmer = tuple(reverse_complement(list(starting_node)))
            contig = list(firstKmer)
            already_visited.add(starting_node)
            node = starting_node
            end = starting_end
            while len(graph_edges_dict[(node, end)]) == 1 and len(graph_edges_dict[list(graph_edges_dict[(node, end)])[0]]) == 1 :
                new_node = list(graph_edges_dict[(node, end)])[0][0]
                new_end = not list(graph_edges_dict[(node, end)])[0][1]
                if new_node == starting_node and new_end == starting_end:
                    break
                
                newkmer = deepcopy(new_node)
                if not new_end:
                    newkmer = tuple(reverse_complement(list(new_node)))
                
                contig += [newkmer[k-1]]

                already_visited.add(new_node)
                node = new_node
                end = new_end
            
            contigs.append(contig)

    return contigs


#Master function of the file
#Input : initial gfa (as a list of segments), a GAF file with long reads mapped to the segments, names (which is an index numbering the contig)
#Output : new gfa (as a list of segments) corrected with long reads, and modified copiesnumber (taking into account contigs that have been duplicated)
def DBG_long_reads(segments, names, copiesnumber, gafFile):

    lines = []
    print("Reading the gaf file...")
    read_GAF(gafFile, 0, 0, lines)
    print("Finished going through the gaf file.")

    paths = []
    for line in lines:
        #split the line on '>' and '<' to get the path
        contigs = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del contigs[0] #because the first element is always ''

        paths.append([(contigs[i], orientations[i]) for i in range(len(contigs))])

        # if contigs[0] == "edge_322@0@0" :
        #     print("contigs : ", contigs, " and orientations : ", orientations, " and line : ", line)
        #     print (paths[-1])

    #Now we have a list of paths, each path being a list of tuples (contig, orientation)
    k = 1 #starts with 1-mers
    go_on = True
    limit_abundance =0
    extra_paths = [] #correspond to contigs built in the previous iteration

    while go_on:

        graph_nodes = set()
        graph_edges = set() #edges are represented by a tuple of tuple ((kmer1, orientation1), (kmer2, orientation2))
        graph_edges_dict = dict() #edges are represented by an assiciation of tuple ((kmer1, orientation1) : (kmer2, orientation2))
        kmer_abundance = dict() #abundance of each kmer

        graph_nodes, graph_edges, graph_edges_dict, kmer_abundance = build_graph(k, paths+extra_paths, graph_nodes, graph_edges, graph_edges_dict, kmer_abundance)

        clean_graph(graph_nodes, graph_edges, graph_edges_dict, kmer_abundance, limit_abundance)

        #now build the contigs from the graph
        contigs = build_contigs(graph_nodes, graph_edges, graph_edges_dict, k)


        # print("\nHere are all the contigs : ", k)
        # for contig in contigs:
        #     print(contig)
        extra_paths = []

        for contig in contigs:
            for copy in range (limit_abundance+1) : #to be sure that they will be considered for the next iteration
                extra_paths.append(contig)
        output_graph(graph_nodes, graph_edges, "tmp/DBG_long_reads.gdf")


        k += 1
        new_graph_nodes = set()

        if k == 3 :
            go_on = False

        #output the graph in gdf format


    sys.exit()

    return segments

