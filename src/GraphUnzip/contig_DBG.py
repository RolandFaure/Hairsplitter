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
                # if kmer_foward[0] == ("edge_31@3@0",">") or kmer_reverse[0] == ("edge_31@3@0",">")  :
                #     print("for k = ", k,  " kmer : ", kmer, " and orientation : ", orientation)
                #     print("from path : ", path)
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

#function that builds contig from the DBG
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

#function that builds contigs from the graph and rename them
def build_final_graph(graph_nodes, graph_edges, graph_edges_dict, k):
    contigs = []
    already_visited = set()
    end_kmers_of_contig = {} #associates a contig with its end kmers
    end_kmers_to_contig = {} #associates an end kmer with its contig and an end

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


                # print("qdoicopiuc ", newkmer, " ", k)
                contig += [newkmer[k-1]]

                already_visited.add(new_node)
                node = new_node
                end = new_end
                        
            end_kmers_of_contig[tuple(contig)] = ((starting_node, not starting_end), (node, end))
            end_kmers_to_contig[(starting_node, not starting_end)] = (contig, 0)
            end_kmers_to_contig[(node, end)] = (contig, 1)

            contigs.append(contig)

    #add the edges
    edges = set()
    tuple_contigs = set()
    for contig in contigs:
        tuple_contigs.add(tuple(contig))
        end_kmers = end_kmers_of_contig[tuple(contig)]
        end = 0
        for end_kmer in end_kmers:
            for edge in graph_edges_dict[end_kmer]:
                new_edge = ((tuple(contig),end), (tuple(end_kmers_to_contig[edge][0]), end_kmers_to_contig[edge][1]))
                # print("New edge: ", new_edge)
                edges.add(new_edge)
            end += 1  

    output_graph(contigs, edges, "tmp/old_graph.gdf")

    #find the max chunk of each contig
    max_chunk_of_contig = {} #associates a contig with its max chunk
    for contig in contigs:
        for subcontig in contig :
            nb = int(subcontig[0].split("_")[-1])
            co = "_".join(subcontig[0].split("_")[:-1])
            if co not in max_chunk_of_contig :
                max_chunk_of_contig[co] = nb
            elif nb > max_chunk_of_contig[co] :
                max_chunk_of_contig[co] = nb

    #rename the contigs
    print("Compressssing the contigs...")
    old_contig_to_new_contig = {}
    newcontigs = []
    for contig in range(len(contigs)):
        # print("Conyigs sbefore compression : ", contigs[contig])
        newcontig = []
        last_edge = ""
        seen_last_chunk = False
        seen_first_chunk = False
        nb = 0
        
        
        for edge in range(len(contigs[contig])):
            name = "_".join(contigs[contig][edge][0].split("_")[:-1])
            nb = int(contigs[contig][edge][0].split("_")[-1])

            if name != last_edge :
                # if (seen_first_chunk == False or seen_last_chunk == False) and edge != 0:
                #     previous_name = "_".join(contigs[contig][edge-1][0].split("_")[:-1])
                #     print("Was not going to take ", last_edge, " anyway ", contigs[contig], " ", max_chunk_of_contig[previous_name], " ", int(contigs[contig][edge-1][0].split("_")[-1]))
                if edge != 0 :
                    seen_last_chunk = False
                    seen_first_chunk = False
                last_edge = name
                newcontig += [(name,contigs[contig][edge][1])]

            # if nb == max_chunk_of_contig[name] :
            #     seen_last_chunk = True

            # if nb == 0 :
            #     seen_first_chunk = True
            
            # if seen_first_chunk and seen_last_chunk: #both ends of a contig need to be seen
            #     newcontig += [(name,contigs[contig][edge][1])]

        newcontigs.append(newcontig)
        old_contig_to_new_contig[tuple(contigs[contig])] = tuple(newcontig)

    #rename the edges
    newedges = set()
    for edge in edges:
        newedges.add(((old_contig_to_new_contig[edge[0][0]], edge[0][1]), (old_contig_to_new_contig[edge[1][0]], edge[1][1])))
    
    output_graph(newcontigs, newedges, "tmp/new_graph.gdf")

    return newcontigs, newedges


#Master function of the file
#Input : initial gfa (as a list of segments), a GAF file with long reads mapped to the segments, names (which is an index numbering the contig)
#Output : new gfa (as a list of segments) corrected with long reads, and modified copiesnumber (taking into account contigs that have been duplicated)
def DBG_long_reads(segments, names, copiesnumber, gafFile):

    lines = []
    print("Reading the gaf file...")
    read_GAF(gafFile, 0, 0, lines)
    print("Finished going through the gaf file.")

    size_of_chunks = 1000

    paths = []
    for line in lines:
        #split the line on '>' and '<' to get the path
        contigs = re.split('[><]' , line)
        orientations = "".join(re.findall("[<>]", line))
        del contigs[0] #because the first element is always ''

        path = []
        for contig in range(len(contigs)):
            length_of_contig = segments[names[contigs[contig]]].get_length()
            #a contig is a series of chunks
            #we need to split the contig in chunks of size size_of_chunks
            if orientations[contig] == ">":
                for chunk in range(int(length_of_contig/size_of_chunks)+1):
                    path.append((contigs[contig] + "_" + str(chunk), orientations[contig]))
            else:
                for chunk in range(int(length_of_contig/size_of_chunks), -1, -1):
                    path.append((contigs[contig] + "_" + str(chunk), orientations[contig]))
                
        paths.append(path)
        # paths.append([(contigs[i], orientations[i]) for i in range(len(contigs))])

        # if contigs[0] == "edge_322@0@0" :
        #     print("contigs : ", contigs, " and orientations : ", orientations, " and line : ", line)
        #     print (paths[-1])

    #Now we have a list of paths, each path being a list of tuples (contig, orientation)
    k = 1 #starts with 1-mers
    go_on = True
    limit_abundance =0
    extra_paths = [] #correspond to contigs built in the previous iteration

    while go_on:

        print("GraphUnzip: unzipping with long reads, k=", k)


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

        if k == 10 :
            go_on = False

        #output the graph in gdf format

    #create the final assembly graph
    contigs, edges = build_final_graph(graph_nodes, graph_edges, graph_edges_dict, k-1)

    #create the contigs first
    new_segments = []
    end_of_contigs = {}
    for contig in contigs:
        sc = 0
        previous_end = 0
        previous_seg = sg.Segment([], [], [], [])
        end = 0
        for subcontig in contig :

            s = segments[names[subcontig[0]]]
            # print("sffdsq : ", s.names, " ", s.insideCIGARs)
            newseg = sg.Segment(s.names, s.orientations, s.lengths , segInsideCIGARs = s.insideCIGARs)

            if subcontig[1] == ">" :
                end = 0
            else :
                end = 1

            if sc > 0 : #link to previous subcontigs
                #find the CIGAR
                index = sg.find_this_link(previous_seg, previous_end, s.links[end], s.otherEndOfLinks[end] , warning=True)
                cigar = s.CIGARs[end][index]

                sg.add_link(newseg, end, new_segments[-1], previous_end, CIGAR = cigar)

            new_segments.append(newseg)
            previous_end = 1 - end
            previous_seg = s

            if sc == 0 :
                end_of_contigs[tuple(contig)] = [(len(new_segments)-1, end)]
            if sc == len(contig) -1 :
                end_of_contigs[tuple(contig)] += [(len(new_segments)-1, previous_end)]

            sc += 1

    #now add all the edges in the final graph
    for edge in edges :
        #compute the length of the overlap
        contig1 = edge[0][0]
        if edge[0][1] == 0:
            contig1 = reverse_complement(contig1)
        contig2 = edge[1][0]
        if edge[1][1] == 1 :
            contig2 = reverse_complement(contig2)

        length_of_overlap = 0
        for number_of_contigs in range(min(len(contig1), len(contig2)), 0, -1):
            length_of_overlap = 0
            good_overlap = False
            for c in range(number_of_contigs):
                if contig1[len(contig1) - number_of_contigs + c] != contig2[c] :
                    good_overlap = False
                    break
                length_of_overlap += segments[names[contig2[c][0]]].length

            if good_overlap :
                break
            
            number_of_contigs += 1
        
        # if length_of_overlap == 0 or True :
            # print("edge qdsc: ", edge, " ", contig1, " ", contig2, " ", length_of_overlap)
            # sys.exit()
        
        sg.add_link(new_segments[end_of_contigs[edge[0][0]][edge[0][1]][0]], end_of_contigs[edge[0][0]][edge[0][1]][1], \
                    new_segments[end_of_contigs[edge[1][0]][edge[1][1]][0]], end_of_contigs[edge[1][0]][edge[1][1]][1],
                    CIGAR = str(length_of_overlap)+'M')


    return new_segments

