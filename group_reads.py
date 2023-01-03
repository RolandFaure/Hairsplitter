#!/usr/bin/env python3
# this can be used alone, but is also used in HairSplitter 

import argparse
import re
import pickle
import numpy as np
import os
import sys
import time

#write a function to parse command line arguments
def parse_args():
    parser=argparse.ArgumentParser(description="cluster reads")
    parser.add_argument("--alignment", "-a", required=False, default=None, help="all-vs-all mappings in sam format (or dumped as a pickle in previous run)(if not provided, will be computed from the reads with minimap2)")
    parser.add_argument("--reads", "-r", required=True, help="reads file")
    parser.add_argument("--outgraph", "-o", required=False, default="plasmid_graph.gdf", help="output graph file")
    parser.add_argument("--consensus_out", "-c", required=False, default="consensus.fasta", help="output consensus sequences file")
    return parser.parse_args()

#function to compure the distance between two reads from the sam record
def compute_distance(line):
    #parse the cigar string: M count as 0, I and D and S count as 1
    cigar=line[5]
    #split the cigar on the letters
    splitcigar = re.split("(\D)", cigar)
    #get rid of the empty strings in splitcigar
    splitcigar = [x for x in splitcigar if x != ""]
    
    distCIGAR = 0
    totalCIGAR = 0
    lastchar = ""
    for i in splitcigar:
        if i in ["I", "D", "S", "H", "P", "N", "X"]:
            lastint = int(lastchar)
            distCIGAR += lastint
            totalCIGAR += lastint
            lastchar = ""
        elif i == "M" or i == "=" :
            lastint = int(lastchar)
            totalCIGAR += lastint
            lastchar=""
        else:
            lastchar += i
    if totalCIGAR > 0:
        distSup = float(distCIGAR)/totalCIGAR
    else:
        distSup = 1

    #find the field that starts with de:f:
    for field in line:
        if field[:5]=="de:f:":
            dist=float(field.split(":")[2]) + distSup
            #if distance is 0, print the sam line
            if dist==0:
                print("distance between two reads is 0 : " + "\t".join([line[i] for i in range(3)]), " ", field)
            return dist

    return -1

#function to parse the sam file and return a matrix of distance between reads
def parse_sam(sam_file):

    #open the sam file
    sam=open(sam_file, 'r')

    #go through the sam file and inventory the reads
    read_dict={}
    read_names={}
    read_dict_size=0
    for line in sam:
        #skip the header
        if line[0]=="@":
            continue
        #split the line
        line=line.split()
        #get the read name
        read_name=line[0]
        #get the reference name
        ref_name=line[2]
        #add the read to the dictionary
        if read_name not in read_dict:
            read_dict[read_name]=read_dict_size
            read_names[read_dict_size]=read_name
            read_dict_size+=1
        #add the reference to the dictionary
        if ref_name not in read_dict:
            read_dict[ref_name]=read_dict_size
            read_names[read_dict_size]=ref_name
            read_dict_size+=1
    
    #create a matrix of distances between all the reads
    dist_matrix=[[1 for i in range(read_dict_size)] for j in range(read_dict_size)]
    #restart from the beginning of the sam file

    #count the total number of lines in the sam file
    sam.seek(0)
    total_lines= sum([1 for line in sam])
    sam.seek(0)

    line_number=0
    #go through the sam file and fill the matrix
    for line2 in sam:

        line_number += 1
        #print the progress if the line is a multiple of 10000
        if line_number%1000==0:
            print(int(line_number/total_lines*100), "%", end="\r")
        #skip the header
        if line2[0]=="@":
            continue
        #split the line
        line=line2.split()
        #get the read name
        read_name=line[0]
        #get the reference name
        ref_name=line[2]
        #get the distance in the de:f: field
        #try to get the distance if the read is not mapping to itself
        if read_name!=ref_name and ref_name!="*":
            dist = 1
            if line[1] != "0":
               dist = compute_distance(line)
            #if the distance is not found, leave it to 1
            else :
                print("unmatched: "+line[0])
            #add the distance to the matrix
            if dist < dist_matrix[read_dict[read_name]][read_dict[ref_name]] :
                dist_matrix[read_dict[read_name]][read_dict[ref_name]]=dist
                dist_matrix[read_dict[ref_name]][read_dict[read_name]]=dist

    #close the sam file
    sam.close()

    #export the matrix
    with open("tmp/dist_matrix.pickle", "wb") as f:
        pickle.dump((dist_matrix, read_names, read_dict), f)

    #return the matrix
    return dist_matrix, read_names, read_dict

#function to parse the paf and directly return a graph
def parse_paf(paf_file):

    #open the paf file
    paf=open(paf_file, 'r')

    #go through the paf file and inventory the reads
    read_dict={}
    #dictionary associate the index to the read name
    read_names={}

    read_dict_size=0
    for line in paf:
        #split the line
        line=line.split()
        #get the read name
        read_name=line[0]
        #get the reference name
        ref_name=line[5]
        #add the read to the dictionary
        if read_name not in read_dict:
            read_dict[read_name]=read_dict_size
            read_names[read_dict_size]=read_name
            read_dict_size+=1
        #add the reference to the dictionary
        if ref_name not in read_dict:
            read_dict[ref_name]=read_dict_size
            read_names[read_dict_size]=ref_name
            read_dict_size+=1
    
    #create the graph
    graph=[[0 for i in range(read_dict_size)] for j in range(read_dict_size)]

    #restart from the beginning of the paf file
    paf.seek(0)

    #go through the paf file and fill the graph
    for line in paf:
        #split the line
        line=line.split()
        #get the read name
        read_name=line[0]
        #get the reference name
        ref_name=line[5]
        #fill the matrix with a 1 if the read is mapped to the reference
        graph[read_dict[read_name]][read_dict[ref_name]]=1
        graph[read_dict[ref_name]][read_dict[read_name]]=1

    #close the paf file
    paf.close()

    #return the graph and the read dictionary
    return graph, read_names

#function to process the distance matrix and return a graph
def process_matrix(dist_matrix):

    #define the graph, which is the same size as the distance matrix
    graph=[[0 for i in range(len(dist_matrix))] for j in range(len(dist_matrix))]

    #go through the distance matrix and fill the graph
    #put 1 for the 5 best scores of each row
    best_scores = [1 for i in range(5)]
    for i in range(len(dist_matrix)):
        #get the 5 best scores
        best_scores=sorted(dist_matrix[i])[:5]
        #take out of best scores all the 0 and and all the 1
        best_scores=[x for x in best_scores if x!=1]
        #go through the scores and fill the graph
        for j in range(len(dist_matrix[i])):
            if dist_matrix[i][j] in best_scores:
            # if dist_matrix[i][j] < 0.2:
                graph[i][j]=1

    

    return graph


#function to write the graph to a file in CSV format
def write_graph(graph, edgenames, clusters, outgraph):

    #open the output file
    out=open(outgraph, 'w')

    out.write("nodedef>name VARCHAR,label VARCHAR, plasmid INT, cluster INT\n")
    for e in edgenames.keys():
        out.write(edgenames[e] + ","+str(e)+ "," + edgenames[e].split("_")[0]+","+str(clusters[e])+"\n")

    out.write("edgedef>node1 VARCHAR,node2 VARCHAR")

    #go through the graph and write it to the file
    line_number=0
    for i in range(len(graph)):
        for j in range(0, len(graph[i])):
            if graph[i][j]!=0:
                out.write(edgenames[i]+","+edgenames[j])
                #out.write(","+str(graph[i][j]))
                out.write("\n")
                line_number+=1
    print("number of edges: ", line_number)


#function to cluster the adjacency matrix using the Chinese Whispers algorithm
def cluster_graph_cw(graph):

    #convert the graph to a numpy array
    adjacency_matrix=np.array(graph)

    changed = True
    communities = [i for i in range(len(adjacency_matrix))]
    while changed:
        changed = False
        for i in range(len(adjacency_matrix)):
            #get the neighbors of the node
            neighbors = [j for j in range(len(adjacency_matrix[i])) if adjacency_matrix[i][j] != 0]
            #get the community of the neighbors
            neighbors_communities = [communities[j] for j in neighbors]
            #get the most frequent community
            if len(neighbors_communities)==0:
                continue
            most_frequent_community = max(set(neighbors_communities), key=neighbors_communities.count)
            #if the most frequent community is different from the current community, change it
            if most_frequent_community != communities[i]:
                communities[i] = most_frequent_community
                changed = True

    #print(communities)
    
    return communities

#function to compute the consensus sequence of a cluster
def consensus_sequence(readFile, reads_names, communities, outfile):

    #go through the read file and load all the reads in a dictionary
    reads = {}
    r = open(readFile, "r")
    read_name = ""
    next = False
    for line in r:
        if line[0]==">" or line[0]=="@":
            read_name = line[1:].split()[0]
            reads[read_name] = ""
            next = True
        elif read_name in reads and next:
            reads[read_name] += line.strip()
            next = False


    #create a dictionary associating a list of sequence to each community
    community_sequences = {}
    for i in range(len(communities)):
        if communities[i] not in community_sequences:
            community_sequences[communities[i]] = []
        if i in reads_names: #to avoid cases where the read is not in the sam file
            community_sequences[communities[i]].append(reads[reads_names[i]])

    #go through the communities and compute the consensus sequence using racon
    consensus_sequences = {}
    os.system("mkdir tmp/")
    for c in community_sequences:
        if len(community_sequences[c])>5:
            #create a fasta file with the sequences of the community
            fo = open("tmp/read0_"+str(c)+".fasta", "w")
            fo.write(">"+str(c)+"\n")
            fo.write(community_sequences[c][0]+"\n")
            fo.close()
            f = open("tmp/cluster_"+str(c)+".fasta", "w")
            for i in range(len(community_sequences[c])):
                f.write(">read_"+str(i)+"\n")
                f.write(community_sequences[c][i]+"\n")
            f.close()

            print("Creating consensus for cluster "+str(c)+"...")
            commandMinimap = "minimap2 -ax map-pb tmp/read0_"+str(c)+".fasta tmp/cluster_"+str(c)+".fasta > tmp/mapping_"+str(c)+".sam" + " 2> /tmp/minimap.log"
            os.system(commandMinimap)
            command = "racon -t 1 tmp/cluster_"+str(c)+".fasta tmp/mapping_"+str(c)+".sam tmp/read0_"+str(c)+".fasta > tmp/consensus_"+str(c)+".fasta 2> /tmp/racon.log"
            code = os.system(command)
            if code != 0:
                print("RACON ERROR while running command: "+command)
                sys.exit()
            print("Done polishing the consensus sequence of cluster "+str(c)+".")
    
    #concatenate all the consensus sequences
    os.system("cat tmp/consensus_*.fasta > tmp/consensuses.fasta")

    #now check if you can merge the consensus sequences
    #perform all-vs-all alignment of the consensus sequences
    command = "minimap2 -ax asm5 tmp/consensuses.fasta tmp/consensuses.fasta > tmp/ava_consensuses.sam 2> /tmp/minimap.log"
    os.system(command)
    #parse the sam file
    fs = open("tmp/ava_consensuses.sam", "r")
    
    #set of consensus sequences to discard
    to_discard = set()
    #go through the sam file to see if two consensus sequences can be merged
    for line2 in fs:
        #skip the header
        if line2[0]=="@":
            continue
        #split the line
        line=line2.split()
        #get the read name
        read_name=line[0]
        #get the reference name
        ref_name=line[2]
        #get the distance in the de:f: field
        if read_name!=ref_name and ref_name!="*":
            dist = 1
            if line[1] != "0":
                dist = compute_distance(line)
                #if the distance is less than 2% of the length of the consensus sequence, merge the two consensus sequences
                if dist < 0.02:
                    if read_name not in to_discard and ref_name not in to_discard:
                        to_discard.add(read_name)
                        print("Merging "+read_name+" and "+ref_name+"...")
    fs.close()

    #output the finished consensus sequences
    f=open(outfile, "w")
    #go thorugh tmp/consensuses.fasta and output the sequences that are not in to_discard
    fs = open("tmp/consensuses.fasta", "r")
    writenext = True
    for line in fs:
        if line[0]==">":
            read_name = line[1:].split()[0]
            if read_name not in to_discard:
                f.write(line)
                writenext = True
            else:
                writenext = False
        elif writenext:
            f.write(line)
    f.close()



#main function
def main():
    #get the command line arguments
    args=parse_args()

    alignment = args.alignment

    #mkdir tmp/
    os.system("mkdir tmp/")

    #check the alignment file format
    #if a distance matrix is provided, process it
    graph = []
    #if the alignment file is a fasta, run minimap2
    if args.alignment==None:
        #run minimap2
        command = "minimap2 -ax ava-ont "+args.reads+" "+args.reads+" > tmp/ava.sam 2> tmp/minimap.log"
        print("Fasta file provided: I run minimap2 to compute the all-vs-all alignment of the reads.")
        os.system(command)
        alignment = "tmp/ava.sam"
        print("Finished aligning. The alignment is in tmp/ava.sam. (if you need to rerun the program, you can use this file as input with the -a option)")

    if alignment[-7:]==".pickle":
        dist_matrix, edge_names, read_ids  = pickle.load(open(alignment, "rb"))
        print(len(dist_matrix), len(dist_matrix[0]))
        #process the distance matrix
        graph=process_matrix(dist_matrix)
    elif alignment[-3:]=="sam":
        #parse the sam file
        dist_matrix, edge_names, read_ids = parse_sam(alignment)
        #process the distance matrix
        graph=process_matrix(dist_matrix)
        

    #cluster the graph
    communities = cluster_graph_cw(graph)
    #write the graph to a file
    write_graph(graph, edge_names, communities, args.outgraph)

    #take out all the "*" in the edge names
    edge_names = {k:v for k,v in edge_names.items() if v!="*"}
        
    #compute the consensus sequence of each cluster
    consensus_sequence(args.reads, edge_names, communities, args.consensus_out)

    #remove the tmp/ directory
    os.system("rm -r tmp/")


#run the main function from command line
if __name__=="__main__":
    main()



