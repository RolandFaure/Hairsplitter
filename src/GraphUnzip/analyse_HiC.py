import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from transform_gfa import load_gfa
import basic_functions as bf


def short_distance_interactions(fragcontacts, fraglist):

    scores = []
    scoreFar = []
    for i in fragcontacts:
        if i[0] + 1 == i[1]:
            scores += [i[2]]
        elif i[1] - i[0] > 10:
            scoreFar += [i[2]]
    plt.hist(scores, bins=2000)
    plt.xlim([0, 70])
    plt.ylim([0, 1000])
    plt.show()

    plt.hist(scoreFar, bins=2000)
    plt.xlim([0, 70])
    plt.ylim([0, 1000])
    plt.show()

def distance_law(hiccontactsfile, fragmentList, header = True):
    tableDistance = (
        []
    )  # we're going to visualize distance law by looking at the inside of contigs
    tableIntensity = []
    
    f = open(hiccontactsfile,'r')
    
    for line in f:
        if not header :
            l = line.strip('\n').split('\t')
            i = [int(l[0]), int(l[1]), int(l[2])]
            if (fragmentList[i[0]][0] == fragmentList[i[1]][0] and i[2] > 0 ):  # i.e. we are in the same contig
    
                distance = fragmentList[i[1]][1] - fragmentList[i[0]][1]
                tableDistance += [distance]
                tableIntensity += [i[2]] #division to normalize to the probability of interaction per bp
                #/ fragmentList[i[0]][3] / fragmentList[i[1]][3]
        else :
            header = False
            
    plt.scatter(tableDistance, tableIntensity, alpha=0.2)
    plt.xlabel("Distance (bp)")
    plt.ylabel("Intensity of contact")
    plt.xlim([0, 200000])
    plt.ylim([0, 25])
    plt.show()

#    bf.export_to_csv([tableDistance, tableIntensity], 'listsPython/distanceIntensite.csv')

def with_how_many_contig_does_one_contig_interact(hiccontactsfile, fragmentList):

    with open(hiccontactsfile) as f:

        what_contigs_interact_with_this_one = [
            [] for i in range(fragmentList[-1][0] + 1)
        ]
        for line in f:

            line = line.strip("\n")
            line = line.split("\t")

            if line[0] != "487796":  # because the first line is a header
                contact = [int(line[0]), int(line[1]), int(line[2])]

                if (
                    contact[2] > 4
                ):  # this diminishes a lot the number of contig one contig interacts with, probably filtering out errors and rare events

                    contig1 = fragmentList[contact[0]][0]
                    contig2 = fragmentList[contact[1]][0]

                    if contig2 not in what_contigs_interact_with_this_one[contig1]:
                        what_contigs_interact_with_this_one[contig1] += [contig2]
                    if contig1 not in what_contigs_interact_with_this_one[contig2]:
                        what_contigs_interact_with_this_one[contig2] += [contig1]

        how_many_contigs_interact_with_this_one = [
            len(x) for x in what_contigs_interact_with_this_one
        ]
        plt.hist(how_many_contigs_interact_with_this_one, bins=20)
        # plt.xlim([0,200])
        plt.xlabel(
            "Number of other contigs interacting with one at least by 4 contacts"
        )
        plt.ylabel("Number of contig interacting with x others")
        plt.show()

# breaking down contigs to see how much HiC contact have fragments that are actually touching
def testHiC_vs_GFA(hiccontacts, info_contigs):

    contactNumber = 0
    score = []

    for contig in range(len(info_contigs)):

        start_frag = info_contigs[contig][3]
        end_frag = info_contigs[contig][3] + info_contigs[contig][2]

        if info_contigs[contig][2] > 29:
            cut = random.randint(
                15, info_contigs[contig][2] - 14
            )  # we cut the contig in two random parts

            score += [0]

            if contig < 5:
                print(contig, start_frag, end_frag, cut)

            while hiccontacts[contactNumber][0] < start_frag + cut:

                if (
                    hiccontacts[contactNumber][0] > start_frag
                    and hiccontacts[contactNumber][1] >= start_frag + cut
                    and hiccontacts[contactNumber][1] < end_frag
                ):
                    score[-1] += hiccontacts[contactNumber][2]

                contactNumber += 1

        while (
            hiccontacts[contactNumber][0] < end_frag
            and contactNumber < len(hiccontacts) - 1
        ):
            contactNumber += 1

    print(score)
    plt.hist(score)


#hiccontacts = bf.read_abs_fragments_contact_weighted('data/results/abs_fragments_contacts_weighted.txt')
# hiccontacts = import_from_csv('listsPython/hiccontacts.csv')
# print(hiccontacts[:20])
# print(hiccontacts[:100])
#fragmentList = bf.read_fragment_list("data/results/fragments_list.txt")
# print(fragmentList[:100])
# infcontigs = read_info_contig('data/results/info_contigs.txt')
# links = gfa_to_python(1312)
# short_distance_interactions(hiccontacts, fragmentList)

# links = import_from_csv('listsPython/links.csv')
# links = [[int(i) for i in j] for j in links]

# hiccontacts = import_from_csv('listsPython/hiccontacts.csv')

# confirmationOfLinks = HiC_vs_GFAtwo('data/results/abs_fragments_contacts_weighted.txt', links, fragmentList)

# distance_law(hiccontacts, fragmentList)
# testHiC_vs_GFA(hiccontacts, infcontigs)
# determine_HiC_coverage(hiccontacts, infcontigs, fragmentList)
# confirmationOfLinks = import_from_csv('listsPython/confirmationsDeslinks.csv')
# print(confirmationOfLinks[:19])

# check_links(links)
# coverage = determine_HiC_coverage(hiccontacts, infcontigs, fragmentList)

# coverage = bf.import_from_csv("listsPython/HiCcoverage.csv")
# coverage = [x[0] for x in coverage]

# conf, confweight = HiC_vs_GFAtwo('data/results/abs_fragments_contacts_weighted.txt', links, fragmentList, coverage)
# print(conf[:20],confweight[:20])
# print(links[:20])

# with_how_many_contig_does_one_contig_interact('data/results/abs_fragments_contacts_weighted.txt', fragmentList)

# interaction_Matrix = interactionMatrix(
#    "data/results/abs_fragments_contacts_weighted.txt", fragmentList, coverage
# )
# im = sp.lil_matrix(interaction_Matrix)
# pickle.dump(im, "listsPython/interactionMatrix.pickle")
# bf.export_to_csv(interaction_Matrix, "listsPython/interactionMatrix.csv") 
# print(interaction_Matrix[217][323], interaction_Matrix[217][359])

#distance_law('data/results/abs_fragments_contacts_weighted.txt', fragmentList)

#print("Finished")

