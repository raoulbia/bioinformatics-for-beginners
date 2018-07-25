#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
skew arrays are used for finding DnaA boxes in the region of the E. coli genome
hypothesized by the minimum skew as ori.
see here for an example: https://stepik.org/lesson/23062/step/1?unit=6794
Note: hidden messages have a tendency to cluster within a genome, and most of them have nothing to do with replication.
"""
def SkewArray(Genome):

    Genome = '_' + Genome
    skewArray = [0 for i in Genome]

    n = len(Genome)
    for i in range(n):

        if Genome[i] in ('A', 'T'):
            skewArray[i] = skewArray[i-1]
        elif Genome[i] in ('G'):
            skewArray[i] = skewArray[i-1]+1
        elif Genome[i] in ('C'):
            skewArray[i] = skewArray[i-1]-1

    return skewArray

def MinimumSkew(Genome):

    positions = []
    skewArray = SkewArray(Genome=Genome)

    # find the minimum value of all values in the skew array
    min_val = min(skewArray)

    # range over the length of the skew array and add all positions achieving the min to positions
    for i in range(len(skewArray)):
        if skewArray[i] == min_val:
            positions.append(i)
    return positions

text = 'CATGGGCATCGGCCATACGCC'
text = 'AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT'
text = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
text = 'CATTCCAGTACTTCGATGATGGCGTGAAGA'
res = SkewArray(Genome=text)
print('SkewArray: {}'.format(" ".join(str(x) for x in res)))

res2 = MinimumSkew(Genome=text)
print(res2)