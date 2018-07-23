#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
res = SkewArray(Genome=text)
print('SkewArray: {}'.format(" ".join(str(x) for x in res)))

res2 = MinimumSkew(Genome=text)
print(res2)